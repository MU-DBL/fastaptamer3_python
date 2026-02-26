from services.file_service import read_file, save_sequences
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from pathlib import Path
from typing import List, Optional, Dict
import pandas as pd
import re
import os
from services.constants import ColumnName

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))

class MotifTrackerInput(BaseModel):
    input_paths: List[str]  # List of FASTA file paths
    population_names: List[str]
    query_list: List[str]  # List of motif patterns
    query_aliases: Optional[List[str]] = None
    motif_type: str = "Nucleotide"
    output_format: str = "csv"


IUPAC_MAP = {
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]",
    "D": "[AGT]", "H": "[ACT]", "V": "[ACG]",
    "N": "[ACGT]",
}

def format_motif(motif: str, motif_type: str) -> str:
    motif = motif.upper()
    if motif_type == "Nucleotide":
        for k, v in IUPAC_MAP.items():
            motif = motif.replace(k, v)
    return motif

@router.post("/motif-tracker")
async def motif_tracker(params: MotifTrackerInput):
    validate_inputs(params)

    dfs, total_reads = load_fasta_files(params)

    motifs = [format_motif(m, params.motif_type) for m in params.query_list]
    pattern_map = dict(zip(params.query_list, motifs))

    rows = []

    for pop_idx, df in enumerate(dfs, start=1):
        for motif, regex in pattern_map.items():
            matches = df[df[ColumnName.SEQUENCES].str.contains(regex, regex=True, na=False)]
            if matches.empty:
                continue

            matches = matches.assign(
                **{
                    ColumnName.MOTIF: motif,
                    ColumnName.POPULATION_NUMBER: pop_idx,
                    ColumnName.POPULATION_NAME: params.population_names[pop_idx - 1],
                }
            )
            rows.append(matches)

    if rows:
        result_df = (
            pd.concat(rows)
            .groupby(
                [ColumnName.POPULATION_NUMBER, ColumnName.POPULATION_NAME, ColumnName.MOTIF],
                as_index=False,
            )[ColumnName.READS]
            .sum()
            .rename(columns={ColumnName.READS: ColumnName.TOTAL_READS})
        )

        result_df[ColumnName.PERCENTAGE] = result_df.apply(
            lambda r: round(
                r[ColumnName.TOTAL_READS] / total_reads[r[ColumnName.POPULATION_NUMBER] - 1] * 100,
                2,
            ) if total_reads[r[ColumnName.POPULATION_NUMBER] - 1] > 0 else 0.0,
            axis=1,
        )
    else:
        result_df = pd.DataFrame(
            columns=[
                ColumnName.POPULATION_NUMBER,
                ColumnName.POPULATION_NAME,
                ColumnName.MOTIF,
                ColumnName.TOTAL_READS,
                ColumnName.PERCENTAGE,
            ]
        )

    result_df = add_aliases(result_df, params, ColumnName.MOTIF)

    base = Path(params.input_paths[0]).stem
    output = UPLOAD_DIR / f"{base}_motif_tracker.csv"
    enrichment = UPLOAD_DIR / f"{base}_motif_tracker_enrichment.csv"

    save_sequences(result_df, output, "csv")
    tracker_enrichment(result_df, enrichment, "csv")

    return {
        "status": "ok",
        "result": output.name,
        "enrichment": enrichment.name,
    }


@router.post("/sequence-tracker")
async def sequence_tracker(params: MotifTrackerInput):
    validate_inputs(params)

    dfs, _ = load_fasta_files(params)

    rows = []

    for pop_idx, df in enumerate(dfs, start=1):
        # Filter sequences that exactly match the query list
        matches = df[df[ColumnName.SEQUENCES].isin(params.query_list)]
        if matches.empty:
            continue

        matches = matches.assign(
            **{
                ColumnName.POPULATION_NUMBER: pop_idx,
                ColumnName.POPULATION_NAME: params.population_names[pop_idx - 1],
            }
        )
        rows.append(matches)

    if rows:
        # Concatenate all matches - do NOT aggregate, keep individual sequences
        result_df = pd.concat(rows)
        # Select only relevant columns and arrange by Sequences
        result_df = result_df[[
            ColumnName.POPULATION_NUMBER,
            ColumnName.POPULATION_NAME,
            ColumnName.SEQUENCES,
            ColumnName.RANK,
            ColumnName.READS,
            ColumnName.RPU
        ]].sort_values(ColumnName.SEQUENCES)
    else:
        result_df = pd.DataFrame(
            columns=[
                ColumnName.POPULATION_NUMBER,
                ColumnName.POPULATION_NAME,
                ColumnName.SEQUENCES,
                ColumnName.RANK,
                ColumnName.READS,
                ColumnName.RPU,
            ]
        )

    result_df = add_aliases(result_df, params, ColumnName.SEQUENCES)

    base = Path(params.input_paths[0]).stem
    output = UPLOAD_DIR / f"{base}_sequence_tracker.csv"
    enrichment = UPLOAD_DIR / f"{base}_sequence_tracker_enrichment.csv"

    save_sequences(result_df, output, "csv")
    tracker_enrichment(result_df, enrichment, "csv")

    return {
        "status": "ok",
        "result": output.name,
        "enrichment": enrichment.name,
    }
    

def tracker_enrichment(df: pd.DataFrame, output_path: Path, output_format: str):
    if df.empty:
        save_sequences(pd.DataFrame(), output_path, output_format)
        return

    query_col = (
        ColumnName.MOTIF if ColumnName.MOTIF in df.columns else ColumnName.SEQUENCES
    )
    value_col = (
        ColumnName.PERCENTAGE if ColumnName.PERCENTAGE in df.columns else ColumnName.RPU
    )

    records = []

    for query, sub in df.groupby(query_col):
        sub = sub.sort_values(ColumnName.POPULATION_NUMBER)

        for i in range(1, len(sub)):
            prev, curr = sub.iloc[i - 1], sub.iloc[i]
            
            # Calculate enrichment with proper handling of edge cases
            if prev[value_col] > 0:
                enrichment = round(curr[value_col] / prev[value_col], 2)
            elif curr[value_col] > 0:
                # If previous is 0 but current is not, show as "inf" or very high number
                enrichment = "inf"
            else:
                # Both are 0, no enrichment
                enrichment = 0.0

            records.append(
                {
                    ColumnName.COMPARISON: f"{curr[ColumnName.POPULATION_NAME]} : {prev[ColumnName.POPULATION_NAME]}",
                    ColumnName.QUERY: query,
                    ColumnName.ALIAS: curr.get(ColumnName.ALIAS, "-"),
                    ColumnName.ENRICHMENT: enrichment,
                }
            )

    save_sequences(pd.DataFrame(records), output_path, output_format)
    

def validate_inputs(params: MotifTrackerInput):
    if not params.input_paths:
        raise HTTPException(400, "input_paths is required")

    if len(params.input_paths) != len(params.population_names):
        raise HTTPException(400, "input_paths must match population_names")

    if params.query_aliases and len(params.query_aliases) != len(params.query_list):
        raise HTTPException(400, "query_aliases must match query_list")

def add_aliases(df, params, key_col):
    if not params.query_aliases or df.empty:
        return df

    alias_df = pd.DataFrame(
        {key_col: params.query_list, ColumnName.ALIAS: params.query_aliases}
    )
    return df.merge(alias_df, on=key_col, how="left")

def load_fasta_files(params: MotifTrackerInput):
    dfs = []
    total_reads = []

    for path in params.input_paths:
        filepath = UPLOAD_DIR / path
        if not filepath.exists():
            raise HTTPException(404, f"File not found: {path}")

        df = read_file(filepath)
        dfs.append(df)
        total_reads.append(df[ColumnName.READS].sum())

    return dfs, total_reads