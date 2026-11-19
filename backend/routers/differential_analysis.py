from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from pathlib import Path
from typing import List
import asyncio
import tempfile
import pandas as pd
import numpy as np
import os
from services.file_service import read_file, save_sequences
from services.constants import ColumnName
import state

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))
EDGER_SCRIPT = Path(__file__).parent.parent / "scripts" / "edger_qlf.R"

class EdgeRPairTestInput(BaseModel):
    cond1_paths: List[str]  # Condition 1 FASTA file paths
    cond2_paths: List[str]  # Condition 2 FASTA file paths
    lfc_cutoff: float = 1.0  # |log2FC| threshold for significance
    logcpm_cutoff: float = 0.0  # minimum logCPM threshold for significance (0 = off)
    pvalue_cutoff: float = 1.0  # adjusted p-value threshold (1.0 = disabled)
    output_format: str = "csv"

class DifferentialExpressionResponse(BaseModel):
    status: str
    result: str

@router.post("/differential-expression", response_model=DifferentialExpressionResponse)
async def differential_expression(params: EdgeRPairTestInput):
    """
    Perform differential expression analysis between two conditions.
    Uses edgeR (scripts/edger_qlf.R): filterByExpr, TMM normalization, quasi-likelihood
    F-test with Benjamini-Hochberg FDR correction (tested only with >= 2 replicates per condition).
    
    Args:
        params: Input parameters with condition file paths and p-value cutoff
        
    Returns:
        Results with log fold changes and adjusted p-values
    """
    # Validate inputs
    if not params.cond1_paths or not params.cond2_paths:
        raise HTTPException(status_code=400, detail="Both conditions must have at least one file")

    if len(params.cond1_paths) > 26 or len(params.cond2_paths) > 26:
        raise HTTPException(status_code=400, detail="Maximum 26 replicates per condition")

   # -----------------------
    # LOAD CONDITION 1
    # -----------------------
    fadf_cond1 = []
    original_lib_sizes_cond1 = []

    for i, file_path in enumerate(params.cond1_paths):
        df = read_file(UPLOAD_DIR / file_path)

        if ColumnName.SEQUENCES not in df or ColumnName.READS not in df:
            raise HTTPException(status_code=400, detail=f"{file_path} missing required columns")

        original_lib_sizes_cond1.append(df[ColumnName.READS].sum())

        suffix = chr(ord("A") + i)
        fadf_cond1.append(
            df[[ColumnName.SEQUENCES, ColumnName.READS]]
            .rename(columns={ColumnName.READS: f"Reads.{suffix}"})
        )

    # -----------------------
    # LOAD CONDITION 2
    # -----------------------
    fadf_cond2 = []
    original_lib_sizes_cond2 = []

    for i, file_path in enumerate(params.cond2_paths):
        df = read_file(UPLOAD_DIR / file_path)

        if ColumnName.SEQUENCES not in df or ColumnName.READS not in df:
            raise HTTPException(status_code=400, detail=f"{file_path} missing required columns")

        original_lib_sizes_cond2.append(df[ColumnName.READS].sum())

        suffix = chr(ord("a") + i)
        fadf_cond2.append(
            df[[ColumnName.SEQUENCES, ColumnName.READS]]
            .rename(columns={ColumnName.READS: f"Reads.{suffix}"})
        )

    # -----------------------
    # MERGE ALL FILES
    # -----------------------
    # Outer join: a sequence missing from a file has 0 reads there (edgeR models zero counts)
    cond_merged = fadf_cond1[0]
    for df in fadf_cond1[1:] + fadf_cond2:
        cond_merged = cond_merged.merge(df, on=ColumnName.SEQUENCES, how="outer")

    # -----------------------
    # BUILD COUNT MATRIX
    # -----------------------
    cond1_cols = [c for c in cond_merged.columns if c.startswith("Reads.") and c[6].isupper()]
    cond2_cols = [c for c in cond_merged.columns if c.startswith("Reads.") and c[6].islower()]

    count_matrix = cond_merged[cond1_cols + cond2_cols].fillna(0).to_numpy(dtype=np.int32)
    sequences = cond_merged[ColumnName.SEQUENCES].to_numpy()
    del cond_merged

    # filterByExpr always drops sequences with < 15 total reads; dropping them here only
    # shrinks what is passed to R (library sizes are passed separately, so results are unchanged)
    sequences_before_filter = len(sequences)
    prefilter = count_matrix.sum(axis=1) >= 15
    count_matrix = count_matrix[prefilter]
    sequences = sequences[prefilter]

    if len(sequences) < 10:
        raise HTTPException(
            status_code=400,
            detail=f"Too few sequences with >= 15 total reads (found {len(sequences)} of "
                   f"{sequences_before_filter} unique sequences). Your FASTA files may have low read counts."
        )

    lib_sizes = original_lib_sizes_cond1 + original_lib_sizes_cond2
    n_cond1, n_cond2 = len(cond1_cols), len(cond2_cols)

    # -----------------------
    # EDGER (filterByExpr, TMM, quasi-likelihood F-test, BH)
    # -----------------------
    with tempfile.TemporaryDirectory() as tmp:
        counts_file = Path(tmp) / "counts.bin"
        out_file = Path(tmp) / "edger.bin"
        np.ascontiguousarray(count_matrix.T).tofile(counts_file)  # column-major for R
        try:
            proc = await asyncio.create_subprocess_exec(
                "Rscript", str(EDGER_SCRIPT), str(counts_file), str(len(count_matrix)),
                str(n_cond1), str(n_cond2), ",".join(str(int(s)) for s in lib_sizes), str(out_file),
                stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE,
            )
        except FileNotFoundError:
            raise HTTPException(status_code=500, detail="Rscript not found: R and edgeR must be installed to run differential analysis")
        # Register so /cancel (or the next request) can kill R
        state.current_process = proc
        try:
            _, stderr = await proc.communicate()
        finally:
            if state.current_process is proc:
                state.current_process = None
        if proc.returncode < 0:
            raise HTTPException(status_code=500, detail="edgeR was cancelled")
        if proc.returncode != 0:
            raise HTTPException(status_code=500, detail=f"edgeR failed: {stderr.decode(errors='replace')[-1000:]}")
        edger_out = np.fromfile(out_file, dtype=np.float64).reshape(4, len(count_matrix))

    kept = edger_out[0] == 1
    log_fc, log_cpm, fdr = edger_out[1][kept], edger_out[2][kept], edger_out[3][kept]
    sequences = sequences[kept]

    if len(sequences) < 10:
        raise HTTPException(
            status_code=400,
            detail=f"Too few sequences passed edgeR's filterByExpr (found {len(sequences)} of "
                   f"{sequences_before_filter} unique sequences). Your FASTA files may have low read counts."
        )

    # -----------------------
    # CLASSIFY
    # -----------------------
    lfc_sig = (np.abs(log_fc) >= params.lfc_cutoff) & (log_cpm >= params.logcpm_cutoff)

    if min(n_cond1, n_cond2) < 2:
        # No test without replicates; classify by logFC/logCPM thresholds only
        p_col = np.full(len(sequences), np.nan)
        p_class = np.where(lfc_sig, "Sig.", "Not Sig.")
    else:
        # PValue column holds the BH-adjusted p-value (FDR)
        p_col = np.round(fdr, 6)
        pvalue_sig = p_col <= params.pvalue_cutoff
        p_class = np.where(lfc_sig & pvalue_sig, "Sig.", "Not Sig.")

    # -----------------------
    # RESULTS TABLE
    # -----------------------
    d_test = pd.DataFrame({
        "logFC": np.round(log_fc, 3),
        "logCPM": np.round(log_cpm, 3),
        "PValue": p_col,
        "PClass": p_class,
        "Sequence": sequences,
    }).sort_values("logFC", key=np.abs, ascending=False)

    # -----------------------
    # SAVE OUTPUT
    # -----------------------
    base = Path(params.cond1_paths[0]).stem
    output = UPLOAD_DIR / f"{base}_diff_expr.{params.output_format}"
    save_sequences(d_test, output, params.output_format)

    return DifferentialExpressionResponse(
        status="ok",
        result=output.name,
    )
