from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field
from pathlib import Path
from typing import List
from collections import Counter
import os
import numpy as np
import pandas as pd
from scipy.stats import norm

from services.file_service import read_file, save_sequences
from services.constants import ColumnName

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))

DNA_ALPHABET = "ACGT"
PROTEIN_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

class MotifDiscoveryInput(BaseModel):
    input_path: str
    min_reads: int = Field(default=10, ge=1)
    length_range: List[int] = Field(default_factory=lambda: [5, 10])
    output_format: str = "csv"
    alphabet: str = "dna"  # "dna" or "protein"


def calculate_subsequence_frequencies(sequences: List[str]) -> Counter:
    freq = Counter()
    for seq in sequences:
        n = len(seq)
        for i in range(n):
            for j in range(i + 1, n + 1):
                freq[seq[i:j]] += 1
    return freq

def calculate_base_frequencies(sequences: List[str], alphabet: str) -> dict:
    counts = Counter("".join(sequences))
    total = sum(counts[b] for b in alphabet)
    uniform = 1.0 / len(alphabet)

    return {
        b: counts.get(b, 0) / total if total else uniform
        for b in alphabet
    }

def expected_motif_probability(motif: str, base_freqs: dict) -> float:
    prob = 1.0
    for b in motif:
        prob *= base_freqs.get(b, 0.25)
    return prob


def discover_enriched_motifs(
    sequences: List[str],
    subseq_freq: Counter,
    base_freqs: dict,
    lmin: int,
    lmax: int
) -> pd.DataFrame:

    results = []
    seq_lengths = [len(s) for s in sequences]

    for motif, obs in subseq_freq.items():
        L = len(motif)
        if not (lmin <= L <= lmax):
            continue

        p = expected_motif_probability(motif, base_freqs)
        total_positions = sum(max(0, n - L + 1) for n in seq_lengths)
        exp = p * total_positions

        if exp <= 0:
            continue

        R = obs / exp
        var = exp * (1 - p)
        Z = (obs - exp) / np.sqrt(var) if var > 0 else 0.0
        ZZ = Z / np.sqrt(L)
        P = 1 - norm.cdf(Z)

        results.append({
            ColumnName.MOTIF: motif,
            ColumnName.RATIO: R,
            ColumnName.P_VALUE: P,
            ColumnName.Z_SCORE: Z,
            ColumnName.ZZ_SCORE: ZZ,
            ColumnName.MOTIF_LENGTH: L,
            ColumnName.OBSERVED_COUNT: obs,
            ColumnName.EXPECTED_COUNT: exp,
        })

    if not results:
        return pd.DataFrame()

    df = pd.DataFrame(results)
    df = df.sort_values(ColumnName.ZZ_SCORE, ascending=False).reset_index(drop=True)
    df[ColumnName.RANK] = np.arange(1, len(df) + 1)

    return df


@router.post("/motif-discovery")
async def motif_discovery(params: MotifDiscoveryInput):

    if len(params.length_range) != 2:
        raise HTTPException(400, "length_range must be [min, max]")

    lmin, lmax = params.length_range
    if lmin > lmax:
        raise HTTPException(400, "length_range min must <= max")

    filepath = UPLOAD_DIR / params.input_path
    if not filepath.exists():
        raise HTTPException(404, f"File not found: {params.input_path}")

    try:
        df = read_file(filepath)
        if df.empty:
            raise HTTPException(400, "FASTA file is empty")

        df = df[df[ColumnName.READS] >= params.min_reads]
        if df.empty:
            raise HTTPException(
                400, f"No sequences with reads >= {params.min_reads}"
            )

        sequences = df[ColumnName.SEQUENCES].tolist()

        alphabet = PROTEIN_ALPHABET if params.alphabet == "protein" else DNA_ALPHABET
        subseq_freq = calculate_subsequence_frequencies(sequences)
        base_freqs = calculate_base_frequencies(sequences, alphabet)

        motif_df = discover_enriched_motifs(
            sequences,
            subseq_freq,
            base_freqs,
            lmin,
            lmax,
        )

        if motif_df.empty:
            output_df = pd.DataFrame(columns=[
                ColumnName.MOTIF,
                ColumnName.P_VALUE,
                ColumnName.ZZ_SCORE,
                ColumnName.MOTIF_LENGTH,
                ColumnName.RANK,
                ColumnName.SEQ_COUNT,
            ])
        else:
            motif_df[ColumnName.SEQ_COUNT] = motif_df[ColumnName.MOTIF].apply(
                lambda m: sum(m in s for s in sequences)
            )

            output_df = motif_df[[
                ColumnName.MOTIF,
                ColumnName.P_VALUE,
                ColumnName.ZZ_SCORE,
                ColumnName.MOTIF_LENGTH,
                ColumnName.RANK,
                ColumnName.SEQ_COUNT,
            ]].round({
                ColumnName.P_VALUE: 5,
                ColumnName.ZZ_SCORE: 3,
            })

        base = Path(params.input_path).stem
        output_path = UPLOAD_DIR / f"{base}_motif_discovery.{params.output_format}"

        save_sequences(output_df, str(output_path), params.output_format)

        return {
            "status": "ok",
            "result": output_path.name,
        }

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(
            500, f"Motif discovery failed ({type(e).__name__}): {str(e)}"
        )
