from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from pathlib import Path
from typing import List
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import os
from services.file_service import read_file, save_sequences
from services.constants import ColumnName

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))

class EdgeRPairTestInput(BaseModel):
    cond1_paths: List[str]  # Condition 1 FASTA file paths
    cond2_paths: List[str]  # Condition 2 FASTA file paths
    lfc_cutoff: float = 1.0  # |log2FC| threshold for significance
    output_format: str = "csv"

class DifferentialExpressionResponse(BaseModel):
    status: str
    result: str

@router.post("/differential-expression", response_model=DifferentialExpressionResponse)
async def differential_expression(params: EdgeRPairTestInput):
    """
    Perform differential expression analysis between two conditions.
    Uses Mann-Whitney U test (or t-test) with Benjamini-Hochberg FDR correction.
    
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
    cond1_merged = fadf_cond1[0]
    for df in fadf_cond1[1:]:
        cond1_merged = cond1_merged.merge(df, on=ColumnName.SEQUENCES, how="inner")

    cond2_merged = fadf_cond2[0]
    for df in fadf_cond2[1:]:
        cond2_merged = cond2_merged.merge(df, on=ColumnName.SEQUENCES, how="inner")

    cond_merged = cond1_merged.merge(cond2_merged, on=ColumnName.SEQUENCES, how="inner")

    if len(cond_merged) < 10:
        raise HTTPException(
            status_code=400, 
            detail=f"Too few sequences present in ALL files for DE analysis (found {len(cond_merged)}). "
                   "edgeR requires sequences to be present in all replicates. "
                   "This is normal if your replicates have different sequence compositions."
        )

    # -----------------------
    # BUILD COUNT MATRIX
    # -----------------------
    cond1_cols = [c for c in cond_merged.columns if c.startswith("Reads.") and c[6].isupper()]
    cond2_cols = [c for c in cond_merged.columns if c.startswith("Reads.") and c[6].islower()]

    count_matrix = cond_merged[cond1_cols + cond2_cols].to_numpy(dtype=float)
    sequences = cond_merged[ColumnName.SEQUENCES].to_numpy()

    # -----------------------
    # CPM NORMALIZATION (SAFE)
    # -----------------------
    lib_sizes = np.asarray(original_lib_sizes_cond1 + original_lib_sizes_cond2, dtype=float)

    cpm_matrix = np.divide(
        count_matrix,
        lib_sizes,
        out=np.zeros_like(count_matrix),
        where=lib_sizes != 0,
    ) * 1e6

    # -----------------------
    # FILTER LOW-COUNT SEQUENCES (CRITICAL)
    # -----------------------
    min_count = 5
    min_samples = 2
    
    sequences_before_filter = len(sequences)

    keep = (count_matrix >= min_count).sum(axis=1) >= min_samples

    count_matrix = count_matrix[keep]
    cpm_matrix = cpm_matrix[keep]
    sequences = sequences[keep]

    # NEW: Check if filtering left any sequences
    if len(count_matrix) == 0:
        raise HTTPException(
            status_code=400, 
            detail=f"No sequences passed filtering criteria (min_count={min_count} reads in at least {min_samples} samples). "
                f"Your files contained {sequences_before_filter} unique sequences before filtering. "
                "This suggests your FASTA files have low read counts. Please check your input files."
        )
    
    # -----------------------
    # LOGFC & LOGCPM
    # -----------------------
    n_cond1 = len(cond1_cols)
    pseudocount = 1.0

    cond1_cpm = cpm_matrix[:, :n_cond1]
    cond2_cpm = cpm_matrix[:, n_cond1:]

    mean1 = cond1_cpm.mean(axis=1) + pseudocount
    mean2 = cond2_cpm.mean(axis=1) + pseudocount

    log_fc = np.log2(mean2 / mean1)
    # logCPM: average log2 CPM across ALL samples (matching edgeR behavior)
    log_cpm = np.log2(cpm_matrix.mean(axis=1) + pseudocount)

    # -----------------------
    # STATISTICAL TESTING (CPM-NORMALIZED)
    # -----------------------
    insufficient_replicates = n_cond1 < 3 or len(params.cond2_paths) < 3

    lfc_sig = np.abs(log_fc) >= params.lfc_cutoff

    if insufficient_replicates:
        # Statistical testing requires at least 3 replicates per condition.
        # Classify by logFC threshold only.
        p_col = np.full(len(count_matrix), np.nan)
        p_class = np.where(lfc_sig, "Sig.", "Not Sig.")
    else:
        # Use CPM for testing so that replicate library-size differences don't
        # inflate variance and wash out real signal.
        cond1_cpm = cpm_matrix[:, :n_cond1]
        cond2_cpm = cpm_matrix[:, n_cond1:]

        p_values = []
        for i in range(len(count_matrix)):
            _, p = stats.ttest_ind(cond1_cpm[i], cond2_cpm[i], equal_var=False)
            if np.isnan(p):
                p = 0.0 if cond1_cpm[i].mean() != cond2_cpm[i].mean() else 1.0
            p_values.append(p)

        p_values = np.array(p_values)
        _, p_adj, _, _ = multipletests(p_values, alpha=0.05, method="fdr_bh")
        p_col = np.round(p_adj, 4)
        # Classify by logFC threshold only; p-value is shown as reference
        p_class = np.where(lfc_sig, "Sig.", "Not Sig.")

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
