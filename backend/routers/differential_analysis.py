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
    p_cutoff: float = 0.1
    output_format: str = "csv"

class DifferentialExpressionResponse(BaseModel):
    status: str
    result: str
    total_sequences: int
    significant_sequences: int
    upregulated: int
    downregulated: int
    p_cutoff: float

def normalize_counts(count_matrix, lib_sizes):
    """
    Normalize counts to CPM (Counts Per Million).
    
    Args:
        count_matrix: numpy array of shape (n_sequences, n_samples)
        lib_sizes: list of library sizes for each sample
        
    Returns:
        Normalized CPM matrix
    """
    lib_sizes = np.array(lib_sizes)
    cpm = (count_matrix / lib_sizes) * 1e6
    return cpm

def calculate_log_cpm(cpm_matrix, prior_count=2):
    """
    Calculate log2 CPM with prior count (similar to edgeR's cpm function).
    
    Args:
        cpm_matrix: CPM normalized counts
        prior_count: Prior count to add before log transformation
        
    Returns:
        log2 CPM values
    """
    return np.log2(cpm_matrix + prior_count)

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
        raise HTTPException(
            status_code=400,
            detail="Both conditions must have at least one file"
        )
    
    if len(params.cond1_paths) > 26 or len(params.cond2_paths) > 26:
        raise HTTPException(
            status_code=400,
            detail="Maximum 26 replicates per condition (A-Z suffixes)"
        )
    
    if not (0 < params.p_cutoff < 1):
        raise HTTPException(
            status_code=400,
            detail="p_cutoff must be between 0 and 1"
        )
    
    ## STEP 1: LOAD AND FORMAT CONDITION 1
    fadf_cond1 = []
    original_lib_sizes_cond1 = []
    
    for i, file_path in enumerate(params.cond1_paths):
        input_file = UPLOAD_DIR / file_path
        if not input_file.exists():
            raise HTTPException(
                status_code=404,
                detail=f"File not found: {file_path}"
            )
        
        df = read_file(str(input_file))
        
        # Validate required columns
        if ColumnName.SEQUENCES not in df.columns or ColumnName.READS not in df.columns:
            raise HTTPException(
                status_code=400,
                detail=f"File {file_path} must contain 'Sequences' and 'Reads' columns"
            )
        
        # Store original library size
        original_lib_sizes_cond1.append(df[ColumnName.READS].sum())
        
        # Create suffix (A, B, C, ...)
        suffix = chr(ord('A') + i)
        
        # Select and rename columns
        df_formatted = df[[ColumnName.SEQUENCES, ColumnName.READS]].copy()
        df_formatted = df_formatted.rename(columns={
            ColumnName.READS: f"Reads.{suffix}"
        })
        
        fadf_cond1.append(df_formatted)
    
    ## STEP 2: LOAD AND FORMAT CONDITION 2
    fadf_cond2 = []
    original_lib_sizes_cond2 = []
    
    for i, file_path in enumerate(params.cond2_paths):
        input_file = UPLOAD_DIR / file_path
        if not input_file.exists():
            raise HTTPException(
                status_code=404,
                detail=f"File not found: {file_path}"
            )
        
        df = read_file(str(input_file))
        
        if ColumnName.SEQUENCES not in df.columns or ColumnName.READS not in df.columns:
            raise HTTPException(
                status_code=400,
                detail=f"File {file_path} must contain 'Sequences' and 'Reads' columns"
            )
        
        # Store original library size
        original_lib_sizes_cond2.append(df[ColumnName.READS].sum())
        
        # Create suffix (a, b, c, ...)
        suffix = chr(ord('a') + i)
        
        # Select and rename columns
        df_formatted = df[[ColumnName.SEQUENCES, ColumnName.READS]].copy()
        df_formatted = df_formatted.rename(columns={
            ColumnName.READS: f"Reads.{suffix}"
        })
        
        fadf_cond2.append(df_formatted)
    
    ## STEP 3: MERGE CONDITION 1 (INNER JOIN)
    cond1_merged = fadf_cond1[0]
    for i in range(1, len(fadf_cond1)):
        cond1_merged = pd.merge(
            cond1_merged,
            fadf_cond1[i],
            on=ColumnName.SEQUENCES,
            how='inner'
        )
    
    ## STEP 4: MERGE CONDITION 2 (INNER JOIN)
    cond2_merged = fadf_cond2[0]
    for i in range(1, len(fadf_cond2)):
        cond2_merged = pd.merge(
            cond2_merged,
            fadf_cond2[i],
            on=ColumnName.SEQUENCES,
            how='inner'
        )
    
    ## STEP 5: MERGE BOTH CONDITIONS (INNER JOIN)
    cond_merged = pd.merge(
        cond1_merged,
        cond2_merged,
        on=ColumnName.SEQUENCES,
        how='inner'
    )
    
    if len(cond_merged) < 10:
        raise HTTPException(
            status_code=400,
            detail=f"Too few shared sequences ({len(cond_merged)}). Need at least 10 for differential analysis."
        )
    
    ## STEP 6: CREATE COUNT MATRIX
    # Get read columns
    cond1_cols = [f"Reads.{chr(ord('A') + i)}" for i in range(len(fadf_cond1))]
    cond2_cols = [f"Reads.{chr(ord('a') + i)}" for i in range(len(fadf_cond2))]
    all_read_cols = cond1_cols + cond2_cols
    
    # Extract count matrix (sequences x samples)
    count_matrix = cond_merged[all_read_cols].values
    sequences = cond_merged[ColumnName.SEQUENCES].values
    
    ## STEP 7: NORMALIZE COUNTS TO CPM
    all_lib_sizes = original_lib_sizes_cond1 + original_lib_sizes_cond2
    cpm_matrix = normalize_counts(count_matrix, all_lib_sizes)
    
    # Split into conditions
    n_cond1 = len(fadf_cond1)
    cond1_cpm = cpm_matrix[:, :n_cond1]
    cond2_cpm = cpm_matrix[:, n_cond1:]
    
    ## STEP 8: CALCULATE LOG FOLD CHANGE
    # Calculate mean CPM for each condition (with pseudocount)
    pseudocount = 1
    cond1_mean_cpm = cond1_cpm.mean(axis=1) + pseudocount
    cond2_mean_cpm = cond2_cpm.mean(axis=1) + pseudocount
    
    # Calculate log2 fold change (cond2 vs cond1)
    log_fc = np.log2(cond2_mean_cpm / cond1_mean_cpm)
    
    # Calculate average log CPM
    log_cpm = np.log2((cond1_mean_cpm + cond2_mean_cpm) / 2)
    
    ## STEP 9: STATISTICAL TESTING
    p_values = []
    
    # Choose test based on number of replicates
    if n_cond1 >= 3 and len(fadf_cond2) >= 3:
        # Use t-test for 3+ replicates
        for i in range(len(count_matrix)):
            _, p_val = stats.ttest_ind(cond1_cpm[i], cond2_cpm[i])
            p_values.append(p_val)
    else:
        # Use Mann-Whitney U test for fewer replicates (more robust)
        for i in range(len(count_matrix)):
            try:
                _, p_val = stats.mannwhitneyu(
                    cond1_cpm[i], 
                    cond2_cpm[i],
                    alternative='two-sided'
                )
                p_values.append(p_val)
            except ValueError:
                # Handle cases where test cannot be computed
                p_values.append(1.0)
    
    p_values = np.array(p_values)
    
    ## STEP 10: MULTIPLE TESTING CORRECTION (BENJAMINI-HOCHBERG)
    # Replace NaN p-values with 1.0
    p_values = np.nan_to_num(p_values, nan=1.0)
    
    # Apply Benjamini-Hochberg FDR correction
    reject, p_adjusted, _, _ = multipletests(
        p_values,
        alpha=params.p_cutoff,
        method='fdr_bh'
    )
    
    ## STEP 11: FORMAT RESULTS
    d_test = pd.DataFrame({
        'logFC': np.round(log_fc, 3),
        'logCPM': np.round(log_cpm, 3),
        'PValue': np.round(p_adjusted, 3),
        'PClass': ['Sig.' if p < params.p_cutoff else 'Not Sig.' for p in p_adjusted],
        'Sequence': sequences
    })
    
    # Sort by p-value
    d_test = d_test.sort_values('PValue').reset_index(drop=True)
    
    ## STEP 12: SAVE RESULTS
    base = Path(params.cond1_paths[0]).stem
    output = UPLOAD_DIR / f"{base}_diff_expr.{params.output_format}"
    
    save_sequences(d_test, output, params.output_format)
    
    
    return DifferentialExpressionResponse(
        status="ok",
        result=output.name,
        p_cutoff=params.p_cutoff
    )