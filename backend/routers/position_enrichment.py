from typing import Dict, List, Any
from services.msa_service import perform_msa
from services.file_service import read_file, save_sequences
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from pathlib import Path
import os
import pandas as pd
import numpy as np
import Levenshtein
from services.constants import ColumnName
from numba import jit, prange
from typing import Any, Literal, Optional, List, Tuple
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import multiprocessing as mp


router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))

class PosEnrichInput(BaseModel):
    fadf_recluster_path: str = ""
    cluster_selection: int = 1
    seq_type: Literal["dna", "protein"] = "dna"
    output_format: str = "csv"
    max_sequences: int = 500  # Add this to control sampling

class PosEnrichResponse(BaseModel):
    status: str
    result: str  # filename
    num_sequences_analyzed: int
    sequence_length: int
    num_positions: int
    enrichment_range: Dict[str, float]
    # Matrix data for heatmap
    matrix: List[List[Optional[float]]]  # Change to Optional[float] to allow None
    residues: List[str]         # row labels
    positions: List[int]        # column labels

@router.post("/position-enrichment", response_model=PosEnrichResponse)
async def fa_pos_enrich_endpoint(params: PosEnrichInput):
    """
    Calculate positional enrichment analysis for a specific cluster.
    Returns:
    - CSV/TSV file with results
    - Matrix format for heatmap visualization
    - Summary statistics
    """
    print("ENDPOINT HIT: /position-enrichment")
     
    if not params.fadf_recluster_path:
        raise HTTPException(
            status_code=400,
            detail="fadf_recluster_path is required"
        )
    
    print("Starting positional enrichment analysis...")
    filepath = UPLOAD_DIR / params.fadf_recluster_path
    output_path = UPLOAD_DIR / f"pos_enrich_cluster_{params.cluster_selection}.{params.output_format}"
    
    try:
        print(f"[STEP 1] Reading file: {filepath}")
        # Read reclustered data
        fadf_recluster = read_file(filepath)
        print(f"[STEP 1] ✓ Loaded {len(fadf_recluster)} sequences")
        
        # Clean enrichment column if present
        if "Enrichment" in fadf_recluster.columns:
            fadf_recluster["Enrichment"] = pd.to_numeric(
                fadf_recluster["Enrichment"], 
                errors='coerce'
            )
        
        print(f"[STEP 2] Starting positional enrichment analysis...")
        # Perform positional enrichment analysis with sampling
        result_df = fa_pos_enrich(
            fadf_recluster=fadf_recluster,
            cluster_selection=params.cluster_selection,
            seq_type=params.seq_type,
            max_sequences=params.max_sequences
        )
        print(f"[STEP 2] ✓ Analysis complete. Result shape: {result_df.shape}")
        
        # Ensure AvEnrich is numeric
        result_df['AvEnrich'] = pd.to_numeric(result_df['AvEnrich'], errors='coerce')
        
        print(f"[STEP 3] Saving output to: {output_path}")
        # Save output file
        save_sequences(result_df, str(output_path), params.output_format)
        print(f"[STEP 3] ✓ File saved")
        
        print(f"[STEP 4] Creating pivot table for heatmap...")
        # Convert to matrix format for heatmap
        matrix_df = result_df.pivot(
            index='Residue', 
            columns='Position', 
            values='AvEnrich'
        )
        print(f"[STEP 4] ✓ Pivot table created. Shape: {matrix_df.shape}")
        
        # Extract matrix data with proper type handling
        residues = matrix_df.index.tolist()
        positions = [int(p) for p in matrix_df.columns.tolist()]
        
        # CRITICAL FIX: Replace NaN/Inf with None for JSON compliance
        print(f"[STEP 5] Cleaning matrix for JSON serialization...")
        matrix_clean = matrix_df.replace([np.nan, np.inf, -np.inf], None).values.tolist()
        
        # Alternative if you prefer 0 instead of null in frontend:
        # matrix_clean = matrix_df.fillna(0).replace([np.inf, -np.inf], 0).values.tolist()
        
        print(f"[STEP 5] ✓ Matrix cleaned: {len(residues)} residues x {len(positions)} positions")
        
        print(f"[STEP 6] Calculating summary statistics...")
        # Calculate summary statistics on valid values only
        enrichment_values = result_df["AvEnrich"].replace([np.inf, -np.inf], np.nan).dropna()
        
        if len(enrichment_values) > 0:
            enrichment_range = {
                "min": float(enrichment_values.min()),
                "max": float(enrichment_values.max()),
                "mean": float(enrichment_values.mean()),
                "median": float(enrichment_values.median())
            }
        else:
            enrichment_range = {
                "min": 0.0,
                "max": 0.0,
                "mean": 0.0,
                "median": 0.0
            }
        print(f"[STEP 6] ✓ Statistics: {enrichment_range}")
        
        print("=" * 60)
        print("REQUEST COMPLETED SUCCESSFULLY")
        print("=" * 60)
        
        return PosEnrichResponse(
            status="ok",
            result=output_path.name,
            num_sequences_analyzed=len(result_df["Position"].unique()) if len(result_df) > 0 else 0,
            sequence_length=int(result_df["Position"].max()) if len(result_df) > 0 else 0,
            num_positions=len(result_df["Position"].unique()) if len(result_df) > 0 else 0,
            enrichment_range=enrichment_range,
            matrix=matrix_clean,  # Use cleaned matrix
            residues=residues,
            positions=positions
        )
    
    except Exception as e:
        print("=" * 60)
        print(f"ERROR: {type(e).__name__}")
        print(f"Message: {str(e)}")
        import traceback
        traceback.print_exc()
        print("=" * 60)
        
        raise HTTPException(
            status_code=500,
            detail=f"Positional enrichment analysis failed: {str(e)}"
        )


def fa_pos_enrich(
    fadf_recluster: pd.DataFrame,
    cluster_selection: int = 1,
    seq_type: Literal["dna", "protein"] = "dna",
    max_sequences: int = 2000
) -> pd.DataFrame:
    """
    Calculate positional enrichment for sequences in a specific cluster.
    
    This function:
    1. Filters sequences from a specific cluster
    2. Removes sequences with ambiguous characters
    3. Performs Multiple Sequence Alignment (MSA) with sampling if needed
    4. Calculates average enrichment for each residue type at each position
    
    Args:
        fadf_recluster: DataFrame with reclustered sequences and enrichment scores
                       Must contain columns: Sequences, Cluster, Enrichment
        cluster_selection: Cluster ID to analyze
        seq_type: Type of sequences ("dna" or "protein")
        max_sequences: Maximum sequences for MSA (prevents timeout)
    
    Returns:
        DataFrame with columns:
        - Residue: The nucleotide/amino acid character
        - Position: Position in aligned sequence (1-indexed)
        - AvEnrich: Average enrichment for that residue at that position
    """
    print(f"\n{'='*60}")
    print(f"INSIDE fa_pos_enrich()")
    print(f"Input shape: {fadf_recluster.shape}")
    print(f"Cluster selection: {cluster_selection}")
    print(f"Sequence type: {seq_type}")
    print(f"Max sequences: {max_sequences}")
    print(f"{'='*60}\n")
    
    # Define acceptable alphabets for each sequence type
    seq_alphabets = {
        "dna": ["A", "C", "G", "T", "-"],
        "protein": ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
                   "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-"]
    }
    seq_alphabet = seq_alphabets[seq_type]
    
    # Create regex pattern for filtering out ambiguous characters
    alphabet_pattern = f"[^{''.join(seq_alphabet)}]"
    
    # If DNA sequences, convert U (uracil) to T (thymine)
    if seq_type == "dna":
        print("[SUBSTEP 1] Converting U to T for DNA sequences...")
        fadf_recluster = fadf_recluster.copy()
        fadf_recluster[ColumnName.SEQUENCES] = fadf_recluster[ColumnName.SEQUENCES].str.replace(
            "U", "T", regex=False
        )
        print("[SUBSTEP 1] ✓ Conversion complete")
    
    # Filter sequences based on three criteria:
    print("[SUBSTEP 2] Filtering sequences...")
    filtered_df = fadf_recluster[
        (fadf_recluster[ColumnName.CLUSTER] == cluster_selection) &
        (~fadf_recluster[ColumnName.SEQUENCES].str.contains(alphabet_pattern, regex=True, na=False)) &
        (fadf_recluster["Enrichment"].notna())
    ].copy()
    
    num_sequences = len(filtered_df)
    print(f"[SUBSTEP 2] Found {num_sequences} valid sequences for cluster {cluster_selection}")
    
    # Check if we have any sequences to analyze
    if num_sequences == 0:
        raise ValueError(
            f"No valid sequences found for cluster {cluster_selection}. "
            f"Check that: (1) cluster exists, (2) sequences have no ambiguous characters, "
            f"(3) sequences have enrichment scores."
        )
    
    # Sample if too many sequences
    if num_sequences > max_sequences:
        print(f"[WARNING] {num_sequences} sequences exceed limit of {max_sequences}")
        print(f"[INFO] Sampling top {max_sequences} sequences by enrichment...")
        filtered_df = filtered_df.nlargest(max_sequences, "Enrichment", keep='first')
        print(f"[INFO] ✓ Sampled {len(filtered_df)} sequences")
    
    # Ensure ID column exists (required by perform_msa function)
    print("[SUBSTEP 3] Checking ID column...")
    if ColumnName.ID not in filtered_df.columns:
        filtered_df[ColumnName.ID] = [f"seq_{i}" for i in range(len(filtered_df))]
        print(f"[SUBSTEP 3] Created IDs for {len(filtered_df)} sequences")
    print("[SUBSTEP 3] ✓ ID column ready")
    
    # Check if MSA is needed
    seq_lengths = filtered_df[ColumnName.SEQUENCES].str.len()
    all_same_length = seq_lengths.nunique() == 1
    
    if all_same_length:
        print(f"[INFO] All sequences have uniform length ({seq_lengths.iloc[0]}). Skipping MSA.")
        aligned_df = filtered_df.copy()
    else:
        # Perform Multiple Sequence Alignment
        print(f"[SUBSTEP 4] Starting Multiple Sequence Alignment...")
        print(f"  - Aligning {len(filtered_df)} sequences")
        print(f"  - This may take a while...")
        print(f"  - Start time: {pd.Timestamp.now()}")
        
        aligned_df = perform_msa(filtered_df)
        
        print(f"[SUBSTEP 4] ✓ MSA complete! End time: {pd.Timestamp.now()}")
        print(f"  - Aligned sequences shape: {aligned_df.shape}")
    
    # Ensure Enrichment column is present after MSA
    print("[SUBSTEP 5] Ensuring Enrichment column...")
    if "Enrichment" not in aligned_df.columns:
        aligned_df = aligned_df.merge(
            filtered_df[[ColumnName.ID, "Enrichment"]], 
            on=ColumnName.ID, 
            how='left'
        )
    print("[SUBSTEP 5] ✓ Enrichment column ready")
    
    # Extract only the columns we need for analysis
    print("[SUBSTEP 6] Processing aligned sequences...")
    posenrich_df = aligned_df[[ColumnName.SEQUENCES, "Enrichment"]].copy()
    posenrich_df["Enrichment"] = pd.to_numeric(posenrich_df["Enrichment"], errors='coerce')
    
    # Convert each sequence string into a list of characters
    seq_matrix = posenrich_df[ColumnName.SEQUENCES].apply(lambda x: list(x))
    
    # Get the length of aligned sequences
    seq_length = len(seq_matrix.iloc[0])
    print(f"[SUBSTEP 6] Sequence length: {seq_length}")
    
    # Initialize list to store results
    position_data = []
    
    print("[SUBSTEP 7] Calculating positional enrichment...")
    # Iterate through each position in the aligned sequences
    for pos_idx in range(seq_length):
        # Extract the character at this position for all sequences
        chars_at_pos = seq_matrix.apply(lambda x: x[pos_idx])
        
        # Create temporary DataFrame pairing each character with its enrichment
        temp_df = pd.DataFrame({
            'char': chars_at_pos.values,
            'enrichment': posenrich_df["Enrichment"].values
        })
        
        # Group by character and calculate mean enrichment
        avg_enrich = temp_df.groupby('char')['enrichment'].mean()
        
        # For each possible residue in the alphabet, store the result
        for residue in seq_alphabet:
            enrich_value = avg_enrich.get(residue, np.nan)
            position_data.append({
                'Residue': residue,
                'Position': pos_idx + 1,  # Convert to 1-indexed
                'AvEnrich': round(float(enrich_value), 3) if pd.notna(enrich_value) else np.nan
            })
    
    print(f"[SUBSTEP 7] ✓ Calculated enrichment for {seq_length} positions")
    
    # Create final results DataFrame
    pe_results = pd.DataFrame(position_data)
    
    # Ensure AvEnrich is float type
    pe_results['AvEnrich'] = pe_results['AvEnrich'].astype(float)
    
    print(f"✓ fa_pos_enrich complete. Result shape: {pe_results.shape}\n")
    
    return pe_results