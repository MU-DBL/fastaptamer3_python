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
    cluster_selection: Optional[int] = None
    cluster_column: Optional[str] = None   # e.g. "Cluster" or "Cluster.a"; None = auto-detect
    seq_type: Literal["dna", "protein"] = "dna"
    output_format: str = "csv"
    max_sequences: int = 2000

class PosEnrichResponse(BaseModel):
    status: str
    result: str  # filename
    enrichment_matrix: List[List[Optional[float]]]  # Change to Optional[float] to allow None
    residues: List[str]         # row labels
    positions: List[int]        # column labels
    avg_enrichment: List[Optional[float]]  # average enrichment per position
            
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
    cluster_label = str(params.cluster_selection) if params.cluster_selection is not None else "all"
    output_path = UPLOAD_DIR / f"pos_enrich_cluster_{cluster_label}.{params.output_format}"
    
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
            cluster_column=params.cluster_column,
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
        
        # print(f"[STEP 4] Creating pivot table for heatmap...")
        # Convert to matrix format for heatmap
        matrix_df = result_df.pivot(
            index='Residue', 
            columns='Position', 
            values='AvEnrich'
        )
        # print(f"[STEP 4] ✓ Pivot table created. Shape: {matrix_df.shape}")
        
        # Extract matrix data with proper type handling
        residues = matrix_df.index.tolist()
        positions = [int(p) for p in matrix_df.columns.tolist()]       
        
        clean_df = (
            matrix_df
            .replace([np.inf, -np.inf], np.nan)  
            .fillna(0)                            
        )

        matrix_clean = clean_df.values.tolist()
        avg_enrichment = clean_df.mean(axis=0).tolist()

        return PosEnrichResponse(
            status="ok",
            result=output_path.name,
            residues=residues,
            positions=positions,
            avg_enrichment = avg_enrichment,
            enrichment_matrix=matrix_clean
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
            detail=f"Positional enrichment analysis failed ({type(e).__name__}): {str(e)}"
        )


def fa_pos_enrich(
    fadf_recluster: pd.DataFrame,
    cluster_selection: Optional[int] = None,
    cluster_column: Optional[str] = None,
    seq_type: Literal["dna", "protein"] = "dna",
    max_sequences: int = 2000
) -> pd.DataFrame:
    """
    Calculate positional enrichment for sequences.

    Supports recluster CSV (Cluster column), sequence enrichment CSV (Cluster.a),
    or any file with an Enrichment column (no cluster filter when cluster_selection is None).
    """
    
    # Define acceptable alphabets for each sequence type, including non-canonical/ambiguous
    # residues (IUPAC nucleotide ambiguity codes; extended amino acid codes) so the heatmap's
    # y-axis can show them instead of silently dropping sequences that contain them.
    seq_alphabets = {
        "dna": ["A", "C", "G", "T",
                "R", "Y", "W", "S", "M", "K", "B", "D", "H", "V", "N", "-"],
        "protein": ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                   "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
                   "B", "Z", "X", "J", "U", "O", "*", "-"]
    }
    seq_alphabet = seq_alphabets[seq_type]
    
    # Create regex pattern for filtering out ambiguous characters
    alphabet_pattern = f"[^{''.join(seq_alphabet)}]"
    
    # If DNA sequences, convert U (uracil) to T (thymine)
    if seq_type == "dna":
        fadf_recluster = fadf_recluster.copy()
        fadf_recluster[ColumnName.SEQUENCES] = fadf_recluster[ColumnName.SEQUENCES].str.replace(
            "U", "T", regex=False
        )
    
    # Auto-detect cluster column if not specified
    if cluster_column is None:
        for candidate in [ColumnName.CLUSTER, "Cluster.a"]:
            if candidate in fadf_recluster.columns:
                cluster_column = candidate
                break

    # Build filter mask
    valid_seqs = ~fadf_recluster[ColumnName.SEQUENCES].str.contains(alphabet_pattern, regex=True, na=False)
    has_enrichment = fadf_recluster["Enrichment"].notna()

    if cluster_column is not None and cluster_selection is not None:
        in_cluster = fadf_recluster[cluster_column] == cluster_selection
        filtered_df = fadf_recluster[in_cluster & valid_seqs & has_enrichment].copy()
    else:
        filtered_df = fadf_recluster[valid_seqs & has_enrichment].copy()

    num_sequences = len(filtered_df)

    # Check if we have any sequences to analyze
    if num_sequences == 0:
        cluster_info = f"cluster {cluster_selection}" if cluster_selection is not None else "the dataset"
        raise ValueError(
            f"No valid sequences found for {cluster_info}. "
            f"Check that sequences have no ambiguous characters and enrichment scores are present."
        )
    
    # Sample if too many sequences
    if num_sequences > max_sequences:
        filtered_df = filtered_df.nlargest(max_sequences, "Enrichment", keep='first')
    
    # Ensure ID column exists (required by perform_msa function)
    if ColumnName.ID not in filtered_df.columns:
        filtered_df[ColumnName.ID] = [f"seq_{i}" for i in range(len(filtered_df))]
    
    # Check if MSA is needed
    seq_lengths = filtered_df[ColumnName.SEQUENCES].str.len()
    all_same_length = seq_lengths.nunique() == 1
    
    if all_same_length:
        aligned_df = filtered_df.copy()
    else:
        # Perform Multiple Sequence Alignment
        aligned_df = perform_msa(filtered_df)

    
    # Ensure Enrichment column is present after MSA
    if "Enrichment" not in aligned_df.columns:
        aligned_df = aligned_df.merge(
            filtered_df[[ColumnName.ID, "Enrichment"]], 
            on=ColumnName.ID, 
            how='left'
        )
    
    # Extract only the columns we need for analysis
    posenrich_df = aligned_df[[ColumnName.SEQUENCES, "Enrichment"]].copy()
    posenrich_df["Enrichment"] = pd.to_numeric(posenrich_df["Enrichment"], errors='coerce')
    
    # Convert each sequence string into a list of characters
    seq_matrix = posenrich_df[ColumnName.SEQUENCES].apply(lambda x: list(x))
    
    # Get the length of aligned sequences
    seq_length = len(seq_matrix.iloc[0])
    
    # Initialize list to store results
    position_data = []
    
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
    
    
    # Create final results DataFrame
    pe_results = pd.DataFrame(position_data)
    
    # Ensure AvEnrich is float type
    pe_results['AvEnrich'] = pe_results['AvEnrich'].astype(float)    
    return pe_results
