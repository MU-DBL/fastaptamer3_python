import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import tempfile
import os
from pathlib import Path
from routers.position_enrichment import perform_msa
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from services.constants import ColumnName
from services.file_service import save_sequences,read_file
from fastapi import APIRouter, HTTPException
from pathlib import Path
from typing import List, Dict, Any
import numpy as np
from sklearn.metrics import mutual_info_score
import subprocess

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))

class ClusterMSAInput(BaseModel):
    input_path: str = ""
    output_format: str = "fasta"
    seq_type: str = "dna"
    cluster_selected: int = 1

@router.post("/clustermsa")
async def cluster_msa_endpoint(params: ClusterMSAInput):
    """
    Perform MSA on a selected cluster from clustered sequences.
    """
    if not params.input_path:
        raise HTTPException(status_code=400, detail="input_path is required")
    
    filepath = UPLOAD_DIR / params.input_path
    base_name = filepath.stem
    output_path = UPLOAD_DIR / f"{base_name}_cluster{params.cluster_selected}_msa.{params.output_format}"
    
    try:
        # Read clustered data
        df_cluster = read_file(filepath) 
        
        # Filter for selected cluster
        cluster_df = df_cluster[df_cluster[ColumnName.CLUSTER] == params.cluster_selected].copy()
        
        if len(cluster_df) == 0:
            raise ValueError(f"No sequences found for cluster {params.cluster_selected}")
        
        # Single sequence - no alignment needed
        if len(cluster_df) == 1:
            result_df = cluster_df
        else:
            # Perform MSA
            result_df = perform_msa(cluster_df)
            
        # Relocate Sequences column to end and sort
        cols = [col for col in result_df.columns if col != ColumnName.SEQUENCES]
        cols.append(ColumnName.SEQUENCES)
        result_df = result_df[cols]
        
        if ColumnName.RANK in result_df.columns:
            result_df = result_df.sort_values(ColumnName.RANK).reset_index(drop=True)
        
        # Save output
        save_sequences(result_df, str(output_path), params.output_format)
        
        return {
            "status": "ok",
            "result": output_path.name,
            "num_sequences": len(result_df)
        }
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Cluster MSA failed: {str(e)}")



class MSAAnalysisInput(BaseModel):
    input_path: str = ""

class EntropyResponse(BaseModel):
    status: str
    positions: List[int]
    entropy_values: List[float]
    total_positions: int

class MutualInfoResponse(BaseModel):
    status: str
    mutual_info_matrix: List[List[float]]
    positions: List[int]
    total_positions: int

@router.post("/cluster-msa-entropy", response_model=EntropyResponse)
async def calculate_msa_entropy(params: MSAAnalysisInput):

    if not params.input_path:
        raise HTTPException(status_code=400, detail="input_path is required")
    
    filepath = UPLOAD_DIR / params.input_path
    
    if not filepath.exists():
        raise HTTPException(status_code=404, detail=f"File not found: {params.input_path}")
    
    try:
        df = read_file(filepath)
        
        # Convert sequences to matrix
        msa_matrix = sequences_to_matrix(df)
        
        # Calculate entropy for each position (column)
        entropy_values = []
        for col_idx in range(msa_matrix.shape[1]):
            column = msa_matrix[:, col_idx]
            # Calculate Shannon entropy in nats (natural log)
            col_entropy = calculate_shannon_entropy(column)
            entropy_values.append(round(col_entropy, 6))
        
        positions = list(range(1, len(entropy_values) + 1))
        
        return {
            "status": "ok",
            "positions": positions,
            "entropy_values": entropy_values,
            "total_positions": len(positions)
        }
    
    except ValueError as ve:
        raise HTTPException(status_code=400, detail=str(ve))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Entropy calculation failed: {str(e)}")


# ============= Mutual Information API =============
@router.post("/cluster-msa-mutinfo", response_model=MutualInfoResponse)
async def calculate_msa_mutual_info(params: MSAAnalysisInput):

    if not params.input_path:
        raise HTTPException(status_code=400, detail="input_path is required")
    
    filepath = UPLOAD_DIR / params.input_path
    
    if not filepath.exists():
        raise HTTPException(status_code=404, detail=f"File not found: {params.input_path}")
    
    try:
        # Read MSA file
        df = read_file(filepath)
        
        # Convert sequences to matrix
        msa_matrix = sequences_to_matrix(df)
        
        n_positions = msa_matrix.shape[1]
        
        # Calculate pairwise mutual information
        mi_matrix = np.zeros((n_positions, n_positions))
        
        for i in range(n_positions):
            for j in range(n_positions):
                if i == j:
                    # Self-information equals entropy
                    mi_matrix[i, j] = calculate_shannon_entropy(msa_matrix[:, i])
                else:
                    # Mutual information between positions
                    mi_matrix[i, j] = calculate_mutual_information(
                        msa_matrix[:, i], 
                        msa_matrix[:, j]
                    )
        
        # Round values for readability
        mi_matrix = np.round(mi_matrix, 6)
        
        positions = list(range(1, n_positions + 1))
        
        return {
            "status": "ok",
            "mutual_info_matrix": mi_matrix.tolist(),
            "positions": positions,
            "total_positions": n_positions
        }
    
    except ValueError as ve:
        raise HTTPException(status_code=400, detail=str(ve))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Mutual information calculation failed: {str(e)}")


def sequences_to_matrix(df: pd.DataFrame):
    sequences = df[ColumnName.SEQUENCES].tolist()
    
    # Check all sequences have same length (requirement for MSA)
    seq_lengths = [len(seq) for seq in sequences]
    if len(set(seq_lengths)) > 1:
        raise ValueError("Sequences have different lengths - not a valid MSA")
    
    # Convert to matrix
    msa_matrix = np.array([list(seq) for seq in sequences])
    
    return msa_matrix


def calculate_shannon_entropy(column: np.ndarray) -> float:
    # Get character frequencies
    unique, counts = np.unique(column, return_counts=True)
    probabilities = counts / len(column)
    
    # Calculate entropy (using natural log for nats)
    # Filter out zero probabilities to avoid log(0)
    probabilities = probabilities[probabilities > 0]
    ent = -np.sum(probabilities * np.log(probabilities))
    
    return float(ent)


def calculate_mutual_information(col1: np.ndarray, col2: np.ndarray) -> float:
    # Create contingency table
    unique_col1 = np.unique(col1)
    unique_col2 = np.unique(col2)
    
    # Map characters to indices
    map1 = {char: idx for idx, char in enumerate(unique_col1)}
    map2 = {char: idx for idx, char in enumerate(unique_col2)}
    
    col1_encoded = np.array([map1[char] for char in col1])
    col2_encoded = np.array([map2[char] for char in col2])
    
    # Calculate MI using sklearn
    mi = mutual_info_score(col1_encoded, col2_encoded)
    
    return float(mi)