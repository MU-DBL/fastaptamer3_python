from services.file_service import read_file
from fastapi import APIRouter, HTTPException
from pathlib import Path
import os
import pandas as pd
import numpy as np
from typing import Any, List, Optional
from pydantic import BaseModel
from services.constants import ColumnName
import json

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))

class PHMMInput(BaseModel):
    input_path: str = ""
    num_sequences: int = 10
    sequence_length: int = 70
    output_format_phmm: str = "txt"
    output_format_simulation: str = "fasta"
    pseudocount_method: str = "laplace" 

class PHMMResponse(BaseModel):
    status: str
    simulated_sequences: List[str]
    num_unique_sequences: int
    output_file_phmm: str    
    output_file_simulation: str
    phmm_stats: dict

@router.post("/cluster-phmm-simulate", response_model=PHMMResponse)
async def phmm_create_and_simulate(params: PHMMInput):
    if not params.input_path:
        raise HTTPException(status_code=400, detail="input_path is required")
    
    filepath = UPLOAD_DIR / params.input_path
    
    if not filepath.exists():
        raise HTTPException(status_code=404, detail=f"File not found: {params.input_path}")
    
    try:
        # Step 1: Read MSA file
        df_msa = read_file(filepath)
        
        # Step 2: Create PHMM from MSA
        phmm, phmm_stats = create_phmm_from_msa(
            df_msa, 
            pseudocount_method=params.pseudocount_method
        )
        
        # Step 3: Simulate sequences using PHMM
        simulated_sequences = simulate_sequences_from_phmm(
            phmm=phmm,
            num_sequences=params.num_sequences,
            sequence_length=params.sequence_length
        )
        
        # Initialize variables to ensure they exist even if if-statements fail
        output_file_sim_path = None
        output_file_phmm_path = None
        base_name = filepath.stem

        # Logic for Saving Simulation
        if params.output_format_simulation:
            sim_filename = f"{base_name}_phmm_simulated.{params.output_format_simulation}"
            output_path_sim = UPLOAD_DIR / sim_filename
            
            saved_file_sim = save_simulated_sequences(
                sequences=simulated_sequences,
                output_path=output_path_sim,
                output_format=params.output_format_simulation
            )
            # Store the name/path string, not the file object
            output_file_sim_path = saved_file_sim.name if hasattr(saved_file_sim, 'name') else str(output_path_sim)

        # Logic for Saving PHMM (Optional addition based on your Input Model)
        if params.output_format_phmm:
            phmm_filename = f"{base_name}_phmm_model.{params.output_format_phmm}"
            output_path_phmm = UPLOAD_DIR / phmm_filename
            
            # Assuming a save_phmm function exists
            save_phmm(phmm, output_path_phmm) 
            output_file_phmm_path = str(phmm_filename)

        return {
            "status": "ok",
            "simulated_sequences": simulated_sequences,
            "num_unique_sequences": len(set(simulated_sequences)), # Ensure this is calculated on the unique set
            "output_file_phmm": output_file_phmm_path,
            "output_file_simulation": output_file_sim_path,
            "phmm_stats": phmm_stats
        }
    
    except ValueError as ve:
        raise HTTPException(status_code=400, detail=str(ve))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"PHMM simulation failed: {str(e)}")

def create_phmm_from_msa(df_msa: pd.DataFrame, pseudocount_method: str = "laplace") -> tuple:

    sequences = df_msa[ColumnName.SEQUENCES].tolist()
    
    # Convert sequences to matrix
    msa_matrix = np.array([list(seq.lower()) for seq in sequences])
    n_seqs, seq_len = msa_matrix.shape
    
    # Identify match columns (columns with < 50% gaps)
    gap_threshold = 0.5
    match_columns = []
    
    for col_idx in range(seq_len):
        column = msa_matrix[:, col_idx]
        gap_fraction = np.sum(column == '-') / len(column)
        if gap_fraction < gap_threshold:
            match_columns.append(col_idx)
    
    n_match_states = len(match_columns)
    
    # Build emission probabilities for match states
    match_emissions = []
    alphabet = ['a', 'c', 'g', 't', 'n']  # DNA alphabet
    
    for match_col in match_columns:
        column = msa_matrix[:, match_col]
        # Remove gaps for emission calculation
        column_no_gaps = column[column != '-']
        
        # Count character frequencies
        emissions = {}
        for char in alphabet:
            count = np.sum(column_no_gaps == char)
            
            # Add pseudocounts based on method
            if pseudocount_method == "laplace":
                count += 1
                total = len(column_no_gaps) + len(alphabet)
            elif pseudocount_method == "background":
                count += 0.25  # Uniform background
                total = len(column_no_gaps) + 1
            else:
                total = len(column_no_gaps)
            
            emissions[char] = count / total if total > 0 else 1.0 / len(alphabet)
        
        match_emissions.append(emissions)
    
    # Build transition probabilities (simplified model)
    # Transitions: Match->Match, Match->Insert, Match->Delete
    transitions = calculate_transition_probabilities(msa_matrix, match_columns)
    
    phmm = {
        "match_states": n_match_states,
        "match_emissions": match_emissions,
        "transitions": transitions,
        "alphabet": alphabet
    }
    
    stats = {
        "num_sequences": n_seqs,
        "alignment_length": seq_len,
        "num_match_states": n_match_states,
        "pseudocount_method": pseudocount_method
    }
    
    return phmm, stats


def calculate_transition_probabilities(msa_matrix: np.ndarray, match_columns: List[int]) -> dict:
    n_seqs = msa_matrix.shape[0]
    
    # Count transitions
    mm_count = 0  # Match to Match
    md_count = 0  # Match to Delete
    mi_count = 0  # Match to Insert
    
    for seq_idx in range(n_seqs):
        for i in range(len(match_columns) - 1):
            curr_col = match_columns[i]
            next_col = match_columns[i + 1]
            
            curr_char = msa_matrix[seq_idx, curr_col]
            next_char = msa_matrix[seq_idx, next_col]
            
            if curr_char != '-':
                if next_char != '-':
                    mm_count += 1
                else:
                    md_count += 1
                
                # Check for insert states between match columns
                if next_col - curr_col > 1:
                    has_insert = any(msa_matrix[seq_idx, curr_col+1:next_col] != '-')
                    if has_insert:
                        mi_count += 1
    
    total = mm_count + md_count + mi_count
    
    if total == 0:
        return {"M->M": 0.9, "M->D": 0.05, "M->I": 0.05}
    
    return {
        "M->M": mm_count / total,
        "M->D": md_count / total,
        "M->I": mi_count / total
    }


def simulate_sequences_from_phmm(phmm: dict, num_sequences: int, sequence_length: int) -> List[str]:

    simulated = set()
    attempts = 0
    max_attempts = num_sequences * 3  # Try up to 3x to get unique sequences
    
    while len(simulated) < num_sequences and attempts < max_attempts:
        sequence = generate_single_sequence(phmm, sequence_length)
        simulated.add(sequence.upper())
        attempts += 1
    
    return sorted(list(simulated))


def generate_single_sequence(phmm: dict, target_length: int) -> str:
    sequence = []
    match_emissions = phmm["match_emissions"]
    transitions = phmm["transitions"]
    
    # Generate sequence by sampling from match states
    state_idx = 0
    while len(sequence) < target_length and state_idx < len(match_emissions):
        emissions = match_emissions[state_idx]
        
        # Sample character from emission probabilities
        chars = list(emissions.keys())
        probs = list(emissions.values())
        
        # Normalize probabilities
        total_prob = sum(probs)
        probs = [p / total_prob for p in probs]
        
        char = np.random.choice(chars, p=probs)
        sequence.append(char)
        
        # Decide next state based on transition probabilities
        rand = np.random.random()
        if rand < transitions["M->I"]:
            # Insert state - add extra character
            if len(sequence) < target_length:
                sequence.append(np.random.choice(chars, p=probs))
        elif rand < transitions["M->I"] + transitions["M->D"]:
            # Delete state - skip next match state
            state_idx += 1
        
        state_idx += 1
    
    # Pad or trim to target length
    while len(sequence) < target_length:
        # Use last match state emissions
        emissions = match_emissions[-1]
        chars = list(emissions.keys())
        probs = [p / sum(emissions.values()) for p in emissions.values()]
        sequence.append(np.random.choice(chars, p=probs))
    
    return ''.join(sequence[:target_length])


def save_simulated_sequences(sequences: List[str], output_path: Path, output_format: str) -> Path:
    if output_format == "fasta":
        with open(output_path, 'w') as f:
            for idx, seq in enumerate(sequences, 1):
                f.write(f">simulated_seq_{idx}\n{seq}\n")
    
    else:
        raise ValueError(f"Unsupported output format: {output_format}")
    
    return output_path

def save_phmm(phmm: Any, output_path: Path) -> Path:
    with open(output_path, "w") as f:
        f.write(str(phmm))

    return output_path