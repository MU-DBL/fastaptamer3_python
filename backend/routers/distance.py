from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from pathlib import Path
from typing import Optional
import pandas as pd
from Levenshtein import distance as levenshtein_distance
import os
from services.file_service import read_file, save_sequences
from services.constants import ColumnName

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))

class DistanceInput(BaseModel):
    input_path: str  # FASTA file path
    query_sequence: str
    output_format: str = "csv"

@router.post("/sequence-distance")
async def sequence_distance(params: DistanceInput):
    # Validate inputs
    if not params.query_sequence:
        raise HTTPException(status_code=400, detail="Query sequence is required")
    
    # Load FASTA file
    input_file = UPLOAD_DIR / params.input_path
    if not input_file.exists():
        raise HTTPException(status_code=404, detail=f"File not found: {params.input_path}")
    
    fa_df = read_file(str(input_file))
    
    # Validate required columns
    required_cols = [ColumnName.ID, ColumnName.RANK, ColumnName.READS, 
                     ColumnName.RPU, ColumnName.SEQUENCES]
    missing_cols = [col for col in required_cols if col not in fa_df.columns]
    if missing_cols:
        raise HTTPException(
            status_code=400, 
            detail=f"Missing required columns: {', '.join(missing_cols)}"
        )
    
    # Normalize sequences to uppercase
    query_upper = params.query_sequence.upper()
    
    # Calculate edit distance for each sequence
    fa_df[ColumnName.DISTANCE] = fa_df[ColumnName.SEQUENCES].apply(
        lambda seq: levenshtein_distance(seq.upper(), query_upper)
    )
    
    # Select and reorder columns
    result_df = fa_df[[
        ColumnName.ID,
        ColumnName.RANK,
        ColumnName.READS,
        ColumnName.RPU,
        ColumnName.SEQUENCES,
        ColumnName.DISTANCE
    ]].copy()
    
    # Sort by distance (ascending)
    result_df = result_df.sort_values(by=ColumnName.DISTANCE).reset_index(drop=True)
    
    # Save output
    base = Path(params.input_path).stem
    output = UPLOAD_DIR / f"{base}_distance.{params.output_format}"
    
    save_sequences(result_df, output, params.output_format)
    
    return {
        "status": "ok",
        "result": output.name
    }