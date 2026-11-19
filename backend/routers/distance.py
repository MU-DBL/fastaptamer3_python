from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from pathlib import Path
from typing import Optional
import numpy as np
import pandas as pd
from Levenshtein import distance as levenshtein_distance
from rapidfuzz import process
from rapidfuzz.distance import Levenshtein
import os
from services.file_service import read_file, save_sequences, format_rpu, scaling_factor_of
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
    
    fa_df = read_file(input_file)
    
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
    
    # Calculate edit distance for each sequence in one multithreaded call
    # (cdist splits work across threads by query row, so the sequences are the queries)
    fa_df[ColumnName.DISTANCE] = process.cdist(
        fa_df[ColumnName.SEQUENCES].str.upper().tolist(), [query_upper],
        scorer=Levenshtein.distance, workers=-1
    )[:, 0].astype(np.int64)

    # Sort by distance (ascending)
    result_df = fa_df.sort_values(by=ColumnName.DISTANCE).reset_index(drop=True)

    # U of the input's RPU, kept in the IDs (None for inputs that do not record it)
    u = scaling_factor_of(result_df[ColumnName.U]) if ColumnName.U in result_df.columns else None

    # Build ID column with all metadata including Distance (column-wise; a row-wise apply is slow)
    result_df[ColumnName.ID] = (
        f"{ColumnName.RANK}=" + result_df[ColumnName.RANK].astype(np.int64).astype(str)
        + f";{ColumnName.READS}=" + result_df[ColumnName.READS].astype(np.int64).astype(str)
        + f";{ColumnName.RPU}=" + format_rpu(result_df[ColumnName.RPU])
        + (f";{ColumnName.U}={u}" if u is not None else "")
        + f";{ColumnName.DISTANCE}=" + result_df[ColumnName.DISTANCE].astype(str)
    )
    
    # Select and reorder columns
    result_df = result_df[[
        ColumnName.ID,
        ColumnName.RANK,
        ColumnName.READS,
        ColumnName.RPU,
        ColumnName.SEQUENCES,
        ColumnName.DISTANCE
    ]].copy()
    
    # Save output
    base = Path(params.input_path).stem
    output = UPLOAD_DIR / f"{base}_distance.{params.output_format}"
    
    save_sequences(result_df, output, params.output_format)
    
    return {
        "status": "ok",
        "result": output.name
    }