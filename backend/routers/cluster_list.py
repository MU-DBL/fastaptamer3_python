from fastapi import APIRouter, HTTPException, UploadFile, File
from pathlib import Path
import os
import pandas as pd
from typing import List
from pydantic import BaseModel
from services.constants import ColumnName
from services.file_service import read_file

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))

class ClusterListInput(BaseModel):
    input_path: str = ""


@router.post("/cluster-list")
async def get_cluster_list(params: ClusterListInput):
    if not params.input_path:
        raise HTTPException(status_code=400, detail="input_path is required")
    
    filepath = UPLOAD_DIR / params.input_path
    
    if not filepath.exists():
        raise HTTPException(status_code=404, detail=f"File not found: {params.input_path}")
    
    try:
        df = read_file(filepath)
        if ColumnName.CLUSTER not in df.columns:
            raise ValueError(f"Column '{ColumnName.CLUSTER}' not found in file")
        
        # Extract unique cluster numbers and sort
        clusters = sorted(df[ColumnName.CLUSTER].dropna().unique().astype(int).tolist())
        
        if not clusters:
            raise ValueError("No valid cluster numbers found in file")
        
        return {"clusters": clusters}
    
    except ValueError as ve:
        raise HTTPException(status_code=400, detail=str(ve))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to read clusters: {str(e)}")
