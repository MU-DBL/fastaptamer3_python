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

        # Try candidate cluster columns in priority order
        cluster_column = None
        for candidate in [ColumnName.CLUSTER, "Cluster.a"]:
            if candidate in df.columns:
                cluster_column = candidate
                break

        if cluster_column is None:
            # No cluster column found — caller will treat all sequences as one group
            return {"clusters": [], "cluster_column": None}

        # Extract unique cluster numbers and sort
        clusters = sorted(df[cluster_column].dropna().unique().astype(int).tolist())

        return {"clusters": clusters, "cluster_column": cluster_column}

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to read clusters ({type(e).__name__}): {str(e)}")
