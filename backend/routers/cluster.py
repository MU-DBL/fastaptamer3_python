import mimetypes
from fastapi import APIRouter, File, UploadFile, HTTPException
from pathlib import Path
import os
from fastapi.responses import FileResponse
from datetime import datetime
from pydantic import BaseModel
from services import cluster_service

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))

class ClusterInput(BaseModel):
    input_path: str = ""
    output_format: str = "fasta"
    min_reads: int = 10
    max_led: int = 7
    total_clusters: int = 30
    keep_nc: bool = True

@router.post("/clusterled")
async def cluster(params:ClusterInput):
    if not params.input_path:
        raise HTTPException(status_code=400, detail="input_path is required")
     
    filepath = f"{UPLOAD_DIR}/{params.input_path}"
    base_name = os.path.splitext(os.path.basename(params.input_path))[0]
    output_format=params.output_format
    output_path=f"{UPLOAD_DIR}/{base_name}_cluster.{output_format}"

    try:

        output_path = cluster_service.run_cluster_led(
                input_path=str(filepath),
                output_format=output_format,
                output_path=str(output_path),
                min_reads=params.min_reads,
                max_led=params.max_led,
                total_clusters=params.total_clusters,
                keep_nc=params.keep_nc
            )
        
        return {"status": "ok", "result": os.path.basename(output_path)}

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Clustering failed: {str(e)}")
