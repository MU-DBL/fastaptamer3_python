from fastapi import APIRouter, HTTPException
from pathlib import Path
import os
from pydantic import BaseModel
from services import cluster_service

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))

class DiversityInput(BaseModel):
    input_path: str = ""
    output_format: str = "csv"

@router.post("/clusterdiversity")
async def cluster_diversity(params: DiversityInput):
    if not params.input_path:
        raise HTTPException(status_code=400, detail="input_path is required")
    
    filepath = UPLOAD_DIR / params.input_path
    try:
        # Generate output path
        base_name = os.path.splitext(os.path.basename(params.input_path))[0]
        output_format = params.output_format.lower()
        output_path = UPLOAD_DIR / f"{base_name}_diversity.{output_format}"
        
        # Run diversity analysis
        output_path = cluster_service.run_cluster_diversity(
            inputpath=str(filepath),
            output_format=output_format,
            output_path=str(output_path)
        )
        
        return {
            "status": "ok", 
            "result": os.path.basename(output_path)
        }
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Diversity analysis failed: {str(e)}")
