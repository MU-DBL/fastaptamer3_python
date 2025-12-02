import mimetypes
from fastapi import APIRouter, File, UploadFile, HTTPException
from pathlib import Path
import os
from fastapi.responses import FileResponse
from datetime import datetime
from pydantic import BaseModel
from services import count_service

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))

class CountInput(BaseModel):
    input_path: str=""
    reverseComplement: bool=False
    scaling_factor: float = 1e6
    output_format: str=""

@router.post("/count")
async def count(params:CountInput):
    filepath = f"{UPLOAD_DIR}/{params.input_path}"
    base_name = os.path.splitext(os.path.basename(params.input_path))[0]
    output_format=params.output_format
    output_path=f"{UPLOAD_DIR}/{base_name}_count.{output_format}"

    output_path = count_service.run_count(
        inputpath=filepath,  
        reverseComplement=params.reverseComplement, 
        scaling_factor=params.scaling_factor, 
        output_format=output_format, 
        output_path=output_path)
    
    return {"status": "ok", "result": os.path.basename(output_path)}
