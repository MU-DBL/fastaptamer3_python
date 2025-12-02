import mimetypes
from fastapi import APIRouter, File, UploadFile, HTTPException
from pathlib import Path
import os
from fastapi.responses import FileResponse
import shutil
from datetime import datetime
from pydantic import BaseModel
from services import recount_service

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))

class RecountInput(BaseModel):
    input_path_1: str="" 
    input_path_2: str=""
    scaling_factor: float = 1e6
    output_format: str=""

@router.post("/recount")
async def recount(params:RecountInput):
    filepath_1 = f"{UPLOAD_DIR}/{params.input_path_1}"
    filepath_2 = f"{UPLOAD_DIR}/{params.input_path_2}"
    base_name_1 = os.path.splitext(os.path.basename(params.input_path_1))[0]
    base_name_2 = os.path.splitext(os.path.basename(params.input_path_2))[0]
    output_format=params.output_format
    output_path=f"{UPLOAD_DIR}/combine_{base_name_1}_{base_name_2}.{output_format}"
    
    output_path = recount_service.run_recount(
        input_path_1=filepath_1,
        input_path_2=filepath_2,
        output_path=output_path,
        output_format=output_format,
        scaling_factor=params.scaling_factor
    )

    return {"status": "ok", "result": os.path.basename(output_path)}

