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
    input_paths: list[str]
    scaling_factor: float = 1e6
    output_format: str=""

@router.post("/recount")
async def recount(params: RecountInput):
    if len(params.input_paths) < 2 or len(params.input_paths) > 5:
        raise HTTPException(
            status_code=400,
            detail="Recount requires between 2 and 5 input files."
        )
    filepaths = [f"{UPLOAD_DIR}/{p}" for p in params.input_paths]
    base_names = [os.path.splitext(os.path.basename(p))[0] for p in params.input_paths]
    output_format = params.output_format
    combined_name = "_".join(base_names)
    output_path = f"{UPLOAD_DIR}/combine_{combined_name}.{output_format}"

    output_path = recount_service.run_recount(
        input_paths=filepaths,
        output_path=output_path,
        output_format=output_format,
        scaling_factor=params.scaling_factor
    )

    return {"status": "ok", "result": os.path.basename(output_path)}

