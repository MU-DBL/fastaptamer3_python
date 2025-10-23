import os
from pathlib import Path
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Optional
from services import preprocess_service

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))

class PreprocessInput(BaseModel):
    input_path: str | None = None
    const5p: str = ""
    const3p: str = ""
    min_length: int = 0
    max_length: int = 100
    max_error: float = 0.005

@router.post("/preprocess")
async def preprocess(params: PreprocessInput):
    filepath = f"{UPLOAD_DIR}/{params.input_path}"
    base_name = os.path.splitext(os.path.basename(params.input_path))[0]
    output_format='fasta'
    output_path=f"{UPLOAD_DIR}/{base_name}_preprocess.{output_format}"

    output_path = preprocess_service.run_preprocess(
        input_path=filepath,
        const5p=params.const5p,
        const3p=params.const3p,
        length_range=[params.min_length, params.max_length], 
        max_error=params.max_error,
        output_path=output_path,
        output_format=output_format
    )

    return {"status": "ok", "result": os.path.basename(output_path)}
