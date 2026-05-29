import asyncio
import os
from pathlib import Path
import queue
import uuid
from fastapi import APIRouter, BackgroundTasks
from pydantic import BaseModel
from services import preprocess_service
import json
import state
from state import progress_queues

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))

class PreprocessInput(BaseModel):
    input_path: str | None = None
    const5p: str = ""
    const3p: str = ""
    trim5_fixed: int = 0
    trim3_fixed: int = 0
    min_length: int = 0
    max_length: int = 100
    max_error: float = 0.005
    adapter_error_rate: float = 0.1
    output_format: str = "fasta"

@router.post("/preprocess")
async def preprocess(params: PreprocessInput, background_tasks: BackgroundTasks):
    filepath = f"{UPLOAD_DIR}/{params.input_path}"
    p = Path(params.input_path)
    base_name = p.stem if p.suffix.lower() != '.gz' else Path(p.stem).stem
    output_format=params.output_format
    output_path=f"{UPLOAD_DIR}/{base_name}_preprocess.{output_format}"

    job_id = str(uuid.uuid4())  # Simple job ID generation

    async def delayed_start():
        await asyncio.sleep(0.5)  # Give client time to connect
        await preprocess_service.run_preprocess(
            job_id,
            filepath,
            const5p=params.const5p,
            const3p=params.const3p,
            trim5_fixed=params.trim5_fixed,
            trim3_fixed=params.trim3_fixed,
            length_range=[params.min_length, params.max_length],
            max_error=params.max_error,
            adapter_error_rate=params.adapter_error_rate,
            output_path=output_path,
            output_format=output_format
        )
    
    state.current_task = asyncio.create_task(delayed_start())
    # background_tasks.add_task(delayed_start)
    
    return {"status": "ok", "result": job_id, "output_path": os.path.basename(output_path)}

