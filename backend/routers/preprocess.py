import asyncio
import os
from pathlib import Path
import queue
import uuid
from fastapi import APIRouter, BackgroundTasks
from pydantic import BaseModel
from services import preprocess_service
import json
from state import progress_queues

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))

class PreprocessInput(BaseModel):
    input_path: str | None = None
    const5p: str = ""
    const3p: str = ""
    min_length: int = 0
    max_length: int = 100
    max_error: float = 0.005
    output_format: str = "fasta"

@router.post("/preprocess")
async def preprocess(params: PreprocessInput, background_tasks: BackgroundTasks):
    filepath = f"{UPLOAD_DIR}/{params.input_path}"
    base_name = os.path.splitext(os.path.basename(params.input_path))[0]
    output_format=params.output_format
    output_path=f"{UPLOAD_DIR}/{base_name}_preprocess.{output_format}"

    job_id = str(uuid.uuid4())  # Simple job ID generation
    progress_queues[job_id] = queue.Queue()

    async def delayed_start():
        await asyncio.sleep(0.5)  # Give client time to connect
        preprocess_service.run_preprocess(
            job_id,
            filepath,
            const5p=params.const5p,
            const3p=params.const3p,
            length_range=[params.min_length, params.max_length],
            max_error=params.max_error,
            output_path=output_path,
            output_format=output_format
        )
    
    background_tasks.add_task(delayed_start)
    
    return {"status": "ok", "result": job_id, "output_path": output_path}
