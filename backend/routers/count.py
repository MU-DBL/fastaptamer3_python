import asyncio
import mimetypes
from fastapi import APIRouter, File, UploadFile, HTTPException
from pathlib import Path
import os
from fastapi.responses import FileResponse
from datetime import datetime
from pydantic import BaseModel
from typing import List, Optional
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
    p = Path(params.input_path)
    base_name = p.stem if p.suffix.lower() != '.gz' else Path(p.stem).stem
    output_format=params.output_format
    output_path=f"{UPLOAD_DIR}/{base_name}_count.{output_format}"

    loop = asyncio.get_event_loop()
    output_path = await loop.run_in_executor(
        None,
        lambda: count_service.run_count(
            inputpath=filepath,
            reverseComplement=params.reverseComplement,
            scaling_factor=params.scaling_factor,
            output_format=output_format,
            output_path=output_path,
        )
    )

    return {"status": "ok", "result": os.path.basename(output_path)}


class CountPreviewInput(BaseModel):
    input_path: str = ""
    limit: int = 50000


@router.post("/count-preview")
async def count_preview(params: CountPreviewInput):
    filepath = UPLOAD_DIR / params.input_path
    if not filepath.exists():
        raise HTTPException(status_code=404, detail=f"File not found: {params.input_path}")

    loop = asyncio.get_event_loop()
    result = await loop.run_in_executor(
        None, lambda: count_service.get_table_preview(str(filepath), params.limit)
    )
    return {"status": "ok", **result}


class ReadsPerRankInput(BaseModel):
    input_path: str = ""
    min_reads: int = 0
    max_rank: Optional[int] = None
    metric: str = "reads"  # "reads" or "rpu"


@router.post("/count-reads-per-rank")
async def count_reads_per_rank(params: ReadsPerRankInput):
    filepath = UPLOAD_DIR / params.input_path
    if not filepath.exists():
        raise HTTPException(status_code=404, detail=f"File not found: {params.input_path}")

    loop = asyncio.get_event_loop()
    points = await loop.run_in_executor(
        None,
        lambda: count_service.get_reads_per_rank(
            str(filepath), min_reads=params.min_reads, max_rank=params.max_rank, metric=params.metric
        )
    )
    return {"status": "ok", "data": [{"rank": r, "value": v} for r, v in points]}


class SequenceLengthHistogramInput(BaseModel):
    input_path: str = ""


@router.post("/count-sequence-length-histogram")
async def count_sequence_length_histogram(params: SequenceLengthHistogramInput):
    filepath = UPLOAD_DIR / params.input_path
    if not filepath.exists():
        raise HTTPException(status_code=404, detail=f"File not found: {params.input_path}")

    loop = asyncio.get_event_loop()
    histogram = await loop.run_in_executor(
        None, lambda: count_service.get_sequence_length_histogram(str(filepath))
    )
    return {"status": "ok", **histogram}


class AbundanceBinsInput(BaseModel):
    input_path: str = ""
    breakpoints: List[int] = [10, 100, 1000]
    use_singleton: bool = True


@router.post("/count-abundance")
async def count_abundance(params: AbundanceBinsInput):
    filepath = UPLOAD_DIR / params.input_path
    if not filepath.exists():
        raise HTTPException(status_code=404, detail=f"File not found: {params.input_path}")

    loop = asyncio.get_event_loop()
    bins = await loop.run_in_executor(
        None,
        lambda: count_service.get_abundance_bins(
            str(filepath), breakpoints=params.breakpoints, use_singleton=params.use_singleton
        )
    )
    return {"status": "ok", "data": bins}
