import os
from pathlib import Path
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from services import motif_service

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))


class MotifSearchInput(BaseModel):
    input_path: str = ""
    motif: str = ""
    highlight: bool = False
    partial: bool = False
    motif_type: str = "Nucleotide"
    max_mismatches: int = 0
    max_mismatches_per_motif: str = ""
    output_format: str = "fasta"


@router.post("/motif-search")
async def motif_search(params: MotifSearchInput):
    """
    Search for sequences containing user-defined motifs.

    Args:
        params: MotifSearchInput containing:
            - input_path: Name of the input FASTA file
            - motif: Comma-separated list of motifs
            - highlight: Whether to highlight motifs in output
            - partial: Whether to allow partial matches (OR vs AND)
            - motif_type: Type of motif (Nucleotide, AminoAcid, String)
            - max_mismatches: Mismatch tolerance applied to every motif
            - max_mismatches_per_motif: Optional comma-separated per-motif mismatch tolerances,
              same order and count as `motif`; overrides max_mismatches when provided
            - output_format: Output format (fasta or csv)

    Returns:
        JSON response with status and result filename
    """
    filepath = f"{UPLOAD_DIR}/{params.input_path}"
    base_name = os.path.splitext(os.path.basename(params.input_path))[0]
    output_format = params.output_format
    output_path = f"{UPLOAD_DIR}/{base_name}_motif_search.{output_format}"

    try:
        output_path = motif_service.search_motif(
            fasta_input=filepath,
            motif=params.motif,
            highlight=params.highlight,
            partial=params.partial,
            motif_type=params.motif_type,
            max_mismatches=params.max_mismatches,
            max_mismatches_per_motif=params.max_mismatches_per_motif,
            output_format=output_format,
            output_path=output_path
        )
    except ValueError as ve:
        raise HTTPException(status_code=400, detail=str(ve))

    return {"status": "ok", "result": os.path.basename(output_path)}


class MotifOmitInput(BaseModel):
    input_path: str = ""
    motif: str = ""
    partial: bool = False
    motif_type: str = "Nucleotide"
    max_mismatches: int = 0
    max_mismatches_per_motif: str = ""
    output_format: str = "fasta"


@router.post("/motif-omit")
async def motif_omit(params: MotifOmitInput):
    """
    Omit sequences containing user-defined motifs.

    Args:
        params: MotifOmitInput containing:
            - input_path: Name of the input FASTA file
            - motif: Comma-separated list of motifs
            - partial: Whether to allow partial matches (OR vs AND)
            - motif_type: Type of motif (Nucleotide, AminoAcid, String)
            - max_mismatches: Mismatch tolerance applied to every motif
            - max_mismatches_per_motif: Optional comma-separated per-motif mismatch tolerances,
              same order and count as `motif`; overrides max_mismatches when provided
            - output_format: Output format (fasta or csv)

    Returns:
        JSON response with status and result filename
    """
    filepath = f"{UPLOAD_DIR}/{params.input_path}"
    base_name = os.path.splitext(os.path.basename(params.input_path))[0]
    output_format = params.output_format
    output_path = f"{UPLOAD_DIR}/{base_name}_motif_omit.{output_format}"

    try:
        output_path = motif_service.omit_motif(
            fasta_input=filepath,
            motif=params.motif,
            partial=params.partial,
            motif_type=params.motif_type,
            max_mismatches=params.max_mismatches,
            max_mismatches_per_motif=params.max_mismatches_per_motif,
            output_format=output_format,
            output_path=output_path
        )
    except ValueError as ve:
        raise HTTPException(status_code=400, detail=str(ve))

    return {"status": "ok", "result": os.path.basename(output_path)}
