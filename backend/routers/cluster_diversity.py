from fastapi import APIRouter, HTTPException
from pathlib import Path
import os
from pydantic import BaseModel
from typing import List
from services import cluster_service
from services.kmer_service import analyze_kmer_diversity

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))

class DiversityInput(BaseModel):
    input_path: str = ""
    output_format: str = "csv"

class KmerAnalysisInput(BaseModel):
    input_path: str
    selected_clusters: List[int]
    k_size: int = 4
    method: str = "pca"  # "pca" or "umap"

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

@router.post("/cluster-kmer-analysis")
async def cluster_kmer_analysis(params: KmerAnalysisInput):
    """
    Perform k-mer analysis on selected clusters.
    
    This endpoint implements the R workflow from clusterDiversity_plots.R:
    1. Filters sequences by selected clusters
    2. Generates k-mer count matrix (equivalent to kmer::kcount)
    3. Applies PCA or UMAP dimensionality reduction
    4. Returns coordinates with cluster assignments and colors
    """
    if not params.input_path:
        raise HTTPException(status_code=400, detail="input_path is required")
    
    if not params.selected_clusters:
        raise HTTPException(status_code=400, detail="selected_clusters cannot be empty")
    
    if params.method not in ["pca", "umap"]:
        raise HTTPException(status_code=400, detail="method must be 'pca' or 'umap'")
    
    if params.k_size < 2 or params.k_size > 8:
        raise HTTPException(status_code=400, detail="k_size must be between 2 and 8")
    
    filepath = UPLOAD_DIR / params.input_path
    
    if not filepath.exists():
        raise HTTPException(status_code=404, detail=f"File not found: {params.input_path}")
    
    try:
        # Perform k-mer analysis
        result = analyze_kmer_diversity(
            input_path=str(filepath),
            selected_clusters=params.selected_clusters,
            k_size=params.k_size,
            method=params.method
        )
        
        return {
            "status": "ok",
            "data": result
        }
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"K-mer analysis failed: {str(e)}")
