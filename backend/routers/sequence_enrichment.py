from ast import Dict
from services.file_service import read_file, save_sequences
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from pathlib import Path
import os
import pandas as pd
from services.constants import ColumnName
from numba import jit, prange
from typing import Any, Optional, List, Tuple
import numpy as np

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))

class EnrichInput(BaseModel):
    fadf1_cluster_path: str = ""
    fadf2_cluster_path: str = ""
    keep_na: bool = False
    output_format: str = "csv"

class EnrichResponse(BaseModel):
    status: str
    result: str
    num_sequences: int
    enrichment_stats: Dict[str, Any]

@router.post("/sequence-enrich", response_model=EnrichResponse)
async def fa_enrich_endpoint(params: EnrichInput):
    """
    Calculate enrichment between two populations by comparing RPU values.
    Returns merged data with Enrichment and log2(Enrichment) columns.
    """
    if not params.fadf1_cluster_path or not params.fadf2_cluster_path:
        raise HTTPException(
            status_code=400,
            detail="Both fadf1_cluster_path and fadf2_cluster_path are required"
        )
    
    filepath1 = UPLOAD_DIR / params.fadf1_cluster_path
    filepath2 = UPLOAD_DIR / params.fadf2_cluster_path
    output_path = UPLOAD_DIR / f"enriched_{'full' if params.keep_na else 'inner'}.{params.output_format}"
    
    try:
        # Read clustered data from both populations
        fadf1 = read_file(filepath1)
        fadf2 = read_file(filepath2)
        
        # Perform enrichment analysis
        result_df = fa_enrich(
            fadf1=fadf1,
            fadf2=fadf2,
            keep_na=params.keep_na
        )
        
        # Save output
        save_sequences(result_df, str(output_path), params.output_format)
        
        # Calculate enrichment statistics
        enrichment_stats = calculate_enrichment_stats(result_df)
        
        return EnrichResponse(
            status="ok",
            result=output_path.name,
            num_sequences=len(result_df),
            enrichment_stats=enrichment_stats
        )
    
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Enrichment analysis failed: {str(e)}"
        )


def fa_enrich(
    fadf1: pd.DataFrame,
    fadf2: pd.DataFrame,
    keep_na: bool = False
) -> pd.DataFrame:
    """
    Calculate enrichment between two populations.
    
    Args:
        fadf1: First population dataframe (baseline)
        fadf2: Second population dataframe (enriched)
        keep_na: If True, keeps sequences present in only one population (full join)
                 If False, only keeps sequences present in both (inner join)
    
    Returns:
        DataFrame with enrichment metrics
    """
    # Rename columns for population 1 (add .a suffix)
    population1_rename = {
        ColumnName.ID: "ID.a",
        ColumnName.RANK: "Rank.a",
        ColumnName.READS: "Reads.a",
        ColumnName.RPU: "RPU.a",
        ColumnName.CLUSTER: "Cluster.a",
        ColumnName.RANK_IN_CLUSTER: "RankInCluster.a",
        ColumnName.LED: "LED.a"
    }
    fadf1_renamed = fadf1.rename(columns={
        k: v for k, v in population1_rename.items() if k in fadf1.columns
    })
    
    # Rename columns for population 2 (add .b suffix)
    population2_rename = {
        ColumnName.ID: "ID.b",
        ColumnName.RANK: "Rank.b",
        ColumnName.READS: "Reads.b",
        ColumnName.RPU: "RPU.b",
        ColumnName.CLUSTER: "Cluster.b",
        ColumnName.RANK_IN_CLUSTER: "RankInCluster.b",
        ColumnName.LED: "LED.b"
    }
    fadf2_renamed = fadf2.rename(columns={
        k: v for k, v in population2_rename.items() if k in fadf2.columns
    })
    
    # Merge dataframes
    merge_how = "outer" if keep_na else "inner"
    merge_df = pd.merge(
        fadf1_renamed,
        fadf2_renamed,
        on=ColumnName.SEQUENCES,
        how=merge_how
    )
    
    # Calculate enrichment and log2(enrichment)
    # Handle division by zero and NaN values
    rpu_a = merge_df["RPU.a"].fillna(0)
    rpu_b = merge_df["RPU.b"].fillna(0)
    
    # Calculate enrichment (avoid division by zero)
    enrichment = np.where(
        rpu_a > 0,
        rpu_b / rpu_a,
        np.where(rpu_b > 0, np.inf, 0)  # inf if only in pop2, 0 if both are 0
    )
    
    # Calculate log2 enrichment (handle inf and 0)
    log2_enrichment = np.where(
        enrichment > 0,
        np.log2(enrichment),
        np.where(enrichment == 0, -np.inf, np.inf)
    )
    
    merge_df["Enrichment"] = np.round(enrichment, 3)
    merge_df["log2E"] = np.round(log2_enrichment, 3)
    
    # Replace inf with a large number for practical purposes (optional)
    merge_df["Enrichment"] = merge_df["Enrichment"].replace([np.inf, -np.inf], [999.999, -999.999])
    merge_df["log2E"] = merge_df["log2E"].replace([np.inf, -np.inf], [20.0, -20.0])
    
    # Fill remaining NaN values with 0
    merge_df = merge_df.fillna(0)
    
    # Move Sequences column to first position and sort by Rank.a
    cols = [ColumnName.SEQUENCES] + [col for col in merge_df.columns if col != ColumnName.SEQUENCES]
    merge_df = merge_df[cols]
    
    # Sort by Rank.a (sequences not in pop1 will have Rank.a = 0, so they go to top)
    merge_df = merge_df.sort_values("Rank.a", na_position='last')
    
    return merge_df


def calculate_enrichment_stats(df: pd.DataFrame) -> Dict[str, Any]:
    """Calculate summary statistics for enrichment analysis."""
    # Filter out inf values for statistics
    enrichment_values = df["Enrichment"].replace([np.inf, -np.inf], np.nan).dropna()
    log2e_values = df["log2E"].replace([np.inf, -np.inf], np.nan).dropna()
    
    stats = {
        "enrichment": {
            "min": float(enrichment_values.min()) if len(enrichment_values) > 0 else 0,
            "max": float(enrichment_values.max()) if len(enrichment_values) > 0 else 0,
            "mean": float(enrichment_values.mean()) if len(enrichment_values) > 0 else 0,
            "median": float(enrichment_values.median()) if len(enrichment_values) > 0 else 0,
            "std": float(enrichment_values.std()) if len(enrichment_values) > 0 else 0
        },
        "log2_enrichment": {
            "min": float(log2e_values.min()) if len(log2e_values) > 0 else 0,
            "max": float(log2e_values.max()) if len(log2e_values) > 0 else 0,
            "mean": float(log2e_values.mean()) if len(log2e_values) > 0 else 0,
            "median": float(log2e_values.median()) if len(log2e_values) > 0 else 0,
            "std": float(log2e_values.std()) if len(log2e_values) > 0 else 0
        },
        "sequences_only_in_pop1": int((df["RPU.b"] == 0).sum()),
        "sequences_only_in_pop2": int((df["RPU.a"] == 0).sum()),
        "sequences_in_both": int(((df["RPU.a"] > 0) & (df["RPU.b"] > 0)).sum()),
        "highly_enriched_count": int((df["Enrichment"] > 10).sum()),
        "highly_depleted_count": int((df["Enrichment"] < 0.1).sum())
    }
    
    return stats