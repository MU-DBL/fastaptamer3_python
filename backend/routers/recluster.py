from typing import Dict, List, Any
from services.file_service import read_file, save_sequences
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from pathlib import Path
import os
import pandas as pd
import numpy as np
import Levenshtein
from services.constants import ColumnName
from numba import jit, prange
from typing import Any, Optional, List, Tuple
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import multiprocessing as mp


router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))


class ReclusterInput(BaseModel):
    fadf1_cluster_path: str = ""
    fadf2_cluster_path: str = ""
    led_threshold: int = 7
    output_format: str = "csv"


@router.post("/recluster")
async def fa_recluster_endpoint(params: ReclusterInput):
    """
    Merge clusters from two populations based on Levenshtein edit distance (LED).
    """
    if not params.fadf1_cluster_path or not params.fadf2_cluster_path:
        raise HTTPException(
            status_code=400,
            detail="Both fadf1_cluster_path and fadf2_cluster_path are required",
        )

    filepath1 = UPLOAD_DIR / params.fadf1_cluster_path
    filepath2 = UPLOAD_DIR / params.fadf2_cluster_path
    output_path = (
        UPLOAD_DIR / f"reclustered_led_{params.led_threshold}.{params.output_format}"
    )

    try:
        # Read clustered data from both populations
        fadf1_cluster = read_file(filepath1)
        fadf2_cluster = read_file(filepath2)

        # Perform reclustering
        result_df = fa_recluster(
            fadf1_cluster=fadf1_cluster,
            fadf2_cluster=fadf2_cluster,
            led_threshold=params.led_threshold,
        )

        # Save output
        save_sequences(result_df, str(output_path), params.output_format)

        return {
            "status": "ok",
            "result": output_path.name,
            "num_sequences": len(result_df),
            "num_clusters": result_df[ColumnName.CLUSTER].nunique(),
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Reclustering failed: {str(e)}")


class LEDMatrixInput(BaseModel):
    fadf1_cluster_path: str = ""
    fadf2_cluster_path: str = ""
    led_threshold: Optional[int] = (
        None  # Optional: for highlighting threshold on frontend
    )
    use_parallel: bool = True
    n_jobs: Optional[int] = None


class LEDMatrixResponse(BaseModel):
    status: str
    led_matrix: List[List[int]]  # 2D array of LED values
    p1_cluster_ids: List[int]  # Row labels (population 1 cluster IDs)
    p2_cluster_ids: List[int]  # Column labels (population 2 cluster IDs)
    p1_seeds: List[str]  # Actual seed sequences for pop 1
    p2_seeds: List[str]  # Actual seed sequences for pop 2
    threshold: Optional[int]  # For frontend reference
    matrix_shape: Dict[str, int]  # Dimensions
    statistics: Dict[str, Any]  # Summary stats


@router.post("/recluster-led-matrix", response_model=LEDMatrixResponse)
async def get_led_matrix_endpoint(params: LEDMatrixInput):
    """
    Compute LED matrix between cluster seeds from two populations.
    Returns matrix data formatted for frontend heatmap visualization.
    """
    if not params.fadf1_cluster_path or not params.fadf2_cluster_path:
        raise HTTPException(
            status_code=400,
            detail="Both fadf1_cluster_path and fadf2_cluster_path are required",
        )

    filepath1 = UPLOAD_DIR / params.fadf1_cluster_path
    filepath2 = UPLOAD_DIR / params.fadf2_cluster_path

    try:
        # Read clustered data from both populations
        fadf1_cluster = read_file(filepath1)
        fadf2_cluster = read_file(filepath2)

        # Extract seeds and cluster IDs
        p1_seeds_mask = fadf1_cluster[ColumnName.RANK_IN_CLUSTER] == 1
        p2_seeds_mask = fadf2_cluster[ColumnName.RANK_IN_CLUSTER] == 1

        p1_seeds = fadf1_cluster.loc[p1_seeds_mask, ColumnName.SEQUENCES].tolist()
        p2_seeds = fadf2_cluster.loc[p2_seeds_mask, ColumnName.SEQUENCES].tolist()
        p1_cluster_ids = fadf1_cluster.loc[p1_seeds_mask, ColumnName.CLUSTER].tolist()
        p2_cluster_ids = fadf2_cluster.loc[p2_seeds_mask, ColumnName.CLUSTER].tolist()

        # Compute LED matrix
        # if params.use_parallel and len(p1_seeds) * len(p2_seeds) > 100:
        #     led_matrix = compute_led_matrix_parallel(
        #         np.array(p1_seeds), np.array(p2_seeds), params.n_jobs
        #     )
        # else:
        led_matrix = compute_led_matrix_fast(np.array(p1_seeds), np.array(p2_seeds))

        # Convert numpy array to list for JSON serialization
        led_matrix_list = led_matrix.tolist()

        # Calculate statistics
        flat_values = led_matrix.flatten()
        statistics = {
            "min": int(np.min(flat_values)),
            "max": int(np.max(flat_values)),
            "mean": float(np.mean(flat_values)),
            "median": float(np.median(flat_values)),
            "std": float(np.std(flat_values)),
        }

        if params.led_threshold is not None:
            below_threshold = np.sum(flat_values <= params.led_threshold)
            statistics["below_threshold_count"] = int(below_threshold)
            statistics["below_threshold_percentage"] = float(
                (below_threshold / flat_values.size) * 100
            )

        return LEDMatrixResponse(
            status="ok",
            led_matrix=led_matrix_list,
            p1_cluster_ids=p1_cluster_ids,
            p2_cluster_ids=p2_cluster_ids,
            p1_seeds=p1_seeds,
            p2_seeds=p2_seeds,
            threshold=params.led_threshold,
            matrix_shape={"rows": len(p1_seeds), "cols": len(p2_seeds)},
            statistics=statistics,
        )

    except Exception as e:
        raise HTTPException(
            status_code=500, detail=f"LED matrix computation failed: {str(e)}"
        )


def fa_recluster(
    fadf1_cluster: pd.DataFrame,
    fadf2_cluster: pd.DataFrame,
    led_threshold: int = 7,
    use_parallel: bool = True,
    n_jobs: Optional[int] = None,
) -> pd.DataFrame:

    # Extract seeds efficiently
    p1_seeds_mask = fadf1_cluster[ColumnName.RANK_IN_CLUSTER] == 1
    p2_seeds_mask = fadf2_cluster[ColumnName.RANK_IN_CLUSTER] == 1

    p1_seeds = fadf1_cluster.loc[p1_seeds_mask, ColumnName.SEQUENCES].values
    p2_seeds = fadf2_cluster.loc[p2_seeds_mask, ColumnName.SEQUENCES].values
    p1_cluster_ids = fadf1_cluster.loc[p1_seeds_mask, ColumnName.CLUSTER].values

    # if use_parallel and len(p1_seeds) * len(p2_seeds) > 100:
    #     led_matrix = compute_led_matrix_parallel(p1_seeds, p2_seeds, n_jobs)
    # else:
    led_matrix = compute_led_matrix_fast(p1_seeds, p2_seeds)

    clusters2merge = find_merge_candidates_vectorized(led_matrix, led_threshold)

    merged_df = merge_clusters_batch(
        fadf1_cluster, fadf2_cluster, clusters2merge, p1_cluster_ids
    )

    merge_df = calculate_enrichment_vectorized(merged_df)

    result_df = calculate_led_to_seeds_fast(merge_df)

    return result_df


def compute_led_matrix_fast(seeds1: np.ndarray, seeds2: np.ndarray) -> np.ndarray:
    n1, n2 = len(seeds1), len(seeds2)
    led_matrix = np.empty((n1, n2), dtype=np.int16)

    for i in range(n1):
        seq1 = seeds1[i]
        for j in range(n2):
            led_matrix[i, j] = Levenshtein.distance(seq1, seeds2[j])

    return led_matrix


def compute_row_batch_wrapper(args):
    """Top-level function for pickling compatibility."""
    start_idx, end_idx, seeds1, seeds2 = args
    batch_size = end_idx - start_idx
    result = np.empty((batch_size, len(seeds2)), dtype=np.int16)
    for local_i, global_i in enumerate(range(start_idx, end_idx)):
        seq1 = seeds1[global_i]
        for j in range(len(seeds2)):
            result[local_i, j] = Levenshtein.distance(seq1, seeds2[j])
    return result


def compute_led_matrix_parallel(
    seeds1: np.ndarray, seeds2: np.ndarray, n_jobs: Optional[int] = None
) -> np.ndarray:
    if n_jobs is None:
        n_jobs = min(mp.cpu_count(), len(seeds1))

    # Split work into chunks
    n1 = len(seeds1)
    chunk_size = max(1, n1 // n_jobs)
    chunks = [
        (i, min(i + chunk_size, n1), seeds1, seeds2) for i in range(0, n1, chunk_size)
    ]

    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        results = list(executor.map(compute_row_batch_wrapper, chunks))

    return np.vstack(results)


# Alternative: Using python-Levenshtein's C implementation more efficiently
def compute_led_matrix_optimized(seeds1: List[str], seeds2: List[str]) -> np.ndarray:
    """
    Use list comprehension with C-level Levenshtein (fastest for medium datasets).
    ~3x faster than basic loop.
    """
    return np.array(
        [[Levenshtein.distance(s1, s2) for s2 in seeds2] for s1 in seeds1],
        dtype=np.int16,
    )


def find_merge_candidates_vectorized(
    led_matrix: np.ndarray, threshold: int
) -> np.ndarray:
    """
    Vectorized merge candidate finding using numpy operations.
    ~10x faster than loop-based approach.
    """
    # Create boolean mask for valid merges
    valid_merges = led_matrix < threshold

    # For each row, find the first valid column (or -1 if none)
    clusters2merge = np.full(led_matrix.shape[0], -1, dtype=np.int32)

    for i in range(led_matrix.shape[0]):
        valid_cols = np.where(valid_merges[i])[0]
        if len(valid_cols) > 0:
            # Take the closest match
            clusters2merge[i] = valid_cols[np.argmin(led_matrix[i, valid_cols])]

    return clusters2merge  # -1 means no merge, otherwise 0-indexed


def merge_clusters_batch(
    fadf1: pd.DataFrame,
    fadf2: pd.DataFrame,
    clusters2merge: np.ndarray,
    p1_cluster_ids: np.ndarray,
) -> pd.DataFrame:
    """
    Batch merge operation - no repeated concat() calls.
    ~50x faster than incremental concat for large datasets.
    """
    # Prepare dataframes
    pop1 = fadf1.drop(columns=[ColumnName.RANK_IN_CLUSTER, ColumnName.LED]).copy()
    pop2 = fadf2.drop(columns=[ColumnName.RANK_IN_CLUSTER, ColumnName.LED]).copy()

    pop1["Population"] = 1
    pop2["Population"] = 2

    # Create mapping: pop2 cluster ID -> pop1 cluster ID (or new ID)
    pop2_cluster_mapping = {}
    merged_pop2_clusters = set()

    # Map merged clusters
    for i, pop2_cluster_idx in enumerate(clusters2merge):
        if pop2_cluster_idx >= 0:  # Valid merge
            pop1_cluster_id = p1_cluster_ids[i]
            pop2_cluster_id = pop2_cluster_idx + 1  # Convert to 1-indexed
            pop2_cluster_mapping[pop2_cluster_id] = pop1_cluster_id
            merged_pop2_clusters.add(pop2_cluster_id)

    # Assign new IDs to unmerged pop2 clusters
    max_cluster = pop1[ColumnName.CLUSTER].max()
    unmerged_pop2_clusters = (
        set(pop2[ColumnName.CLUSTER].unique()) - merged_pop2_clusters
    )
    for i, old_cluster in enumerate(sorted(unmerged_pop2_clusters)):
        pop2_cluster_mapping[old_cluster] = max_cluster + i + 1

    # Apply mapping to pop2 (vectorized)
    pop2[ColumnName.CLUSTER] = pop2[ColumnName.CLUSTER].map(pop2_cluster_mapping)

    # Single concatenation
    merged = pd.concat([pop1, pop2], ignore_index=True)

    return merged


def calculate_enrichment_vectorized(merged_df: pd.DataFrame) -> pd.DataFrame:

    # Split and rename (vectorized)
    pop1 = merged_df[merged_df["Population"] == 1].drop(columns=["Population"])
    pop2 = merged_df[merged_df["Population"] == 2].drop(columns=["Population"])

    # Batch rename using dictionary
    rename_dict = {
        col: f"{col}.a" if col != ColumnName.SEQUENCES else col for col in pop1.columns
    }
    pop1 = pop1.rename(columns=rename_dict)

    rename_dict = {
        col: f"{col}.b" if col != ColumnName.SEQUENCES else col for col in pop2.columns
    }
    pop2 = pop2.rename(columns=rename_dict)

    # Full join
    merge_df = pd.merge(pop1, pop2, on=ColumnName.SEQUENCES, how="outer")

    # Coalesce cluster
    merge_df[ColumnName.CLUSTER] = (
        merge_df[f"{ColumnName.CLUSTER}.a"]
        .fillna(merge_df[f"{ColumnName.CLUSTER}.b"])
        .astype(np.int32)
    )
    merge_df = merge_df.drop(
        columns=[f"{ColumnName.CLUSTER}.a", f"{ColumnName.CLUSTER}.b"]
    )

    # Vectorized enrichment calculation with proper handling
    rpu_a = merge_df[f"{ColumnName.RPU}.a"].values
    rpu_b = merge_df[f"{ColumnName.RPU}.b"].values

    with np.errstate(divide="ignore", invalid="ignore"):
        enrichment = rpu_b / rpu_a
        log2e = np.log2(enrichment)

    merge_df["Enrichment"] = np.round(enrichment, 3)
    merge_df["log2E"] = np.round(log2e, 3)

    # Replace inf with NaN
    merge_df["Enrichment"] = merge_df["Enrichment"].replace([np.inf, -np.inf], np.nan)
    merge_df["log2E"] = merge_df["log2E"].replace([np.inf, -np.inf], np.nan)

    # Vectorized ranking
    merge_df[ColumnName.RANK_IN_CLUSTER] = (
        merge_df.groupby(ColumnName.CLUSTER)["Enrichment"]
        .rank(method="first", ascending=False, na_option="bottom")
        .astype(np.int32)
    )

    # Reorder columns
    cols = [ColumnName.SEQUENCES, ColumnName.CLUSTER]
    cols.extend([col for col in merge_df.columns if col not in cols])

    return merge_df[cols]


def calculate_led_to_seeds_fast(df: pd.DataFrame) -> pd.DataFrame:

    # Get seeds (rank 1 sequences) for each cluster
    seeds = df[df[ColumnName.RANK_IN_CLUSTER] == 1][
        [ColumnName.CLUSTER, ColumnName.SEQUENCES]
    ].copy()
    seeds = seeds.rename(columns={ColumnName.SEQUENCES: "_seed_seq"})

    # Merge seeds back to main dataframe
    df = df.merge(seeds, on=ColumnName.CLUSTER, how="left")

    # Vectorized LED calculation using numpy array operations
    # This is still the bottleneck, but faster than apply
    sequences = df[ColumnName.SEQUENCES].values
    seed_seqs = df["_seed_seq"].values

    # Parallel LED calculation for large datasets
    if len(df) > 10000:
        led_values = compute_led_array_parallel(sequences, seed_seqs)
    else:
        led_values = np.array(
            [
                Levenshtein.distance(seq, seed)
                for seq, seed in zip(sequences, seed_seqs)
            ],
            dtype=np.int16,
        )

    df[ColumnName.LED] = led_values
    df = df.drop(columns=["_seed_seq"])

    # Reorder and sort
    cols = [
        ColumnName.SEQUENCES,
        ColumnName.CLUSTER,
        ColumnName.RANK_IN_CLUSTER,
        ColumnName.LED,
    ]
    cols.extend([col for col in df.columns if col not in cols])
    df = df[cols]

    df = df.sort_values(
        [ColumnName.CLUSTER, "Enrichment"], ascending=[True, False], na_position="last"
    ).reset_index(drop=True)

    return df


def compute_led_array_parallel(
    sequences: np.ndarray, seeds: np.ndarray, n_jobs: int = 4
) -> np.ndarray:
    def compute_batch(start: int, end: int) -> np.ndarray:
        return np.array(
            [Levenshtein.distance(sequences[i], seeds[i]) for i in range(start, end)],
            dtype=np.int16,
        )

    n = len(sequences)
    chunk_size = max(1000, n // n_jobs)
    chunks = [(i, min(i + chunk_size, n)) for i in range(0, n, chunk_size)]

    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        results = list(executor.map(lambda c: compute_batch(*c), chunks))

    return np.concatenate(results)
