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
    enrichment_type: str = "avg"       # "avg" (AvgRPU) or "seed" (SeedRPU)
    enrichment_threshold: float = 1.0  # log2E cutoff: >= threshold → enriched, <= -threshold → depleted


class ReclusterMultiInput(BaseModel):
    fadf1_cluster_path: str = ""
    fadf2_cluster_path: str = ""
    fadf3_cluster_path: str = ""
    round1_label: str = "R1"
    round2_label: str = "R2"
    round3_label: str = "R3"
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

    file1_stem = filepath1.stem
    file2_stem = filepath2.stem

    cluster_output_path = UPLOAD_DIR / f"{file1_stem}_{file2_stem}_reclustered_led_{params.led_threshold}_clusters.{params.output_format}"
    sequence_output_path = UPLOAD_DIR / f"{file1_stem}_{file2_stem}_reclustered_led_{params.led_threshold}_sequences.{params.output_format}"

    try:
        # Read clustered data from both populations
        fadf1_cluster = read_file(filepath1)
        fadf2_cluster = read_file(filepath2)

        # Perform reclustering — returns cluster-level summary and sequence-level merged data
        summary_df, sequence_df = fa_recluster(
            fadf1_cluster=fadf1_cluster,
            fadf2_cluster=fadf2_cluster,
            led_threshold=params.led_threshold,
            enrichment_type=params.enrichment_type,
            enrichment_threshold=params.enrichment_threshold,
        )

        # Save both outputs
        save_sequences(summary_df, str(cluster_output_path), params.output_format)
        save_sequences(sequence_df, str(sequence_output_path), params.output_format)

        return {
            "status": "ok",
            "result": cluster_output_path.name,
            "result_sequences": sequence_output_path.name,
            "num_clusters": len(summary_df),
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Reclustering failed ({type(e).__name__}): {str(e)}")


@router.post("/recluster-multi")
async def fa_recluster_multi_endpoint(params: ReclusterMultiInput):
    """
    3-way cluster merge for multi-round SELEX analysis.
    """
    if not params.fadf1_cluster_path or not params.fadf2_cluster_path or not params.fadf3_cluster_path:
        raise HTTPException(status_code=400, detail="All three cluster file paths are required")

    filepath1 = UPLOAD_DIR / params.fadf1_cluster_path
    filepath2 = UPLOAD_DIR / params.fadf2_cluster_path
    filepath3 = UPLOAD_DIR / params.fadf3_cluster_path
    file1_stem = filepath1.stem
    file2_stem = filepath2.stem
    file3_stem = filepath3.stem

    labels = (params.round1_label, params.round2_label, params.round3_label)
    output_path = UPLOAD_DIR / f"{file1_stem}_{file2_stem}_{file3_stem}_recluster_multi_led_{params.led_threshold}.{params.output_format}"

    try:
        fadf1 = read_file(filepath1)
        fadf2 = read_file(filepath2)
        fadf3 = read_file(filepath3)

        result_df = fa_recluster_multi(fadf1, fadf2, fadf3, params.led_threshold, labels)
        save_sequences(result_df, str(output_path), params.output_format)

        return {
            "status": "ok",
            "result": output_path.name,
            "num_clusters": len(result_df),
            "labels": list(labels),
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Multi-round reclustering failed ({type(e).__name__}): {str(e)}")


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
            status_code=500, detail=f"LED matrix computation failed ({type(e).__name__}): {str(e)}"
        )


def fa_recluster(
    fadf1_cluster: pd.DataFrame,
    fadf2_cluster: pd.DataFrame,
    led_threshold: int = 7,
    enrichment_type: str = "avg",
    enrichment_threshold: float = 1.0,
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

    summary_df = summarize_to_cluster_level(result_df, enrichment_type, enrichment_threshold)
    return summary_df, result_df


def summarize_to_cluster_level(
    merged_df: pd.DataFrame,
    enrichment_type: str = "avg",
    enrichment_threshold: float = 1.0,
) -> pd.DataFrame:
    """Aggregate sequence-level merged data to cluster-level summary."""
    results = []
    rpu_a_col = f'{ColumnName.RPU}.a'
    rpu_b_col = f'{ColumnName.RPU}.b'
    oc_a_col = 'OriginalCluster.a'
    oc_b_col = 'OriginalCluster.b'

    for cluster_id, group in merged_df.groupby(ColumnName.CLUSTER):
        seed_row = group[group[ColumnName.RANK_IN_CLUSTER] == 1]
        seed = seed_row[ColumnName.SEQUENCES].iloc[0] if len(seed_row) > 0 else group[ColumnName.SEQUENCES].iloc[0]

        pop1_rows = group[group[oc_a_col].notna()]
        pop2_rows = group[group[oc_b_col].notna()]

        has_pop1 = len(pop1_rows) > 0
        has_pop2 = len(pop2_rows) > 0

        size_pop1 = len(pop1_rows) if has_pop1 else np.nan
        size_pop2 = len(pop2_rows) if has_pop2 else np.nan
        avg_rpu_pop1 = pop1_rows[rpu_a_col].mean() if has_pop1 else np.nan
        avg_rpu_pop2 = pop2_rows[rpu_b_col].mean() if has_pop2 else np.nan

        # Seed RPU: max RPU within each original cluster (cluster seed = highest RPU sequence)
        # RankInCluster is dropped before merging, so we use max RPU per OriginalCluster group
        # If multiple original clusters merged into this super-cluster, take the highest
        if has_pop1:
            seed_rpu_pop1 = pop1_rows.groupby(oc_a_col)[rpu_a_col].max().max()
        else:
            seed_rpu_pop1 = np.nan
        if has_pop2:
            seed_rpu_pop2 = pop2_rows.groupby(oc_b_col)[rpu_b_col].max().max()
        else:
            seed_rpu_pop2 = np.nan

        # Enrichment based on AvgRPU
        if has_pop1 and has_pop2 and not np.isnan(avg_rpu_pop1) and avg_rpu_pop1 > 0:
            raw_enrichment = float(avg_rpu_pop2 / avg_rpu_pop1)
            enrichment = round(raw_enrichment, 3)
            with np.errstate(divide='ignore', invalid='ignore'):
                log2e = round(float(np.log2(raw_enrichment)), 3)
        else:
            enrichment = np.nan
            log2e = np.nan

        # Seed enrichment based on SeedRPU
        if not np.isnan(seed_rpu_pop1) and not np.isnan(seed_rpu_pop2) and seed_rpu_pop1 > 0:
            raw_seed_enrichment = float(seed_rpu_pop2 / seed_rpu_pop1)
            seed_enrichment = round(raw_seed_enrichment, 3)
            with np.errstate(divide='ignore', invalid='ignore'):
                seed_log2e = round(float(np.log2(raw_seed_enrichment)), 3)
        else:
            seed_enrichment = np.nan
            seed_log2e = np.nan

        # Status: compare log2E against the threshold (symmetric: enriched >= t, depleted <= -t)
        metric = seed_log2e if enrichment_type == "seed" else log2e
        if has_pop1 and has_pop2:
            if not np.isnan(metric) and metric >= enrichment_threshold:
                status = 'enriched'
            elif not np.isnan(metric) and metric <= -enrichment_threshold:
                status = 'depleted'
            else:
                status = 'inherited'
        elif has_pop2:
            status = 'emerged'
        else:
            status = 'lost'

        results.append({
            'SuperCluster': int(cluster_id),
            'Seed': seed,
            'Size.Pop1': int(size_pop1) if not np.isnan(size_pop1) else np.nan,
            'Size.Pop2': int(size_pop2) if not np.isnan(size_pop2) else np.nan,
            'AvgRPU.Pop1': round(float(avg_rpu_pop1), 3) if not np.isnan(avg_rpu_pop1) else np.nan,
            'AvgRPU.Pop2': round(float(avg_rpu_pop2), 3) if not np.isnan(avg_rpu_pop2) else np.nan,
            'Enrichment': enrichment,
            'log2E': log2e,
            'SeedRPU.Pop1': round(float(seed_rpu_pop1), 3) if not np.isnan(seed_rpu_pop1) else np.nan,
            'SeedRPU.Pop2': round(float(seed_rpu_pop2), 3) if not np.isnan(seed_rpu_pop2) else np.nan,
            'SeedEnrichment': seed_enrichment,
            'SeedLog2E': seed_log2e,
            'Status': status,
        })

    return pd.DataFrame(results).sort_values('SuperCluster').reset_index(drop=True)


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

    # Coalesce into super-cluster and preserve originals with readable names
    merge_df[ColumnName.CLUSTER] = (
        merge_df[f"{ColumnName.CLUSTER}.a"]
        .fillna(merge_df[f"{ColumnName.CLUSTER}.b"])
        .astype(np.int32)
    )
    merge_df = merge_df.rename(columns={
        f"{ColumnName.CLUSTER}.a": "OriginalCluster.a",
        f"{ColumnName.CLUSTER}.b": "OriginalCluster.b",
    })

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

    led_values = np.array(
        [Levenshtein.distance(seq, seed) for seq, seed in zip(sequences, seed_seqs)],
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


def fa_recluster_multi(
    fadf1: pd.DataFrame,
    fadf2: pd.DataFrame,
    fadf3: pd.DataFrame,
    led_threshold: int = 7,
    labels: tuple = ('R1', 'R2', 'R3'),
) -> pd.DataFrame:
    """
    3-way cluster merge using union-find on pairwise LED comparisons.
    Returns cluster-level summary table with per-round stats.
    """
    def get_seeds_and_ids(df):
        mask = df[ColumnName.RANK_IN_CLUSTER].astype(int) == 1
        return df.loc[mask, ColumnName.SEQUENCES].values, df.loc[mask, ColumnName.CLUSTER].astype(int).values

    seeds1, ids1 = get_seeds_and_ids(fadf1)
    seeds2, ids2 = get_seeds_and_ids(fadf2)
    seeds3, ids3 = get_seeds_and_ids(fadf3)

    led_12 = compute_led_matrix_fast(seeds1, seeds2)
    led_23 = compute_led_matrix_fast(seeds2, seeds3)

    # Union-Find
    parent: dict = {}

    def find(x):
        parent.setdefault(x, x)
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    all_nodes = (
        [f'r1_{i}' for i in range(len(seeds1))] +
        [f'r2_{i}' for i in range(len(seeds2))] +
        [f'r3_{i}' for i in range(len(seeds3))]
    )
    for n in all_nodes:
        find(n)

    for i in range(len(seeds1)):
        cols = np.where(led_12[i] < led_threshold)[0]
        if len(cols) > 0:
            union(f'r1_{i}', f'r2_{cols[np.argmin(led_12[i, cols])]}')

    for j in range(len(seeds2)):
        cols = np.where(led_23[j] < led_threshold)[0]
        if len(cols) > 0:
            union(f'r2_{j}', f'r3_{cols[np.argmin(led_23[j, cols])]}')

    root_to_sc: dict = {}
    sc_counter = 1
    node_to_sc: dict = {}
    for node in all_nodes:
        root = find(node)
        if root not in root_to_sc:
            root_to_sc[root] = sc_counter
            sc_counter += 1
        node_to_sc[node] = root_to_sc[root]

    def build_orig_to_sc(ids, prefix):
        return {int(orig_id): node_to_sc[f'{prefix}_{i}'] for i, orig_id in enumerate(ids)}

    orig_to_sc1 = build_orig_to_sc(ids1, 'r1')
    orig_to_sc2 = build_orig_to_sc(ids2, 'r2')
    orig_to_sc3 = build_orig_to_sc(ids3, 'r3')

    def round_stats(df, orig_to_sc, label):
        df = df.copy()
        df['SuperCluster'] = df[ColumnName.CLUSTER].astype(int).map(orig_to_sc)
        grouped = df.groupby('SuperCluster').agg(
            **{
                f'Size.{label}': (ColumnName.SEQUENCES, 'count'),
                f'AvgRPU.{label}': (ColumnName.RPU, 'mean'),
            }
        ).reset_index()
        grouped[f'AvgRPU.{label}'] = grouped[f'AvgRPU.{label}'].round(3)
        # SeedRPU: max RPU among original cluster seeds (rank 1) per super-cluster
        seed_rows = df[df[ColumnName.RANK_IN_CLUSTER].astype(int) == 1]
        seed_rpu = seed_rows.groupby('SuperCluster')[ColumnName.RPU].max().reset_index()
        seed_rpu = seed_rpu.rename(columns={ColumnName.RPU: f'SeedRPU.{label}'})
        seed_rpu[f'SeedRPU.{label}'] = seed_rpu[f'SeedRPU.{label}'].round(3)
        grouped = grouped.merge(seed_rpu, on='SuperCluster', how='left')
        return grouped

    r1, r2, r3 = labels
    stats1 = round_stats(fadf1, orig_to_sc1, r1)
    stats2 = round_stats(fadf2, orig_to_sc2, r2)
    stats3 = round_stats(fadf3, orig_to_sc3, r3)

    result = stats1.merge(stats2, on='SuperCluster', how='outer')
    result = result.merge(stats3, on='SuperCluster', how='outer')

    # Seed from latest round available
    sc_to_seed: dict = {}
    for i, orig_id in enumerate(ids3):
        sc_to_seed[node_to_sc[f'r3_{i}']] = seeds3[i]
    for i, orig_id in enumerate(ids2):
        sc = node_to_sc[f'r2_{i}']
        if sc not in sc_to_seed:
            sc_to_seed[sc] = seeds2[i]
    for i, orig_id in enumerate(ids1):
        sc = node_to_sc[f'r1_{i}']
        if sc not in sc_to_seed:
            sc_to_seed[sc] = seeds1[i]

    result['Seed'] = result['SuperCluster'].map(sc_to_seed)

    # Enrichment between consecutive rounds (AvgRPU-based)
    e_col_12 = f'E.{r1}.{r2}'
    e_col_23 = f'E.{r2}.{r3}'
    with np.errstate(divide='ignore', invalid='ignore'):
        result[e_col_12] = (result[f'AvgRPU.{r2}'] / result[f'AvgRPU.{r1}']).round(3)
        result[e_col_23] = (result[f'AvgRPU.{r3}'] / result[f'AvgRPU.{r2}']).round(3)
    result[e_col_12] = result[e_col_12].replace([np.inf, -np.inf], np.nan)
    result[e_col_23] = result[e_col_23].replace([np.inf, -np.inf], np.nan)

    # Seed enrichment between consecutive rounds (SeedRPU-based)
    se_col_12 = f'SeedE.{r1}.{r2}'
    se_col_23 = f'SeedE.{r2}.{r3}'
    with np.errstate(divide='ignore', invalid='ignore'):
        result[se_col_12] = (result[f'SeedRPU.{r2}'] / result[f'SeedRPU.{r1}']).round(3)
        result[se_col_23] = (result[f'SeedRPU.{r3}'] / result[f'SeedRPU.{r2}']).round(3)
    result[se_col_12] = result[se_col_12].replace([np.inf, -np.inf], np.nan)
    result[se_col_23] = result[se_col_23].replace([np.inf, -np.inf], np.nan)

    cols = ['SuperCluster', 'Seed',
            f'Size.{r1}', f'Size.{r2}', f'Size.{r3}',
            f'AvgRPU.{r1}', f'AvgRPU.{r2}', f'AvgRPU.{r3}',
            f'SeedRPU.{r1}', f'SeedRPU.{r2}', f'SeedRPU.{r3}',
            e_col_12, e_col_23, se_col_12, se_col_23]
    result = result[cols].sort_values('SuperCluster').reset_index(drop=True)

    return result


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
