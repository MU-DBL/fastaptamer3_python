from typing import Dict, List, Any
from services.file_service import read_file, save_sequences, check_same_scaling_factor, scaling_factor_of, RPU_DECIMALS
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from pathlib import Path
import os
import pandas as pd
import numpy as np
import Levenshtein
from rapidfuzz import process
from rapidfuzz.distance import Levenshtein as RapidLevenshtein
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
    led_threshold: int = 0             # seed LED for merging; 0 = identical seeds only (as in FASTAptameR 2.0)
    output_format: str = "csv"
    enrichment_type: str = "total"     # "total" (TotalRPU), "avg" (AvgRPU) or "seed" (SeedRPU)
    enrichment_threshold: float = 1.0  # log2E cutoff: >= threshold → enriched, <= -threshold → depleted


class ReclusterMultiInput(BaseModel):
    fadf1_cluster_path: str = ""
    fadf2_cluster_path: str = ""
    fadf3_cluster_path: str = ""
    round1_label: str = "R1"
    round2_label: str = "R2"
    round3_label: str = "R3"
    led_threshold: int = 3             # seed LED for linking consecutive rounds
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
        check_same_scaling_factor([fadf1_cluster, fadf2_cluster], [params.fadf1_cluster_path, params.fadf2_cluster_path])

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
        check_same_scaling_factor(
            [fadf1, fadf2, fadf3],
            [params.fadf1_cluster_path, params.fadf2_cluster_path, params.fadf3_cluster_path],
        )

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

        # Extract seeds and cluster IDs (Cluster 0 holds non-clustered sequences and has no seed)
        p1_cluster_ids, p1_seeds = get_cluster_seeds(fadf1_cluster)
        p2_cluster_ids, p2_seeds = get_cluster_seeds(fadf2_cluster)
        p1_cluster_ids, p1_seeds = p1_cluster_ids.tolist(), p1_seeds.tolist()
        p2_cluster_ids, p2_seeds = p2_cluster_ids.tolist(), p2_seeds.tolist()

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
    led_threshold: int = 0,
    enrichment_type: str = "total",
    enrichment_threshold: float = 1.0,
    use_parallel: bool = True,
    n_jobs: Optional[int] = None,
) -> pd.DataFrame:

    # Seeds of the real clusters (Cluster 0 holds non-clustered sequences and has no seed)
    p1_cluster_ids, p1_seeds = get_cluster_seeds(fadf1_cluster)
    p2_cluster_ids, p2_seeds = get_cluster_seeds(fadf2_cluster)

    led_matrix = compute_led_matrix_fast(p1_seeds, p2_seeds)
    assignment = assign_pop2_clusters(led_matrix, p1_cluster_ids, p2_cluster_ids, led_threshold)

    merged_df, pop2_mapping = merge_clusters_batch(fadf1_cluster, fadf2_cluster, assignment)

    # Reference seed of each super-cluster: the Pop1 seed, or the Pop2 seed for clusters found only in Pop2
    ref_seeds = {int(cid): seed for cid, seed in zip(p1_cluster_ids, p1_seeds)}
    for cid, seed in zip(p2_cluster_ids, p2_seeds):
        if int(cid) not in assignment:
            ref_seeds[pop2_mapping[int(cid)]] = seed
    # SeedLED: the largest seed-to-seed LED among the Pop2 clusters merged into each Pop1 cluster
    seed_led = {}
    for pop1_id, led in assignment.values():
        seed_led[pop1_id] = max(seed_led.get(pop1_id, 0), led)

    merge_df = calculate_enrichment_vectorized(merged_df)

    result_df = calculate_led_to_seeds_fast(merge_df)

    summary_df = summarize_to_cluster_level(
        result_df, enrichment_type, enrichment_threshold, ref_seeds=ref_seeds, seed_led=seed_led
    )
    return summary_df, result_df


def summarize_to_cluster_level(
    merged_df: pd.DataFrame,
    enrichment_type: str = "total",
    enrichment_threshold: float = 1.0,
    ref_seeds: Optional[Dict[int, str]] = None,
    seed_led: Optional[Dict[int, int]] = None,
) -> pd.DataFrame:
    """
    Aggregate sequence-level merged data to cluster-level summary.

    ref_seeds maps each super-cluster to its reference seed (original Pop1 seed, or Pop2 seed for
    Pop2-only clusters); member distances are measured from it. seed_led maps merged super-clusters
    to the seed-to-seed LED of their merge.
    """
    ref_seeds = ref_seeds or {}
    seed_led = seed_led or {}
    results = []
    rpu_a_col = f'{ColumnName.RPU}.a'
    rpu_b_col = f'{ColumnName.RPU}.b'
    oc_a_col = 'OriginalCluster.a'
    oc_b_col = 'OriginalCluster.b'

    for cluster_id, group in merged_df.groupby(ColumnName.CLUSTER):
        ref_seed = ref_seeds.get(int(cluster_id))
        if ref_seed is not None:
            seed = ref_seed
            members = group[ColumnName.SEQUENCES].tolist()
            led_to_seed = process.cdist([ref_seed], members, scorer=RapidLevenshtein.distance, dtype=np.int32, workers=-1)[0]
            median_led = float(np.median(led_to_seed))
            max_led = int(led_to_seed.max())
            diameter = cluster_diameter(members)
        else:
            # Non-clustered sequences (Cluster 0) have no seed, so no distance statistics
            seed_row = group[group[ColumnName.RANK_IN_CLUSTER] == 1]
            seed = seed_row[ColumnName.SEQUENCES].iloc[0] if len(seed_row) > 0 else group[ColumnName.SEQUENCES].iloc[0]
            median_led = max_led = diameter = np.nan

        pop1_rows = group[group[oc_a_col].notna()]
        pop2_rows = group[group[oc_b_col].notna()]

        has_pop1 = len(pop1_rows) > 0
        has_pop2 = len(pop2_rows) > 0

        size_pop1 = len(pop1_rows) if has_pop1 else np.nan
        size_pop2 = len(pop2_rows) if has_pop2 else np.nan
        total_rpu_pop1 = pop1_rows[rpu_a_col].sum() if has_pop1 else np.nan
        total_rpu_pop2 = pop2_rows[rpu_b_col].sum() if has_pop2 else np.nan
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

        # Total enrichment based on TotalRPU (the whole family's share of each population)
        if has_pop1 and has_pop2 and total_rpu_pop1 > 0:
            raw_total_enrichment = float(total_rpu_pop2 / total_rpu_pop1)
            total_enrichment = round(raw_total_enrichment, 3)
            with np.errstate(divide='ignore', invalid='ignore'):
                total_log2e = round(float(np.log2(raw_total_enrichment)), 3)
        else:
            total_enrichment = np.nan
            total_log2e = np.nan

        # Status: compare log2E against the threshold (symmetric: enriched >= t, depleted <= -t)
        metric = {"seed": seed_log2e, "total": total_log2e}.get(enrichment_type, log2e)
        if int(cluster_id) == 0:
            status = 'unclustered'
        elif has_pop1 and has_pop2:
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
            'SeedLED': seed_led.get(int(cluster_id), np.nan),
            'MedianLED': median_led,
            'MaxLED': max_led,
            'Diameter': diameter,
            'TotalRPU.Pop1': round(float(total_rpu_pop1), RPU_DECIMALS) if not np.isnan(total_rpu_pop1) else np.nan,
            'TotalRPU.Pop2': round(float(total_rpu_pop2), RPU_DECIMALS) if not np.isnan(total_rpu_pop2) else np.nan,
            'TotalEnrichment': total_enrichment,
            'TotalLog2E': total_log2e,
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

    summary = pd.DataFrame(results).sort_values('SuperCluster').reset_index(drop=True)

    # Show the U behind the RPU columns (both populations share it; checked on input)
    u = scaling_factor_of(*(merged_df[c] for c in (f'{ColumnName.U}.a', f'{ColumnName.U}.b') if c in merged_df.columns))
    if u is not None:
        summary.insert(summary.columns.get_loc('TotalRPU.Pop1'), ColumnName.U, u)
    return summary


def compute_led_matrix_fast(seeds1: np.ndarray, seeds2: np.ndarray) -> np.ndarray:
    # All pairwise LEDs in one multithreaded call
    return process.cdist(list(seeds1), list(seeds2), scorer=RapidLevenshtein.distance, dtype=np.int16, workers=-1)


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


def get_cluster_seeds(df: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
    """Cluster IDs and seed sequences (RankInCluster 1) of a Cluster output, excluding non-clustered sequences (Cluster 0)."""
    mask = (df[ColumnName.RANK_IN_CLUSTER].astype(int) == 1) & (df[ColumnName.CLUSTER].astype(int) != 0)
    return df.loc[mask, ColumnName.CLUSTER].astype(int).to_numpy(), df.loc[mask, ColumnName.SEQUENCES].to_numpy()


def assign_pop2_clusters(
    led_matrix: np.ndarray, p1_cluster_ids: np.ndarray, p2_cluster_ids: np.ndarray, threshold: int
) -> Dict[int, Tuple[int, int]]:
    """
    Assign each Pop2 cluster to the Pop1 cluster whose seed is nearest to its own seed, when
    that seed LED is <= threshold. A Pop1 cluster can take several Pop2 clusters (a family split
    in Pop2); each Pop2 cluster joins at most one. Ties go to the Pop1 cluster listed first.

    Returns {Pop2 cluster ID: (Pop1 cluster ID, seed LED)} for the merged Pop2 clusters.
    """
    if led_matrix.size == 0:
        return {}
    nearest = np.argmin(led_matrix, axis=0)  # first minimum, i.e. the earliest Pop1 cluster on ties
    nearest_led = led_matrix[nearest, np.arange(led_matrix.shape[1])]
    return {
        int(p2_cluster_ids[j]): (int(p1_cluster_ids[nearest[j]]), int(nearest_led[j]))
        for j in np.flatnonzero(nearest_led <= threshold)
    }


def merge_clusters_batch(
    fadf1: pd.DataFrame,
    fadf2: pd.DataFrame,
    assignment: Dict[int, Tuple[int, int]],
) -> Tuple[pd.DataFrame, Dict[int, int]]:
    """
    Combine both populations under super-cluster IDs. Pop1 clusters keep their IDs, merged Pop2
    clusters take their Pop1 cluster's ID, unmerged Pop2 clusters get new IDs after Pop1's, and
    non-clustered sequences (Cluster 0) of both populations stay in Cluster 0.

    Returns the combined data and the Pop2 cluster ID -> super-cluster ID mapping.
    """
    pop1 = fadf1.drop(columns=[ColumnName.RANK_IN_CLUSTER, ColumnName.LED]).copy()
    pop2 = fadf2.drop(columns=[ColumnName.RANK_IN_CLUSTER, ColumnName.LED]).copy()

    pop1["Population"] = 1
    pop2["Population"] = 2

    pop2_ids = pop2[ColumnName.CLUSTER].astype(int)
    pop2_cluster_mapping = {pop2_id: pop1_id for pop2_id, (pop1_id, _) in assignment.items()}
    pop2_cluster_mapping[0] = 0

    next_id = max(int(pop1[ColumnName.CLUSTER].max()), 0) + 1 if len(pop1) else 1
    unmerged_pop2_clusters = sorted(set(pop2_ids.unique()) - set(pop2_cluster_mapping))
    for i, old_cluster in enumerate(unmerged_pop2_clusters):
        pop2_cluster_mapping[old_cluster] = next_id + i

    pop2[ColumnName.CLUSTER] = pop2_ids.map(pop2_cluster_mapping)

    return pd.concat([pop1, pop2], ignore_index=True), pop2_cluster_mapping


def cluster_diameter(sequences: List[str], chunk: int = 2000) -> int:
    """Largest pairwise LED among the sequences, computed in row blocks so memory stays bounded."""
    diameter = 0
    for start in range(0, len(sequences), chunk):
        # Each block only needs the sequences from its own start onwards (the matrix is symmetric)
        block = process.cdist(sequences[start:start + chunk], sequences[start:],
                              scorer=RapidLevenshtein.distance, dtype=np.int32, workers=-1)
        diameter = max(diameter, int(block.max()))
    return diameter


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
    led_threshold: int = 3,
    labels: tuple = ('R1', 'R2', 'R3'),
) -> pd.DataFrame:
    """
    3-round cluster merge, applying the two-population rule to consecutive rounds: each R2
    cluster joins the R1 cluster whose seed is nearest to its seed, and each R3 cluster the
    nearest R2 cluster, when that seed LED is <= led_threshold. Following these links back to
    the earliest round gives each super-cluster, so clusters of the same round are never merged
    with each other. Non-clustered sequences (Cluster 0) stay in super-cluster 0.
    Returns cluster-level summary table with per-round stats.
    """
    ids1, seeds1 = get_cluster_seeds(fadf1)
    ids2, seeds2 = get_cluster_seeds(fadf2)
    ids3, seeds3 = get_cluster_seeds(fadf3)

    link_12 = assign_pop2_clusters(compute_led_matrix_fast(seeds1, seeds2), ids1, ids2, led_threshold)
    link_23 = assign_pop2_clusters(compute_led_matrix_fast(seeds2, seeds3), ids2, ids3, led_threshold)

    # Super-cluster IDs: R1 clusters first, then clusters first seen in R2, then those first seen in R3.
    # ref_seeds holds each super-cluster's seed from the round it first appears in; seed_led the
    # largest seed LED among the links inside it.
    orig_to_sc1, orig_to_sc2, orig_to_sc3 = {0: 0}, {0: 0}, {0: 0}
    ref_seeds: Dict[int, str] = {}
    seed_led: Dict[int, int] = {}
    next_sc = 1
    for cid, seed in zip(ids1, seeds1):
        orig_to_sc1[int(cid)] = next_sc
        ref_seeds[next_sc] = seed
        next_sc += 1
    for orig_to_sc, prev_to_sc, ids, seeds, links in (
        (orig_to_sc2, orig_to_sc1, ids2, seeds2, link_12),
        (orig_to_sc3, orig_to_sc2, ids3, seeds3, link_23),
    ):
        for cid, seed in zip(ids, seeds):
            cid = int(cid)
            if cid in links:
                prev_cid, led = links[cid]
                sc = prev_to_sc[prev_cid]
                orig_to_sc[cid] = sc
                seed_led[sc] = max(seed_led.get(sc, 0), led)
            else:
                orig_to_sc[cid] = next_sc
                ref_seeds[next_sc] = seed
                next_sc += 1

    def round_stats(df, orig_to_sc, label):
        df = df.copy()
        df['SuperCluster'] = df[ColumnName.CLUSTER].astype(int).map(orig_to_sc)
        grouped = df.groupby('SuperCluster').agg(
            **{
                f'Size.{label}': (ColumnName.SEQUENCES, 'count'),
                f'TotalRPU.{label}': (ColumnName.RPU, 'sum'),
                f'AvgRPU.{label}': (ColumnName.RPU, 'mean'),
            }
        ).reset_index()
        grouped[f'TotalRPU.{label}'] = grouped[f'TotalRPU.{label}'].round(RPU_DECIMALS)
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

    # Seed from the round each super-cluster first appears in
    result['Seed'] = result['SuperCluster'].map(ref_seeds)

    # Merge quality: member distances to that seed and the members' diameter (none for Cluster 0)
    members = pd.concat([
        pd.DataFrame({'SuperCluster': df[ColumnName.CLUSTER].astype(int).map(orig_to_sc),
                      ColumnName.SEQUENCES: df[ColumnName.SEQUENCES]})
        for df, orig_to_sc in ((fadf1, orig_to_sc1), (fadf2, orig_to_sc2), (fadf3, orig_to_sc3))
    ]).drop_duplicates()
    distance_rows = []
    for sc, group in members.groupby('SuperCluster'):
        ref_seed = ref_seeds.get(int(sc))
        if ref_seed is None:
            continue
        seqs = group[ColumnName.SEQUENCES].tolist()
        led_to_seed = process.cdist([ref_seed], seqs, scorer=RapidLevenshtein.distance, dtype=np.int32, workers=-1)[0]
        distance_rows.append({'SuperCluster': sc, 'SeedLED': seed_led.get(int(sc), np.nan),
                              'MedianLED': float(np.median(led_to_seed)), 'MaxLED': int(led_to_seed.max()),
                              'Diameter': cluster_diameter(seqs)})
    distance_cols = ['SeedLED', 'MedianLED', 'MaxLED', 'Diameter']
    result = result.merge(pd.DataFrame(distance_rows, columns=['SuperCluster'] + distance_cols),
                          on='SuperCluster', how='left')

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

    # Total enrichment between consecutive rounds (TotalRPU-based)
    te_col_12 = f'TotalE.{r1}.{r2}'
    te_col_23 = f'TotalE.{r2}.{r3}'
    with np.errstate(divide='ignore', invalid='ignore'):
        result[te_col_12] = (result[f'TotalRPU.{r2}'] / result[f'TotalRPU.{r1}']).round(3)
        result[te_col_23] = (result[f'TotalRPU.{r3}'] / result[f'TotalRPU.{r2}']).round(3)
    result[te_col_12] = result[te_col_12].replace([np.inf, -np.inf], np.nan)
    result[te_col_23] = result[te_col_23].replace([np.inf, -np.inf], np.nan)

    cols = ['SuperCluster', 'Seed',
            f'Size.{r1}', f'Size.{r2}', f'Size.{r3}',
            *distance_cols,
            f'TotalRPU.{r1}', f'TotalRPU.{r2}', f'TotalRPU.{r3}',
            f'AvgRPU.{r1}', f'AvgRPU.{r2}', f'AvgRPU.{r3}',
            f'SeedRPU.{r1}', f'SeedRPU.{r2}', f'SeedRPU.{r3}',
            te_col_12, te_col_23, e_col_12, e_col_23, se_col_12, se_col_23]
    result = result[cols].sort_values('SuperCluster').reset_index(drop=True)

    # Show the U behind the RPU columns (all rounds share it; checked on input)
    u = scaling_factor_of(*(df[ColumnName.U] for df in (fadf1, fadf2, fadf3) if ColumnName.U in df.columns))
    if u is not None:
        result.insert(result.columns.get_loc(f'TotalRPU.{r1}'), ColumnName.U, u)

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
