from pathlib import Path
from Bio import SeqIO
import numpy as np
import pandas as pd
import re
from rapidfuzz import process
from rapidfuzz.distance import Levenshtein
import time
from services.file_service import parse_fasta,save_sequences
from services.constants import ColumnName

def run_cluster_led(input_path, min_reads=10, max_led=7, total_clusters=30, keep_nc=True, output_path=None, output_format='fasta'):
    
    df = parse_fasta(input_path)
    n_input = len(df)

    # Filter by read count, initialize cluster-specific columns (sort_values already copies; drop df to free memory)
    cluster_data = df.sort_values(ColumnName.RPU, ascending=False)
    del df
    cluster_data = cluster_data[cluster_data[ColumnName.READS] >= min_reads].reset_index(drop=True)
    if len(cluster_data) == 0:
        print(f"WARNING: No sequences remain after filtering with min_reads={min_reads}")
        print(f"Original data had {n_input} sequences")
        return None
    
    cluster_data[ColumnName.CLUSTER] = pd.NA
    cluster_data[ColumnName.RANK_IN_CLUSTER] = pd.NA
    cluster_data[ColumnName.LED] = pd.NA
    
    # Store original ID for reference
    cluster_data[ColumnName.ORIGINAL_ID] = cluster_data[ColumnName.ID]
    cluster_data = cluster_data.drop(columns=[ColumnName.ID])
    
    # Unclustered sequences, most abundant first
    sequences = np.array(cluster_data[ColumnName.SEQUENCES].tolist(), dtype=object)

    # Create cluster iterator
    cluster_number = 1

    # Sequence -> cluster number, rank in cluster, and LED to seed
    cluster_of, rank_of, led_of = {}, {}, {}

    # Only create a user-specified number of clusters and while you still have sequences
    while cluster_number <= total_clusters and len(sequences) > 0:

        # Get start time
        start_time = time.time()

        # First element of seq. list (most abundant seq.) is seed; remove from seq. list
        seed = sequences[0]
        sequences = sequences[1:]

        # LED between the seed and every remaining sequence in one multithreaded call (cdist splits
        # work across threads by query row, so the remaining sequences are the queries);
        # distances above max_led are returned as max_led + 1
        if len(sequences) > 0:
            led = process.cdist(sequences, [seed], scorer=Levenshtein.distance, score_cutoff=max_led, workers=-1)[:, 0]
        else:
            led = np.empty(0, dtype=np.int32)
        in_cluster = led <= max_led

        # Seed first, then members in abundance order, with their LED to the seed
        cluster_list = [seed, *sequences[in_cluster]]
        seq_led = [0, *led[in_cluster].tolist()]

        # Remove clustered sequences from the unclustered list
        sequences = sequences[~in_cluster]

        # Record cluster information per sequence; written back to cluster_data once after the loop
        for rank, (seq, led) in enumerate(zip(cluster_list, seq_led), start=1):
            cluster_of[seq] = cluster_number
            rank_of[seq] = rank
            led_of[seq] = led

        # Message with number of unique sequences in cluster
        # elapsed_time = time.time() - start_time
        # print(f"Finished cluster {cluster_number}: {len(cluster_list)} unique sequences")
        # print(f"Elapsed time: {elapsed_time:.2f} seconds\n")
        
        # Bump the cluster number
        cluster_number += 1

    # Update cluster_data with cluster information (unclustered sequences stay NA)
    cluster_data[ColumnName.CLUSTER] = cluster_data[ColumnName.SEQUENCES].map(cluster_of)
    cluster_data[ColumnName.RANK_IN_CLUSTER] = cluster_data[ColumnName.SEQUENCES].map(rank_of)
    cluster_data[ColumnName.LED] = cluster_data[ColumnName.SEQUENCES].map(led_of)

    # Keep or omit non-clustered sequences
    if keep_nc:
        # Replace NAs with 0
        cluster_data[[ColumnName.CLUSTER, ColumnName.RANK_IN_CLUSTER, ColumnName.LED]] = cluster_data[[ColumnName.CLUSTER, ColumnName.RANK_IN_CLUSTER, ColumnName.LED]].fillna(0)
        
        # Separate cluster 0
        cluster0 = cluster_data[cluster_data[ColumnName.CLUSTER] == 0].copy()
        cluster_data = cluster_data[cluster_data[ColumnName.CLUSTER] != 0]
        
        # Add rank in cluster (determined by RPU)
        if len(cluster0) > 0:
            cluster0 = cluster0.sort_values(ColumnName.RPU, ascending=False).reset_index(drop=True)
            cluster0[ColumnName.RANK_IN_CLUSTER] = range(1, len(cluster0) + 1)
            
            # Compute LED between the cluster seed and every other sequence
            seed_sequence = cluster0.iloc[0][ColumnName.SEQUENCES]
            cluster0[ColumnName.LED] = process.cdist(cluster0[ColumnName.SEQUENCES].tolist(), [seed_sequence], scorer=Levenshtein.distance, workers=-1)[:, 0]

        # Add cluster 0 back to main data
        cluster_data = pd.concat([cluster_data, cluster0], ignore_index=True)
    else:
        cluster_data = cluster_data.dropna(subset=[ColumnName.CLUSTER, ColumnName.RANK_IN_CLUSTER, ColumnName.LED])
    
    # Ensure numeric types for the ID columns
    cluster_data[ColumnName.CLUSTER] = cluster_data[ColumnName.CLUSTER].astype(int)
    cluster_data[ColumnName.RANK_IN_CLUSTER] = cluster_data[ColumnName.RANK_IN_CLUSTER].astype(int)
    cluster_data[ColumnName.LED] = cluster_data[ColumnName.LED].astype(int)
    
    # Add cluster information to sequence IDs (column-wise; a row-wise apply takes minutes on millions of rows)
    cluster_data[ColumnName.ID] = (
        f"{ColumnName.RANK}=" + cluster_data[ColumnName.RANK].astype(str)
        + f";{ColumnName.READS}=" + cluster_data[ColumnName.READS].astype(str)
        + f";{ColumnName.RPU}=" + cluster_data[ColumnName.RPU].astype(str)
        + f";{ColumnName.CLUSTER}=" + cluster_data[ColumnName.CLUSTER].astype(str)
        + f";{ColumnName.RANK_IN_CLUSTER}=" + cluster_data[ColumnName.RANK_IN_CLUSTER].astype(str)
        + f";{ColumnName.LED}=" + cluster_data[ColumnName.LED].astype(str)
    )
    
    # Reorder columns to put ID first
    cols = [ColumnName.ID, ColumnName.RANK, ColumnName.READS,ColumnName.RPU,ColumnName.CLUSTER,ColumnName.RANK_IN_CLUSTER,ColumnName.SEQUENCES]
    cluster_data = cluster_data[cols]

    save_sequences(cluster_data, output_path, output_format)
    
    return output_path
    
def run_cluster_diversity(
    inputpath,
    output_format='csv',
    output_path=None
):

    # Read and parse clustered FASTA file
    df = parse_fasta(inputpath)
    
    # Verify required columns
    required_cols = [ColumnName.CLUSTER, ColumnName.READS, ColumnName.RPU, ColumnName.LED, ColumnName.RANK_IN_CLUSTER]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Input file missing required columns: {missing_cols}. "
                        f"Please use a clustered FASTA file.")
    
    print(f"Read {len(df)} sequences across {df['Cluster'].nunique()} clusters")
    
    print("\nCalculating diversity statistics...")
    
    # Calculate diversity statistics
    diversity_df = cluster_diversity(df)
    
    print(f"\nDiversity analysis complete!")
    print(f"Analyzed {len(diversity_df)} clusters")
    
    # Generate output path if not provided
    if output_path is None:
        input_path = Path(inputpath)
        output_path = str(input_path.parent / f"{input_path.stem}_diversity.{output_format}")
    
    # Write output
    output_format = output_format.lower()
    if output_format == 'csv':
        diversity_df.to_csv(output_path, index=False)
    else:
        raise ValueError(f"Unsupported output format: {output_format}")
    
    print(f"Written diversity statistics to {output_path}")
    
    return output_path

def cluster_diversity(fa_df_cluster):
  
    # Calculate summary features per cluster
    cluster_stats = fa_df_cluster.groupby(ColumnName.CLUSTER).agg(
        TotalSequences=(ColumnName.CLUSTER, 'size'),
        TotalReads=(ColumnName.READS, 'sum'),
        TotalRPU=(ColumnName.RPU, 'sum'),
        AverageLED=(ColumnName.LED, 'mean')
    ).reset_index()
    
    # Round AverageLED to 2 decimal places
    cluster_stats[ColumnName.AVERAGE_LED] = cluster_stats[ColumnName.AVERAGE_LED].round(2)
    
    # Calculate SID for each cluster
    sid_values = fa_df_cluster.groupby(ColumnName.CLUSTER).apply(calculate_sid).reset_index()
    sid_values.columns = [ColumnName.CLUSTER, ColumnName.SID]
    
    # Merge SID values with cluster stats
    cluster_stats = cluster_stats.merge(sid_values, on=ColumnName.CLUSTER)
    
    # Get cluster seeds (RankInCluster == 1) with their IDs and sequences
    cluster_seeds = fa_df_cluster[fa_df_cluster[ColumnName.RANK_IN_CLUSTER] == 1][[ColumnName.CLUSTER, ColumnName.ID, ColumnName.SEQUENCES]]
    
    # Left join to add cluster seed information
    cluster_stats = cluster_stats.merge(cluster_seeds, on=ColumnName.CLUSTER, how='left')
    
    # Sort by cluster
    cluster_stats = cluster_stats.sort_values(ColumnName.CLUSTER).reset_index(drop=True)
    
    # Reorder columns for better readability
    column_order = [ColumnName.CLUSTER, ColumnName.TOTAL_SEQUENCES, ColumnName.TOTAL_READS, ColumnName.TOTAL_RPU, 
                    ColumnName.AVERAGE_LED, ColumnName.SID, ColumnName.ID, ColumnName.SEQUENCES]
    cluster_stats = cluster_stats[column_order]
    
    return cluster_stats

def calculate_sid(group):

    reads = group[ColumnName.READS].values
    total_reads = reads.sum()
        
    if total_reads <= 1:
        return 0.0
        
    numerator = sum(r * (r - 1) for r in reads)
    denominator = total_reads * (total_reads - 1)
        
    sid = 1 - (numerator / denominator)
    return round(sid, 2)
