from pathlib import Path
from Bio import SeqIO
import pandas as pd
import re
import Levenshtein
import time
from backend.services import file_service
from file_service import parse_fasta
from constants import ColumnName

def run_cluster_led(file_path, min_reads=10, max_led=7, total_clusters=30, keep_nc=True, output_path=None, output_format='fasta'):
    
    df = parse_fasta(file_path)

    # Filter by read count, initialize cluster-specific columns
    cluster_data = df.copy()
    cluster_data = cluster_data.sort_values(ColumnName.RPU, ascending=False)
    cluster_data = cluster_data[cluster_data[ColumnName.READS] >= min_reads].reset_index(drop=True)
    cluster_data[ColumnName.CLUSTER] = pd.NA
    cluster_data[ColumnName.RANK_IN_CLUSTER] = pd.NA
    cluster_data[ColumnName.LED] = pd.NA
    
    # Store original ID for reference
    cluster_data[ColumnName.ORIGINAL_ID] = cluster_data[ColumnName.ID]
    cluster_data = cluster_data.drop(columns=[ColumnName.ID])
    
    # Get sequences as a list
    sequences = cluster_data[ColumnName.SEQUENCES].tolist()
    
    # Create cluster iterator
    cluster_number = 1
    
    # Only create a user-specified number of clusters and while you still have sequences
    while cluster_number <= total_clusters and len(sequences) > 0:
        
        # Get start time
        start_time = time.time()
        
        # Reset iterator (counter)
        counter = 0
        
        # First element of seq. list (most abundant seq.) is seed; remove from seq. list
        cluster_list = [sequences[0]]
        sequences = sequences[1:]
        
        # Make list for string distances
        seq_led = [0]
        
        # Only iterate while you still have sequences
        while len(sequences) > 0 and counter < len(sequences):
            
            # Check LED between seed and current sequence
            if sequences[counter] is not None:
                led_distance = Levenshtein.distance(cluster_list[0], sequences[counter])
                
                if led_distance <= max_led:
                    # Add LED to seq_led list
                    seq_led.append(led_distance)
                    
                    # If LED is <= maxLED, add sequence to cluster list; mark sequence as None
                    cluster_list.append(sequences[counter])
                    sequences[counter] = None
            
            # Check next sequence in list
            counter += 1
        
        # Remove None values from sequence list; these seqs are now clustered
        sequences = [s for s in sequences if s is not None]
        
        # Turn clustered sequences into DataFrame and merge back into cluster_data
        seq_df = pd.DataFrame({
            ColumnName.SEQUENCES: cluster_list,
            ColumnName.CLUSTER: cluster_number,
            ColumnName.RANK_IN_CLUSTER: range(1, len(cluster_list) + 1),
            ColumnName.LED: seq_led
        })
        
        # Update cluster_data with cluster information
        for idx, row in seq_df.iterrows():
            mask = cluster_data[ColumnName.SEQUENCES] == row[ColumnName.SEQUENCES]
            cluster_data.loc[mask, ColumnName.CLUSTER] = row[ColumnName.CLUSTER]
            cluster_data.loc[mask, ColumnName.RANK_IN_CLUSTER] = row[ColumnName.RANK_IN_CLUSTER]
            cluster_data.loc[mask, ColumnName.LED] = row[ColumnName.LED]
        
        # Message with number of unique sequences in cluster
        elapsed_time = time.time() - start_time
        print(f"Finished cluster {cluster_number}: {len(seq_df)} unique sequences")
        print(f"Elapsed time: {elapsed_time:.2f} seconds\n")
        
        # Bump the cluster number
        cluster_number += 1
    
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
            cluster0[ColumnName.SEQUENCES] = cluster0[ColumnName.SEQUENCES].apply(lambda x: Levenshtein.distance(seed_sequence, x))
        
        # Add cluster 0 back to main data
        cluster_data = pd.concat([cluster_data, cluster0], ignore_index=True)
    else:
        cluster_data = cluster_data.dropna(subset=[ColumnName.CLUSTER, ColumnName.RANK_IN_CLUSTER, ColumnName.LED])
    
    # Ensure numeric types for the ID columns
    cluster_data[ColumnName.CLUSTER] = cluster_data[ColumnName.CLUSTER].astype(int)
    cluster_data[ColumnName.RANK_IN_CLUSTER] = cluster_data[ColumnName.RANK_IN_CLUSTER].astype(int)
    cluster_data[ColumnName.LED] = cluster_data[ColumnName.LED].astype(int)
    
    # Add cluster information to sequence IDs
    cluster_data[ColumnName.ID] = cluster_data.apply(
        lambda row: f"{ColumnName.RANK}={row[ColumnName.RANK]};{ColumnName.READS}={row[ColumnName.READS]};{ColumnName.RPU}={row[ColumnName.RPU]};{ColumnName.CLUSTER}={row[ColumnName.CLUSTER]};{ColumnName.RANK_IN_CLUSTER}={row[ColumnName.RANK_IN_CLUSTER]};{ColumnName.LED}={row[ColumnName.LED]}", 
        axis=1
    )
    
    # Reorder columns to put ID first
    cols = [ColumnName.ID, ColumnName.RANK, ColumnName.READS,ColumnName.RPU,ColumnName.CLUSTER,ColumnName.RANK_IN_CLUSTER,ColumnName.SEQUENCES]
    cluster_data = cluster_data[cols]

    file_service.save_sequences(cluster_data, output_path, output_format)
    
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