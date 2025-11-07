"""
K-mer analysis service for sequence similarity analysis.
Implements the R workflow from clusterDiversity_plots.R
"""

import numpy as np
import pandas as pd
from collections import defaultdict, Counter
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import umap
from services.file_service import parse_fasta
from services.constants import ColumnName
import logging

logger = logging.getLogger(__name__)

def generate_kmer_counts(sequences, k=4):
    """
    Generate k-mer counts for a list of sequences.
    Equivalent to R's kmer::kcount() function.
    
    Args:
        sequences (list): List of DNA sequences
        k (int): K-mer size (default: 4)
    
    Returns:
        pd.DataFrame: K-mer count matrix with sequences as rows and k-mers as columns
    """
    # Generate all possible k-mers for DNA sequences
    nucleotides = ['A', 'T', 'G', 'C']
    
    def generate_all_kmers(k):
        """Generate all possible k-mers of length k"""
        if k == 1:
            return nucleotides
        else:
            smaller_kmers = generate_all_kmers(k - 1)
            return [kmer + nucleotide for kmer in smaller_kmers for nucleotide in nucleotides]
    
    all_kmers = generate_all_kmers(k)
    
    # Initialize k-mer count matrix
    kmer_matrix = []
    
    for sequence in sequences:
        sequence = sequence.upper()
        kmer_counts = Counter()
        
        # Count k-mers in this sequence
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            if all(nucleotide in nucleotides for nucleotide in kmer):  # Valid DNA k-mer
                kmer_counts[kmer] += 1
        
        # Create row with counts for all possible k-mers
        row = [kmer_counts.get(kmer, 0) for kmer in all_kmers]
        kmer_matrix.append(row)
    
    return pd.DataFrame(kmer_matrix, columns=all_kmers)

def apply_pca(kmer_matrix, n_components=2):
    """
    Apply PCA to k-mer matrix.
    Equivalent to R's stats::prcomp() function.
    
    Args:
        kmer_matrix (pd.DataFrame): K-mer count matrix
        n_components (int): Number of principal components (default: 2)
    
    Returns:
        tuple: (transformed_data, pca_object, explained_variance_ratio)
    """
    # Standardize the data (R's prcomp with scale=TRUE)
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(kmer_matrix)
    
    # Apply PCA
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(scaled_data)
    
    return pca_result, pca, pca.explained_variance_ratio_

def apply_umap(kmer_matrix, n_components=2, n_neighbors=15, min_dist=0.1):
    """
    Apply UMAP to k-mer matrix.
    Equivalent to R's umap::umap() function.
    
    Args:
        kmer_matrix (pd.DataFrame): K-mer count matrix
        n_components (int): Number of dimensions for UMAP (default: 2)
        n_neighbors (int): Number of neighbors for UMAP (default: 15)
        min_dist (float): Minimum distance for UMAP (default: 0.1)
    
    Returns:
        np.ndarray: UMAP transformed coordinates
    """
    # Standardize the data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(kmer_matrix)
    
    # Apply UMAP
    reducer = umap.UMAP(
        n_components=n_components,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        random_state=42  # For reproducible results
    )
    umap_result = reducer.fit_transform(scaled_data)
    
    return umap_result



def analyze_kmer_diversity(input_path, selected_clusters, k_size=4, method='pca'):
    """
    Perform k-mer analysis on selected clusters following R implementation workflow.
    
    This function implements the exact workflow from clusterDiversity_plots.R:
    1. Filter sequences by selected clusters
    2. Generate k-mer count matrix using kmer::kcount equivalent
    3. Apply dimensionality reduction (PCA or UMAP)
    4. Prepare coordinates with cluster assignments
    
    Args:
        input_path (str): Path to clustered FASTA file
        selected_clusters (list): List of cluster numbers to analyze
        k_size (int): K-mer size for analysis (default: 4)
        method (str): Dimensionality reduction method ('pca' or 'umap')
    
    Returns:
        dict: Analysis results with coordinates, clusters, and metadata
    """
    try:
        logger.info(f"Starting k-mer analysis for clusters: {selected_clusters}")
        
        # Step 1: Read and filter clustered FASTA file
        df = parse_fasta(input_path)
        
        # Verify required columns
        required_cols = [ColumnName.CLUSTER, ColumnName.SEQUENCES]
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Input file missing required columns: {missing_cols}")
        
        # Filter by selected clusters
        filtered_df = df[df[ColumnName.CLUSTER].isin(selected_clusters)].copy()
        
        if len(filtered_df) == 0:
            raise ValueError(f"No sequences found for selected clusters: {selected_clusters}")
        
        logger.info(f"Filtered to {len(filtered_df)} sequences from {len(selected_clusters)} clusters")
        
        # Step 2: Generate k-mer count matrix (equivalent to R's kmer::kcount)
        sequences = filtered_df[ColumnName.SEQUENCES].tolist()
        kmer_matrix = generate_kmer_counts(sequences, k=k_size)
        
        logger.info(f"Generated k-mer matrix: {kmer_matrix.shape[0]} sequences x {kmer_matrix.shape[1]} k-mers")
        
        # Step 3: Apply dimensionality reduction
        if method.lower() == 'pca':
            # Apply PCA (equivalent to R's stats::prcomp)
            coordinates, pca_obj, explained_var = apply_pca(kmer_matrix, n_components=2)
            method_info = {
                'method': 'PCA',
                'explained_variance': {
                    'PC1': round(explained_var[0] * 100, 2),
                    'PC2': round(explained_var[1] * 100, 2)
                }
            }
            axis_labels = {
                'x': f"PC1 ({method_info['explained_variance']['PC1']}%)",
                'y': f"PC2 ({method_info['explained_variance']['PC2']}%)"
            }
        elif method.lower() == 'umap':
            # Apply UMAP (equivalent to R's umap::umap)
            coordinates = apply_umap(kmer_matrix, n_components=2)
            method_info = {'method': 'UMAP'}
            axis_labels = {'x': 'UMAP1', 'y': 'UMAP2'}
        else:
            raise ValueError(f"Unsupported method: {method}. Use 'pca' or 'umap'")
        
        # Step 4: Prepare results (colors handled by frontend)
        unique_clusters = sorted(filtered_df[ColumnName.CLUSTER].unique())
        
        # Prepare coordinate data
        result_data = []
        for i, (x, y) in enumerate(coordinates):
            cluster_id = filtered_df.iloc[i][ColumnName.CLUSTER]
            sequence = filtered_df.iloc[i][ColumnName.SEQUENCES]
            
            result_data.append({
                'x': float(x),
                'y': float(y),
                'cluster': int(cluster_id),
                'sequence': sequence
            })
        
        # Prepare cluster summary
        cluster_summary = []
        for cluster_id in unique_clusters:
            cluster_sequences = filtered_df[filtered_df[ColumnName.CLUSTER] == cluster_id]
            cluster_summary.append({
                'cluster': int(cluster_id),
                'count': len(cluster_sequences)
            })
        
        result = {
            'coordinates': result_data,
            'clusters': cluster_summary,
            'method_info': method_info,
            'axis_labels': axis_labels,
            'parameters': {
                'k_size': k_size,
                'method': method,
                'total_sequences': len(result_data),
                'total_clusters': len(unique_clusters)
            }
        }
        
        logger.info(f"K-mer analysis completed successfully")
        return result
        
    except Exception as e:
        logger.error(f"K-mer analysis failed: {str(e)}")
        raise