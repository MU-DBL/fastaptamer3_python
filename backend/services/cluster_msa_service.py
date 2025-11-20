import pandas as pd
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import tempfile
import os
from Bio.Align.Applications import MuscleCommandline
from services.constants import ColumnName
from services.file_service import parse_fasta,save_sequences

def cluster_msa(df_cluster, output_path, output_format, cluster_selection=1, seq_type="dna"):
    
    if df_cluster is None:
        raise ValueError("df_cluster cannot be None")
    
    # Filter for the selected cluster
    cluster_df = df_cluster[df_cluster[ColumnName.CLUSTER] == cluster_selection].copy()
    
    if len(cluster_df) == 0:
        raise ValueError(f"No sequences found for cluster {cluster_selection}")
    
    # If only one sequence, no alignment needed
    if len(cluster_df) == 1:
        # Relocate Sequences column to the end
        cols = [col for col in cluster_df.columns if col != ColumnName.SEQUENCES]
        cols.append(ColumnName.SEQUENCES)
        cluster_df = cluster_df[cols]
        
        if ColumnName.RANK in cluster_df.columns:
            cluster_df = cluster_df.sort_values(ColumnName.RANK).reset_index(drop=True)
        
        return cluster_df
    
    # Create SeqRecord objects for MSA
    seq_records = [
        SeqRecord(Seq(row[ColumnName.SEQUENCES]), id=row[ColumnName.ID], description="")
        for _, row in cluster_df.iterrows()
    ]
    
    # Perform MSA using MUSCLE
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as input_file:
        input_path = input_file.name
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as output_file:
        output_path = output_file.name
    
    try:
        # Write sequences to input file
        SeqIO.write(seq_records, input_path, "fasta")
        
        # Run MUSCLE
        muscle_cline = MuscleCommandline(
            input=input_path,
            out=output_path
        )
        
        stdout, stderr = muscle_cline()
        
        # Read aligned sequences
        alignment = AlignIO.read(output_path, "fasta")
        
        # Convert to dictionary maintaining order
        msa_dict = {record.id: str(record.seq) for record in alignment}
        
    except Exception as e:
        raise RuntimeError(f"MUSCLE alignment failed for cluster {cluster_selection}: {str(e)}")
    
    finally:
        # Clean up temporary files
        if os.path.exists(input_path):
            os.remove(input_path)
        if os.path.exists(output_path):
            os.remove(output_path)
    
    # Create MSA dataframe
    msa_df = pd.DataFrame({
        ColumnName.ID: list(msa_dict.keys()),
        ColumnName.SEQUENCES: list(msa_dict.values())
    })
    
    # Join with original metadata (excluding original Sequences column)
    metadata_cols = [col for col in df_cluster.columns if col != ColumnName.SEQUENCES]
    metadata_df = df_cluster[metadata_cols]
    
    msa_df = msa_df.merge(metadata_df, on=ColumnName.ID, how='left')
    
    # Relocate Sequences column to the end
    cols = [col for col in msa_df.columns if col != ColumnName.SEQUENCES]
    cols.append(ColumnName.SEQUENCES)
    msa_df = msa_df[cols]
    
    # Sort by Rank if it exists
    if ColumnName.RANK in msa_df.columns:
        msa_df = msa_df.sort_values(ColumnName.RANK).reset_index(drop=True)
    
    save_sequences(msa_df, output_path, output_format)

    return output_path
