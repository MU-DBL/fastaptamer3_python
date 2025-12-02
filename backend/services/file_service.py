from pathlib import Path
import re
from typing import Any
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from services.constants import ColumnName

def save_sequences(seq_df, output_path, output_format='fasta'):
    if output_format in ['fastq', 'fasta']:
        records = []
        for idx, row in seq_df.iterrows():
            seq_record = SeqRecord(
                Seq(row[ColumnName.SEQUENCES]),
                id=row[ColumnName.ID].split()[0],
                description=row[ColumnName.ID] 
            )
            
            # Add quality scores for FASTQ output
            if output_format == 'fastq' and 'Quality' in seq_df.columns:
                seq_record.letter_annotations["phred_quality"] = row['Quality']
            
            records.append(seq_record)
        
        # Write to file
        SeqIO.write(records, output_path, output_format)

    elif output_format == 'csv':
        seq_df.to_csv(output_path, index=False) 
    else:
        raise ValueError(f"Unsupported output format: {output_format}. Choose from: 'fasta', 'fastq', 'csv'")


def parse_fasta(fasta_input):
    
    # Lists to store parsed data
    data = []
    
    # Parse FASTA file using Biopython
    for record in SeqIO.parse(fasta_input, "fasta"):
        # Get the header (description)
        header = record.description
        
        # Get the sequence as string
        sequence = str(record.seq)
        
        # Parse header by splitting on semicolons
        parts = header.split(';')
        
        # Create a dictionary to store key-value pairs
        header_dict = {}
        header_dict[ColumnName.ID] = header
        for part in parts:
            if '=' in part:
                key, value = part.split('=', 1)
                key = key.strip()
                value = value.strip()
                
                # Convert to appropriate type based on key
                # print(key)
                if key == ColumnName.RANK:
                    header_dict[ColumnName.RANK] = int(value)
                elif key == ColumnName.READS:
                    header_dict[ColumnName.READS] = int(value)
                elif key == ColumnName.RPU:
                    header_dict[ColumnName.RPU] = float(value)
                elif key == ColumnName.CLUSTER:
                    header_dict[ColumnName.CLUSTER] = int(value)
                elif key == ColumnName.RANK_IN_CLUSTER:
                    header_dict[ColumnName.RANK_IN_CLUSTER] = int(value)
                elif key == ColumnName.LED:
                    header_dict[ColumnName.LED] = int(value)
        
        # Add sequence
        header_dict[ColumnName.SEQUENCES] = sequence
        
        # Add to data list
        data.append(header_dict)
    
    # Create DataFrame
    fasta_df = pd.DataFrame(data)
    print(fasta_df.columns)
    
    # Ensure we have basic required columns
    if ColumnName.RANK not in fasta_df.columns:
        fasta_df[ColumnName.RANK] = range(1, len(fasta_df) + 1)
    if ColumnName.READS not in fasta_df.columns:
        fasta_df[ColumnName.READS] = 0
    if ColumnName.RPU not in fasta_df.columns:
        fasta_df[ColumnName.RPU] = 0.0
    
    # Determine column order based on what's available
    base_cols = [ColumnName.ID, ColumnName.RANK, ColumnName.READS, ColumnName.RPU]
    cluster_cols = []
    if ColumnName.CLUSTER in fasta_df.columns:
        cluster_cols.append(ColumnName.CLUSTER)
    if ColumnName.RANK_IN_CLUSTER in fasta_df.columns:
        cluster_cols.append(ColumnName.RANK_IN_CLUSTER)
    if ColumnName.LED in fasta_df.columns:
        cluster_cols.append(ColumnName.LED)
    
    column_order = base_cols + cluster_cols + [ColumnName.SEQUENCES]
    fasta_df = fasta_df[column_order]
    
    return fasta_df


def read_file(filepath: Path) -> pd.DataFrame:
    if filepath.suffix.lower() == '.csv':
        df = pd.read_csv(filepath)
    elif filepath.suffix.lower() == '.fasta':
        df = parse_fasta(filepath)
    
    if df.empty:
        raise ValueError("file is empty")
    
    return df


