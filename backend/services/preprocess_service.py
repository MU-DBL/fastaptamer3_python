import pandas as pd
from Bio import SeqIO
import regex
import os
from services import file_service 

def run_preprocess(input_path=None, const5p="", const3p="", 
                  length_range=(10, 100), max_error=0.005,
                  output_path=None, output_format='fasta'):
    """
    Preprocess FASTA/FASTQ files with trimming and filtering
    """
    # Get file extension
    if input_path.endswith(('.fq', '.fastq', '.fq.gz', '.fastq.gz')):
        file_format = 'fastq'
    elif input_path.endswith(('.fa', '.fasta', '.fa.gz', '.fasta.gz')):
        file_format = 'fasta'
    else:
        return None
    
    # Read sequences
    records = list(SeqIO.parse(input_path, file_format))
    
    # Create DataFrame
    seq_df = pd.DataFrame({
        'ID': [rec.description for rec in records],
        'Sequence': [str(rec.seq) for rec in records]
    })
    
    # Add quality scores if FASTQ
    if file_format == 'fastq':
        seq_df['Quality'] = [rec.letter_annotations['phred_quality'] for rec in records]
    
    # Trim 5' constant region
    if const5p != "":
        seq_df = trim_constant_region(seq_df, const5p, file_format)
    
    # Trim 3' constant region
    if const3p != "":
        seq_df = trim_constant_region(seq_df, const3p, file_format)
    
    # Length filtering
    seq_df = seq_df[
        (seq_df['Sequence'].str.len() >= length_range[0]) &
        (seq_df['Sequence'].str.len() <= length_range[1])
    ]
    
    # Quality filtering
    if 'Quality' in seq_df.columns:
        seq_df['Avg_Error'] = seq_df['Quality'].apply(
            lambda q: sum([10**(-score/10) for score in q]) / len(q)
        )
        seq_df = seq_df[seq_df['Avg_Error'] <= max_error]
    
    file_service.save_sequences(seq_df, output_path, output_format)
    
    return output_path


def trim_constant_region(seq_df, const_region, file_format):
    """Trim constant region with fuzzy matching"""
    # max.distance = 0.1
    max_distance = int(len(const_region) * 0.1) + 1
    
    trimmed_sequences = []
    trimmed_qualities = []
    
    for idx, row in seq_df.iterrows():
        sequence = row['Sequence']
        
        # Fuzzy search for pattern
        match = regex.search(
            f"({const_region}){{e<={max_distance}}}", 
            sequence
        )
        
        if match:
            # Remove matched region
            start, end = match.span()
            new_seq = sequence[:start] + sequence[end:]
            trimmed_sequences.append(new_seq)
            
            # Handle quality scores
            if file_format == 'fastq':
                quality = row['Quality']
                new_qual = quality[:start] + quality[end:]
                trimmed_qualities.append(new_qual)
        else:
            trimmed_sequences.append(sequence)
            if file_format == 'fastq':
                trimmed_qualities.append(row['Quality'])
    
    seq_df['Sequence'] = trimmed_sequences
    if file_format == 'fastq':
        seq_df['Quality'] = trimmed_qualities
    
    return seq_df
