import re
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def save_sequences(seq_df, output_path, output_format='fasta'):
    if output_format in ['fastq', 'fasta']:
        records = []
        for idx, row in seq_df.iterrows():
            seq_record = SeqRecord(
                Seq(row['Sequence']),
                id=row['ID'].split()[0],
                description=row['ID'] 
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

def parse_fasta_with_reads(fasta_file):
    sequences = []
    reads = []
    
    for record in SeqIO.parse(fasta_file, 'fasta'):
        sequences.append(str(record.seq))
        header = record.description
        match = re.search(r';read=(\d+)', header)
        if match:
            count = int(match.group(1))
        elif '-' in record.id:
            parts = record.id.split('-')
            count = int(parts[1])
        else:
            count = 1
        
        reads.append(count)
    
    return pd.DataFrame({'Sequence': sequences, 'Read': reads})
