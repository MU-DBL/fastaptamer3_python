import pandas as pd
import gzip
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from services import file_service


def run_count(inputpath=None, reverseComplement=False, scaling_factor=1e6, output_format='fasta', output_path=None):
    
    # Get file extension
    file_path = Path(inputpath)
    ext = file_path.suffix.lower()
    
    # Read sequences based on file type
    sequences = []
    
    if ext == '.gz':
        # Handle gzipped files
        stem_ext = file_path.stem.split('.')[-1].lower()
        
        if stem_ext == 'fastq':
            with gzip.open(inputpath, 'rt') as handle:
                sequences = [str(record.seq) for record in SeqIO.parse(handle, 'fastq')]
        elif stem_ext == 'fasta':
            with gzip.open(inputpath, 'rt') as handle:
                sequences = [str(record.seq) for record in SeqIO.parse(handle, 'fasta')]
        else:
            return None
            
    elif ext == '.fastq':
        with open(inputpath, 'r') as handle:
            sequences = [str(record.seq) for record in SeqIO.parse(handle, 'fastq')]
            
    elif ext in ['.fasta', '.fa']:
        with open(inputpath, 'r') as handle:
            sequences = [str(record.seq) for record in SeqIO.parse(handle, 'fasta')]
            
    else:
        return None

    seq_counts = pd.Series(sequences).value_counts().reset_index()
    seq_counts.columns = ['Sequence', 'Read']
    
    # Sort by reads (descending) and add rank
    seq_counts = seq_counts.sort_values('Read', ascending=False).reset_index(drop=True)
    seq_counts['Rank'] = seq_counts.index + 1

    total_reads = seq_counts['Read'].sum()
    seq_counts['RPU'] = (seq_counts['Read'] / (total_reads / scaling_factor)).round(0).astype(int)
    seq_counts['Length'] = seq_counts['Sequence'].str.len()
    seq_counts['ID'] = 'rank=' + seq_counts['Rank'].astype(str) + ';read=' + seq_counts['Read'].astype(str) + ';RPU=' + seq_counts['RPU'].astype(str)
    seq_counts['Sequence'] = seq_counts['Sequence'].str.replace('\r', '', regex=False)
    
    # Reorder columns
    seq_counts = seq_counts[['ID', 'Rank', 'Read', 'RPU', 'Length', 'Sequence']]
    
    # Optionally make reverse complement
    if reverseComplement:
        seq_counts['Sequence'] = seq_counts['Sequence'].apply(
            lambda x: str(Seq(x).reverse_complement())
        )
    
    file_service.save_sequences(seq_counts, output_path, output_format)
    
    return output_path

