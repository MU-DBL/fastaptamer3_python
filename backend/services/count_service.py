import pandas as pd
import gzip
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from services.constants import ColumnName
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
    seq_counts.columns = [ColumnName.SEQUENCES, ColumnName.READS]
    
    # Sort by reads (descending) and add rank
    seq_counts = seq_counts.sort_values(ColumnName.READS, ascending=False).reset_index(drop=True)
    seq_counts[ColumnName.RANK] = seq_counts.index + 1

    total_reads = seq_counts[ColumnName.READS].sum()
    seq_counts[ColumnName.RPU] = (seq_counts[ColumnName.READS] / (total_reads / scaling_factor)).round(0).astype(int)
    seq_counts[ColumnName.LENGTH] = seq_counts[ColumnName.SEQUENCES].str.len()
    seq_counts[ColumnName.ID] = ColumnName.RANK + '=' + seq_counts[ColumnName.RANK].astype(str) + ';' + ColumnName.READS + '=' + seq_counts[ColumnName.READS].astype(str) + ';' + ColumnName.RPU + '=' + seq_counts[ColumnName.RPU].astype(str)
    seq_counts[ColumnName.SEQUENCES] = seq_counts[ColumnName.SEQUENCES].str.replace('\r', '', regex=False)
    
    # Reorder columns
    seq_counts = seq_counts[[ColumnName.ID, ColumnName.RANK, ColumnName.READS, ColumnName.RPU, ColumnName.LENGTH, ColumnName.SEQUENCES]]
    
    # Optionally make reverse complement
    if reverseComplement:
        seq_counts[ColumnName.SEQUENCES] = seq_counts[ColumnName.SEQUENCES].apply(
            lambda x: str(Seq(x).reverse_complement())
        )
    
    file_service.save_sequences(seq_counts, output_path, output_format)
    
    return output_path

