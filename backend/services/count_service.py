import pandas as pd
import gzip
from collections import Counter
from pathlib import Path
from Bio.Seq import Seq
from services.constants import ColumnName
from services import file_service


def _read_fastq_sequences(handle):
    """Read only sequence lines from FASTQ (every 2nd of 4 lines)."""
    while True:
        header = handle.readline()
        if not header:
            break
        seq = handle.readline().rstrip('\n\r')
        handle.readline()  # +
        handle.readline()  # quality
        yield seq


def _read_fasta_sequences(handle):
    """Read only sequence lines from FASTA (non-header lines)."""
    seq_parts = []
    for line in handle:
        if line.startswith('>'):
            if seq_parts:
                yield ''.join(seq_parts)
                seq_parts = []
        else:
            seq_parts.append(line.rstrip('\n\r'))
    if seq_parts:
        yield ''.join(seq_parts)


def run_count(inputpath=None, reverseComplement=False, scaling_factor=1e6, output_format='fasta', output_path=None):

    file_path = Path(inputpath)
    ext = file_path.suffix.lower()

    opener = gzip.open if ext == '.gz' else open
    inner_ext = file_path.stem.split('.')[-1].lower() if ext == '.gz' else ext.lstrip('.')

    if inner_ext in ('fastq', 'fq'):
        reader = _read_fastq_sequences
    elif inner_ext in ('fasta', 'fa'):
        reader = _read_fasta_sequences
    else:
        return None

    with opener(inputpath, 'rt') as handle:
        counter = Counter(reader(handle))

    seq_counts = pd.DataFrame(counter.most_common(), columns=[ColumnName.SEQUENCES, ColumnName.READS])
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

