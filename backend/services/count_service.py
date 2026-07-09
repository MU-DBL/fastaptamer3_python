import shutil
import subprocess
import pandas as pd
import gzip
from collections import Counter
from pathlib import Path
from services.constants import ColumnName
from services import file_service

_RC_TABLE = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')


def _reverse_complement(seq):
    return seq.translate(_RC_TABLE)[::-1]


def _read_fastq_sequences(handle):
    """Read only sequence lines from FASTQ (every 2nd of 4 lines)."""
    for i, line in enumerate(handle):
        if i % 4 == 1:
            yield line.rstrip('\n\r')


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


def _count_via_sort(inputpath, inner_ext, is_gz):
    """
    Count sequences using sort | uniq -c shell pipeline.
    Uses external merge sort (handles files larger than RAM) at C speed.
    Returns list of (sequence, count) sorted by count descending, or None on failure.
    """
    try:
        # Decompress step
        if is_gz:
            decomp = next((p for p in ('pigz', 'gzip') if shutil.which(p)), None)
            if decomp is None:
                return None
            decompress = f'{decomp} -dc {shutil.quote(inputpath)}'
        else:
            decompress = f'cat {shutil.quote(inputpath)}'

        # Extract sequences only
        if inner_ext in ('fasta', 'fa'):
            extract = f"{decompress} | grep -v '^>'"
        elif inner_ext in ('fastq', 'fq'):
            extract = f"{decompress} | awk 'NR%4==2'"
        else:
            return None

        # Prefer GNU sort (gsort) for --parallel/-S support; fall back to system sort
        sort_bin = 'gsort' if shutil.which('gsort') else 'sort'
        pipeline = f"{extract} | {sort_bin} --parallel=8 -S 2G | uniq -c | sort -rn"

        result = subprocess.run(pipeline, shell=True, capture_output=True, text=True)
        if result.returncode != 0 or not result.stdout:
            return None

        pairs = []
        for line in result.stdout.splitlines():
            line = line.strip()
            if line:
                count_str, seq = line.split(None, 1)
                pairs.append((seq.strip(), int(count_str)))
        return pairs
    except Exception:
        return None


def run_count(inputpath=None, reverseComplement=False, scaling_factor=1e6, output_format='fasta', output_path=None):

    file_path = Path(inputpath)
    ext = file_path.suffix.lower()
    is_gz = ext == '.gz'
    inner_ext = file_path.stem.split('.')[-1].lower() if is_gz else ext.lstrip('.')

    if inner_ext not in ('fasta', 'fa', 'fastq', 'fq'):
        return None

    pairs = _count_via_sort(inputpath, inner_ext, is_gz)

    if pairs is None:
        # Fallback: Python Counter
        opener = gzip.open if is_gz else open
        reader = _read_fastq_sequences if inner_ext in ('fastq', 'fq') else _read_fasta_sequences
        with opener(inputpath, 'rt') as handle:
            counter = Counter(reader(handle))
        pairs = counter.most_common()

    seq_counts = pd.DataFrame(pairs, columns=[ColumnName.SEQUENCES, ColumnName.READS])
    seq_counts[ColumnName.RANK] = seq_counts.index + 1

    total_reads = seq_counts[ColumnName.READS].sum()
    seq_counts[ColumnName.RPU] = (seq_counts[ColumnName.READS] / (total_reads / scaling_factor)).round(0).astype(int)
    seq_counts[ColumnName.LENGTH] = seq_counts[ColumnName.SEQUENCES].str.len()
    seq_counts[ColumnName.ID] = (ColumnName.RANK + '=' + seq_counts[ColumnName.RANK].astype(str) + ';' +
                                  ColumnName.READS + '=' + seq_counts[ColumnName.READS].astype(str) + ';' +
                                  ColumnName.RPU + '=' + seq_counts[ColumnName.RPU].astype(str))

    seq_counts = seq_counts[[ColumnName.ID, ColumnName.RANK, ColumnName.READS, ColumnName.RPU, ColumnName.LENGTH, ColumnName.SEQUENCES]]

    if reverseComplement:
        seq_counts[ColumnName.SEQUENCES] = seq_counts[ColumnName.SEQUENCES].apply(_reverse_complement)

    file_service.save_sequences(seq_counts, output_path, output_format)

    return output_path

