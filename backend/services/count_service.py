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


def run_count(inputpath=None, reverseComplement=False, scaling_factor=1e6, output_format='fasta', output_path=None):

    file_path = Path(inputpath)
    ext = file_path.suffix.lower()
    is_gz = ext == '.gz'
    inner_ext = file_path.stem.split('.')[-1].lower() if is_gz else ext.lstrip('.')

    if inner_ext not in ('fasta', 'fa', 'fastq', 'fq'):
        return None

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


def _iter_count_records(result_path):
    """
    Stream (rank, reads, rpu, length) tuples from a saved count result file
    (CSV or FASTA), without loading the whole file into memory. Records are
    yielded in rank order since that's how run_count wrote them.
    """
    ext = Path(result_path).suffix.lower()

    if ext == '.csv':
        with open(result_path, 'r') as f:
            header = f.readline().rstrip('\n').split(',')
            idx = {name: i for i, name in enumerate(header)}
            for line in f:
                values = line.rstrip('\n').split(',')
                yield (
                    int(values[idx[ColumnName.RANK]]),
                    int(values[idx[ColumnName.READS]]),
                    float(values[idx[ColumnName.RPU]]),
                    int(values[idx[ColumnName.LENGTH]]),
                )
    else:  # fasta
        with open(result_path, 'r') as f:
            rank = reads = rpu = None
            seq_len = 0
            for line in f:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    if rank is not None:
                        yield (rank, reads, rpu, seq_len)
                    parts = dict(p.split('=', 1) for p in line[1:].split(';') if '=' in p)
                    rank = int(parts[ColumnName.RANK])
                    reads = int(parts[ColumnName.READS])
                    rpu = float(parts[ColumnName.RPU])
                    seq_len = 0
                else:
                    seq_len += len(line)
            if rank is not None:
                yield (rank, reads, rpu, seq_len)


def get_reads_per_rank(result_path, min_reads=0, max_rank=None, metric='reads', max_points=3000):
    """
    Reads/RPU per rank, filtered and (if needed) log-spaced downsampled for plotting.
    Stops reading early once max_rank is passed, since records are rank-ordered.
    """
    points = []
    for rank, reads, rpu, _length in _iter_count_records(result_path):
        if max_rank is not None and rank > max_rank:
            break
        value = rpu if metric == 'rpu' else reads
        if value >= min_reads:
            points.append((rank, value))

    if len(points) > max_points:
        step = len(points) / max_points
        idxs = sorted(set(int(i * step) for i in range(max_points)) | {0, len(points) - 1})
        points = [points[i] for i in idxs]

    return points


def _weighted_median(values_and_weights):
    """Median of `values`, each counted `weight` times, without materializing the
    expanded list."""
    total = sum(w for _v, w in values_and_weights)
    if total == 0:
        return 0
    half = total / 2
    cum = 0
    for v, w in sorted(values_and_weights):
        cum += w
        if cum >= half:
            return v
    return values_and_weights[-1][0]


def get_sequence_length_histogram(result_path, outlier_multiple=10):
    """
    Unique-sequence count and total-reads count, grouped by sequence length.

    A malformed record (e.g. a missing FASTA header merging several sequence
    lines together) can produce a length wildly outside the real distribution,
    which then dominates the histogram's axis and makes the real data
    invisible. Lengths more than `outlier_multiple`x the (unique-count-weighted)
    median are reported separately instead of in the main histogram.
    """
    unique_by_length = {}
    reads_by_length = {}
    for _rank, reads, _rpu, length in _iter_count_records(result_path):
        unique_by_length[length] = unique_by_length.get(length, 0) + 1
        reads_by_length[length] = reads_by_length.get(length, 0) + reads

    median_length = _weighted_median(list(unique_by_length.items())) or 1
    threshold = median_length * outlier_multiple

    lengths = sorted(l for l in unique_by_length if l <= threshold)
    outlier_lengths = sorted(l for l in unique_by_length if l > threshold)

    return {
        'lengths': lengths,
        'unique': [unique_by_length[l] for l in lengths],
        'reads': [reads_by_length[l] for l in lengths],
        'excluded_outliers': [
            {'length': l, 'unique': unique_by_length[l], 'reads': reads_by_length[l]}
            for l in outlier_lengths
        ],
    }


def get_abundance_bins(result_path, breakpoints=(10, 100, 1000), use_singleton=True):
    """Bin sequences by read count, matching the FA3 abundance-plot convention."""
    breakpoints = sorted(breakpoints)
    bins = {}

    def bin_label(reads):
        if use_singleton and reads == 1:
            return 'Singleton'
        for i, b in enumerate(breakpoints):
            if i == 0 and reads < b:
                return f'1 < Reads < {b}'
            if reads < b:
                return f'{breakpoints[i - 1]} <= Reads < {b}'
        return f'{breakpoints[-1]} <= Reads'

    for _rank, reads, _rpu, _length in _iter_count_records(result_path):
        label = bin_label(reads)
        entry = bins.setdefault(label, {'count': 0, 'total_reads': 0})
        entry['count'] += 1
        entry['total_reads'] += reads

    total_reads = sum(b['total_reads'] for b in bins.values()) or 1
    return [
        {'bin': label, 'fraction': b['total_reads'] / total_reads, 'unique_count': b['count']}
        for label, b in bins.items()
    ]


def _iter_count_rows(result_path):
    """Stream full row dicts (id, rank, reads, rpm, length, seqs) from a saved
    count result file, without loading the whole file into memory."""
    ext = Path(result_path).suffix.lower()

    if ext == '.csv':
        with open(result_path, 'r') as f:
            header = f.readline().rstrip('\n').split(',')
            idx = {name: i for i, name in enumerate(header)}
            for line in f:
                values = line.rstrip('\n').split(',')
                yield {
                    'id': values[idx[ColumnName.ID]],
                    'rank': int(values[idx[ColumnName.RANK]]),
                    'reads': int(values[idx[ColumnName.READS]]),
                    'rpm': int(float(values[idx[ColumnName.RPU]])),
                    'length': int(values[idx[ColumnName.LENGTH]]),
                    'seqs': values[idx[ColumnName.SEQUENCES]],
                }
    else:  # fasta
        with open(result_path, 'r') as f:
            header_line = None
            seq_parts = []
            for line in f:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    if header_line is not None:
                        yield _fasta_row(header_line, ''.join(seq_parts))
                    header_line = line[1:]
                    seq_parts = []
                else:
                    seq_parts.append(line)
            if header_line is not None:
                yield _fasta_row(header_line, ''.join(seq_parts))


def _fasta_row(header_line, seq):
    parts = dict(p.split('=', 1) for p in header_line.split(';') if '=' in p)
    return {
        'id': header_line,
        'rank': int(parts[ColumnName.RANK]),
        'reads': int(parts[ColumnName.READS]),
        'rpm': int(float(parts[ColumnName.RPU])),
        'length': len(seq),
        'seqs': seq,
    }


def get_table_preview(result_path, limit=50000):
    """
    First `limit` rows for the table, plus totals accumulated over the whole
    file - all in a single streaming pass, without ever holding the full file
    (or a full-size row array) in memory.
    """
    rows = []
    total_records = 0
    total_reads = 0

    for row in _iter_count_rows(result_path):
        total_records += 1
        total_reads += row['reads']
        if total_records <= limit:
            rows.append(row)

    return {
        'rows': rows,
        'total_records': total_records,
        'total_reads': total_reads,
        'truncated': total_records > limit,
    }

