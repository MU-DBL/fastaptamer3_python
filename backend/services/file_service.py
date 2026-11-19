from pathlib import Path
import csv
import io
import re
from typing import Any
import numpy as np
import pandas as pd
from Bio import SeqIO
from services.constants import ColumnName

_FASTA_WRAP = 60
# Records per write, so large outputs are never held in memory as one string
_WRITE_CHUNK = 100_000
# Decimals kept for RPU (as in FASTAptameR 2.0's RPM)
RPU_DECIMALS = 2


def _wrap_seq(seq, width=_FASTA_WRAP):
    if len(seq) <= width:
        return seq
    return '\n'.join(seq[i:i + width] for i in range(0, len(seq), width))


def _sanger_quality_str(qualities):
    return ''.join(chr(min(126, int(round(q)) + 33)) for q in qualities)


def save_sequences(seq_df, output_path, output_format='fasta'):
    # Bypasses Biopython's per-record SeqRecord/SeqIO.write machinery, which is
    # extremely slow (minutes) once seq_df reaches millions of rows.
    if output_format == 'fasta':
        ids = seq_df[ColumnName.ID].to_numpy()
        seqs = seq_df[ColumnName.SEQUENCES].to_numpy()
        with open(output_path, 'w') as f:
            for start in range(0, len(ids), _WRITE_CHUNK):
                end = start + _WRITE_CHUNK
                f.write(''.join(f'>{header}\n{_wrap_seq(seq)}\n' for header, seq in zip(ids[start:end], seqs[start:end])))

    elif output_format == 'fastq':
        if 'Quality' not in seq_df.columns:
            raise ValueError("No suitable quality scores found in data for FASTQ output")
        ids = seq_df[ColumnName.ID].to_numpy()
        seqs = seq_df[ColumnName.SEQUENCES].to_numpy()
        quals = seq_df['Quality'].to_numpy()
        with open(output_path, 'w') as f:
            for start in range(0, len(ids), _WRITE_CHUNK):
                end = start + _WRITE_CHUNK
                f.write(''.join(
                    f'@{header}\n{seq}\n+\n{_sanger_quality_str(qual)}\n'
                    for header, seq, qual in zip(ids[start:end], seqs[start:end], quals[start:end])
                ))

    elif output_format == 'csv':
        seq_df.to_csv(output_path, index=False)
    else:
        raise ValueError(f"Unsupported output format: {output_format}. Choose from: 'fasta', 'fastq', 'csv'")


def compute_rpu(reads, scaling_factor):
    """RPU = reads per U of all reads, kept to RPU_DECIMALS so low-abundance sequences do not round to 0."""
    return (reads / (reads.sum() / scaling_factor)).round(RPU_DECIMALS)


def format_rpu(rpu):
    """Render RPU values for sequence headers with a fixed number of decimals (e.g. 24106.50, 0.16)."""
    return rpu.map(lambda v: f'{v:.{RPU_DECIMALS}f}')


def format_scaling_factor(scaling_factor):
    """Render the RPU scaling factor U for headers, e.g. 1e6 -> '1000000'."""
    value = float(scaling_factor)
    return str(int(value)) if value.is_integer() else str(value)


def scaling_factor_of(*columns):
    """
    The U recorded in the given U columns, formatted for output, or None if none
    records it (e.g. files from before U was recorded). Callers check beforehand
    that the inputs share one U, so the first recorded value stands for all.
    """
    for column in columns:
        values = column.dropna()
        if len(values):
            return format_scaling_factor(values.iloc[0])
    return None


def check_same_scaling_factor(dfs, names):
    """
    Raise ValueError if the inputs' RPU were computed with different U, since
    comparing RPU across them would then be off by the ratio of the U values.
    Inputs without U (e.g. files from before U was recorded) are skipped.
    """
    u_by_file = {}
    for df, name in zip(dfs, names):
        if ColumnName.U in df.columns:
            u_by_file[name] = set(df[ColumnName.U].dropna().astype(float))

    all_u = set().union(*u_by_file.values()) if u_by_file else set()
    if len(all_u) > 1:
        details = ", ".join(
            f"{name}: U={'/'.join(format_scaling_factor(u) for u in sorted(us))}"
            for name, us in u_by_file.items()
        )
        raise ValueError(
            f"Input files use different RPU scaling factors ({details}). "
            f"Re-run Count with the same U for all files before comparing RPU."
        )


# Header keys parsed into columns, with their types
_HEADER_TYPES = {
    ColumnName.RANK: int, ColumnName.READS: int, ColumnName.RPU: float, ColumnName.U: str,
    ColumnName.CLUSTER: int, ColumnName.RANK_IN_CLUSTER: int, ColumnName.LED: int, ColumnName.DISTANCE: int,
}


def _parse_header(header):
    """Parse one 'Key=value;Key=value' header into a dict (unknown keys and parts without '=' are ignored)."""
    header_dict = {ColumnName.ID: header}
    for part in header.split(';'):
        if '=' in part:
            key, value = part.split('=', 1)
            key, value = key.strip(), value.strip()
            if key in _HEADER_TYPES:
                header_dict[key] = _HEADER_TYPES[key](value)
    return header_dict


def _parse_headers_uniform(headers):
    """Parse all headers in one C-level read_csv pass when they share the first header's key layout.
    Returns None if any header differs (the caller then parses headers one by one)."""
    parts = headers[0].split(';')
    keys = [p.split('=', 1)[0] for p in parts]
    if not all('=' in p for p in parts) or len(set(keys)) != len(keys) or any(k != k.strip() for k in keys):
        return None

    # 'Rank=1;Reads=5' -> 'Rank;1;Reads;5': keys in even columns, values in odd columns
    dtypes = {
        2 * i + 1: np.int64 if _HEADER_TYPES.get(k) is int else np.float64 if _HEADER_TYPES.get(k) is float else str
        for i, k in enumerate(keys)
    }
    try:
        kv = pd.read_csv(
            io.StringIO('\n'.join(headers).replace('=', ';')), sep=';', header=None, dtype=dtypes,
            engine='c', quoting=csv.QUOTE_NONE, na_filter=False, skip_blank_lines=False,
        )
    except (ValueError, pd.errors.ParserError):
        return None
    if kv.shape != (len(headers), 2 * len(keys)) or not all((kv[2 * i] == k).all() for i, k in enumerate(keys)):
        return None

    fasta_df = pd.DataFrame({ColumnName.ID: headers})
    for i, k in enumerate(keys):
        if k in _HEADER_TYPES:
            fasta_df[k] = kv[2 * i + 1].str.strip() if _HEADER_TYPES[k] is str else kv[2 * i + 1].to_numpy()
    return fasta_df


def parse_fasta(fasta_input):
    # Read the whole file and split it into records (much faster than Biopython's per-record SeqRecord objects)
    with open(fasta_input) as f:
        text = f.read()
    start = text.find('>')
    if text and start != 0:
        # Same rule as Biopython's 'fasta' parser: nothing (not even blank lines) before the first record
        raise ValueError("This FASTA file contains comments at the beginning of the file, which are not allowed")
    # Tabs are rare: then strip each line like Biopython instead of the fast whole-record replace
    strip_lines = '\t' in text
    records = text[1:].split('\n>') if text else []
    del text

    # Header (without '>') and sequence (lines joined, whitespace removed) of each record
    headers, sequences = [], []
    for record in records:
        header, _, sequence = record.partition('\n')
        headers.append(header.rstrip())
        if strip_lines:
            sequence = ''.join(line.rstrip() for line in sequence.split('\n'))
        sequences.append(sequence.replace('\n', '').replace('\r', '').replace(' ', ''))
    del records

    # Parse headers into typed columns: all at once if they share one layout, otherwise one by one
    fasta_df = _parse_headers_uniform(headers) if headers else None
    if fasta_df is None:
        fasta_df = pd.DataFrame([_parse_header(h) for h in headers])
    fasta_df[ColumnName.SEQUENCES] = sequences

    # Ensure we have basic required columns
    if ColumnName.RANK not in fasta_df.columns:
        fasta_df[ColumnName.RANK] = range(1, len(fasta_df) + 1)
    if ColumnName.READS not in fasta_df.columns:
        fasta_df[ColumnName.READS] = 0
    if ColumnName.RPU not in fasta_df.columns:
        fasta_df[ColumnName.RPU] = 0.0
    
    # Determine column order based on what's available
    base_cols = [ColumnName.ID, ColumnName.RANK, ColumnName.READS, ColumnName.RPU]
    if ColumnName.U in fasta_df.columns:
        base_cols.append(ColumnName.U)
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


