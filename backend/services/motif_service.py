import re
import regex
import pandas as pd
from services.constants import ColumnName
from services import file_service


def _fuzzy_match_series(series: pd.Series, pattern: str, max_mismatches: int) -> pd.Series:
    """Apply fuzzy (approximate) matching to a pandas Series of sequences."""
    fuzzy_pattern = regex.compile(f"(?:{pattern}){{s<={max_mismatches}}}")
    return series.apply(lambda seq: bool(fuzzy_pattern.search(seq)))


def format_motif(motif_str: str, motif_type: str = "Nucleotide") -> str:
    """
    Format a single motif as a Python regular expression.
    
    Args:
        motif_str: A single motif (e.g., "ACT" or "RYN")
        motif_type: One of "Nucleotide" (default), "AminoAcid", or "String"
                   When "Nucleotide" is selected, interprets nucleotide ambiguity codes
    
    Returns:
        A formatted regex pattern for motif matching
    """
    # Remove spaces
    motif = motif_str.replace(" ", "")
    
    if motif_type == "Nucleotide":
        # For nucleotides: convert to uppercase and handle degenerate codes
        motif = motif.upper()
        
        # Convert U to T for DNA sequences
        motif = motif.replace("U", "T")
        
        # Handle degenerate nucleotide codes
        motif = (motif
                 .replace("R", "[AG]")      # puRine
                 .replace("Y", "[CT]")      # pYrimidine
                 .replace("W", "[AT]")      # Weak
                 .replace("S", "[GC]")      # Strong
                 .replace("M", "[AC]")      # aMino
                 .replace("K", "[GT]")      # Keto
                 .replace("B", "[CGT]")     # not A (B comes after A)
                 .replace("D", "[AGT]")     # not C (D comes after C)
                 .replace("H", "[ACT]")     # not G (H comes after G)
                 .replace("V", "[ACG]")     # not T (V comes after T and U)
                 .replace("N", "[ACGT]"))   # aNy
    elif motif_type == "AminoAcid":
        # For amino acids: convert to uppercase, no special codes
        motif = motif.upper()
    # For "String" type: keep original case, no transformations
    return motif


def search_motif(
    fasta_input: str,
    motif: str,
    highlight: bool = False,
    partial: bool = False,
    motif_type: str = "Nucleotide",
    max_mismatches: int = 0,
    output_format: str = "fasta",
    output_path: str = None
) -> str:
    """
    Search for sequences containing user-defined motifs.

    Args:
        fasta_input: Path to input FASTA file
        motif: A comma-separated list of motifs (e.g., "ACT,AAA,ACG")
        highlight: Whether motifs should be placed inside parentheses in output
        partial: Whether partial matches are allowed (True: OR operation, False: AND operation)
        motif_type: Type of motif - "Nucleotide", "AminoAcid", or "String"
        max_mismatches: Maximum number of allowed substitutions (0 = exact match)
        output_format: Output format - "fasta" or "csv"
        output_path: Path for output file

    Returns:
        Path to the output file
    """
    # Parse input FASTA file
    seq_df = file_service.parse_fasta(fasta_input)
    
    # Split motifs and format each one
    motif_list = [m.strip() for m in motif.split(',') if m.strip()]
    
    # Format each motif pattern individually
    formatted_patterns = []
    for m in motif_list:
        formatted_patterns.append(format_motif(m, motif_type))
    
    # Filter sequences based on motif matching
    # Note: Using case-sensitive matching to match R behavior (grepl default)
    # format_motif already handles uppercase conversion for Nucleotide/AminoAcid
    if max_mismatches > 0:
        # Fuzzy matching using the regex module with substitution tolerance
        if partial:
            # OR: any motif fuzzy-matches
            mask = pd.Series([False] * len(seq_df), index=seq_df.index)
            for pattern in formatted_patterns:
                mask = mask | _fuzzy_match_series(seq_df[ColumnName.SEQUENCES], pattern, max_mismatches)
        else:
            # AND: all motifs must fuzzy-match
            mask = pd.Series([True] * len(seq_df), index=seq_df.index)
            for pattern in formatted_patterns:
                mask = mask & _fuzzy_match_series(seq_df[ColumnName.SEQUENCES], pattern, max_mismatches)
    elif partial:
        # Partial filter uses OR operation (any motif matches)
        combined_pattern = "|".join(formatted_patterns)
        print(f"Combined pattern for partial search: {combined_pattern}")
        mask = seq_df[ColumnName.SEQUENCES].str.contains(
            combined_pattern,
            regex=True,
            case=True
        )
    else:
        # Full filter requires ALL motifs to be present (AND operation)
        mask = pd.Series([True] * len(seq_df), index=seq_df.index)
        for pattern in formatted_patterns:
            pattern_mask = seq_df[ColumnName.SEQUENCES].str.contains(
                pattern,
                regex=True,
                case=True
            )
            mask = mask & pattern_mask

    # Apply filter
    filtered_df = seq_df[mask].copy()

    # Highlight motifs if requested by adding parentheses around matches
    if highlight:
        def highlight_sequence(seq):
            result = seq
            for pattern in formatted_patterns:
                if max_mismatches > 0:
                    fuzzy_pat = regex.compile(f"(?:{pattern}){{s<={max_mismatches}}}")
                    result = fuzzy_pat.sub(lambda m: f"({m.group(0)})", result)
                else:
                    result = re.sub(f"({pattern})", r"(\1)", result)
            return result

        filtered_df[ColumnName.SEQUENCES] = filtered_df[ColumnName.SEQUENCES].apply(highlight_sequence)
    
    # Ensure ID column is properly formatted (following count_service pattern)
    # This creates consistent ID format: Rank=X;Reads=Y;RPU=Z for both CSV and FASTA
    filtered_df[ColumnName.ID] = (
        'Rank=' + filtered_df[ColumnName.RANK].astype(str) + ';' +
        'Reads=' + filtered_df[ColumnName.READS].astype(str) + ';' +
        'RPU=' + filtered_df[ColumnName.RPU].astype(str)
    )
    
    # Save results
    if output_format == 'fasta':
        # Write FASTA manually to prevent BioPython's 60-character line wrapping
        # which breaks frontend parsing
        with open(output_path, 'w') as f:
            for _, row in filtered_df.iterrows():
                f.write(f">{row[ColumnName.ID]}\n")
                f.write(f"{row[ColumnName.SEQUENCES]}\n")
    else:
        # For CSV, use standard save_sequences
        file_service.save_sequences(filtered_df, output_path, output_format)
    
    return output_path


def omit_motif(
    fasta_input: str,
    motif: str,
    partial: bool = False,
    motif_type: str = "Nucleotide",
    max_mismatches: int = 0,
    output_format: str = "fasta",
    output_path: str = None
) -> str:
    """
    Omit sequences containing user-defined motifs.

    Args:
        fasta_input: Path to input FASTA file
        motif: A comma-separated list of motifs (e.g., "ACT,AAA,ACG")
        partial: When True (Yes), omits sequences with at least one motif (OR operation - more aggressive)
                 When False (No), omits only sequences with ALL motifs (AND operation - less aggressive)
        motif_type: Type of motif - "Nucleotide", "AminoAcid", or "String"
        max_mismatches: Maximum number of allowed substitutions (0 = exact match)
        output_format: Output format - "fasta" or "csv"
        output_path: Path for output file

    Returns:
        Path to the output file
    """
    # Parse input FASTA file
    seq_df = file_service.parse_fasta(fasta_input)
    
    # Split motifs and format each one
    motif_list = [m.strip() for m in motif.split(',') if m.strip()]
    
    # Format each motif pattern individually
    formatted_patterns = []
    for m in motif_list:
        formatted_patterns.append(format_motif(m, motif_type))
    
    # Filter sequences based on motif matching
    # Note: Using case-sensitive matching to match R behavior (grepl default)
    # format_motif already handles uppercase conversion for Nucleotide/AminoAcid
    # OMIT logic: opposite of search - we keep sequences that DON'T match
    if max_mismatches > 0:
        # Fuzzy matching using the regex module with substitution tolerance
        if partial:
            # OR: omit if any motif fuzzy-matches
            hit_mask = pd.Series([False] * len(seq_df), index=seq_df.index)
            for pattern in formatted_patterns:
                hit_mask = hit_mask | _fuzzy_match_series(seq_df[ColumnName.SEQUENCES], pattern, max_mismatches)
        else:
            # AND: omit only if all motifs fuzzy-match
            hit_mask = pd.Series([True] * len(seq_df), index=seq_df.index)
            for pattern in formatted_patterns:
                hit_mask = hit_mask & _fuzzy_match_series(seq_df[ColumnName.SEQUENCES], pattern, max_mismatches)
        mask = ~hit_mask
    elif partial:
        # partial=True (Yes): Omit sequences with ANY motif (OR operation)
        combined_pattern = "|".join(formatted_patterns)
        mask = ~seq_df[ColumnName.SEQUENCES].str.contains(
            combined_pattern,
            regex=True,
            case=True
        )
    else:
        # partial=False (No): Omit sequences only if they have ALL motifs (AND operation)
        hit_mask = pd.Series([True] * len(seq_df), index=seq_df.index)
        for pattern in formatted_patterns:
            pattern_mask = seq_df[ColumnName.SEQUENCES].str.contains(
                pattern,
                regex=True,
                case=True
            )
            hit_mask = hit_mask & pattern_mask
        mask = ~hit_mask
    
    # Apply filter
    filtered_df = seq_df[mask].copy()
    
    # Ensure ID column is properly formatted (following count_service pattern)
    # This creates consistent ID format: Rank=X;Reads=Y;RPU=Z for both CSV and FASTA
    filtered_df[ColumnName.ID] = (
        'Rank=' + filtered_df[ColumnName.RANK].astype(str) + ';' +
        'Reads=' + filtered_df[ColumnName.READS].astype(str) + ';' +
        'RPU=' + filtered_df[ColumnName.RPU].astype(str)
    )
    
    # Save results (no highlighting for omit functionality)
    if output_format == 'fasta':
        # Write FASTA manually to prevent BioPython's 60-character line wrapping
        # which breaks frontend parsing
        with open(output_path, 'w') as f:
            for _, row in filtered_df.iterrows():
                f.write(f">{row[ColumnName.ID]}\n")
                f.write(f"{row[ColumnName.SEQUENCES]}\n")
    else:
        # For CSV, use standard save_sequences
        file_service.save_sequences(filtered_df, output_path, output_format)
    
    return output_path
