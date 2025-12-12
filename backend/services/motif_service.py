import re
import pandas as pd
from services.constants import ColumnName
from services import file_service


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
    if partial:
        # Partial filter uses OR operation (any motif matches)
        # Combine all patterns with OR
        combined_pattern = "|".join(formatted_patterns)
        mask = seq_df[ColumnName.SEQUENCES].str.contains(combined_pattern, regex=True, case=False)
    else:
        # Full filter requires ALL motifs to be present (AND operation)
        # Each pattern must match independently
        mask = pd.Series([True] * len(seq_df), index=seq_df.index)
        for pattern in formatted_patterns:
            pattern_mask = seq_df[ColumnName.SEQUENCES].str.contains(pattern, regex=True, case=False)
            mask = mask & pattern_mask
    
    # Apply filter
    filtered_df = seq_df[mask].copy()
    
    # Highlight motifs if requested by adding parentheses around matches
    if highlight:
        def highlight_sequence(seq):
            result = seq
            for pattern in formatted_patterns:
                # Use re.sub to add parentheses around matches
                result = re.sub(f"({pattern})", r"(\1)", result, flags=re.IGNORECASE)
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
