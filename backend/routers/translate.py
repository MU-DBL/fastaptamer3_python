from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from pathlib import Path
from typing import Optional, Dict, List
import pandas as pd
import numpy as np
import os
import re
import json
from services.file_service import read_file, save_sequences, scaling_factor_of, format_rpu, RPU_DECIMALS
from services.constants import ColumnName

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))
TRANSLATIONS_FILE = Path(__file__).parent.parent / "data" / "translations.json"

# Nucleotide codes for vectorized codon lookup (codon index = base-4 number); other characters are -1
_NT_CODE = np.full(256, -1, dtype=np.int64)
for _code, _nt in enumerate('ACGT'):
    _NT_CODE[ord(_nt)] = _code
_ACGT_CODONS = [a + b + c for a in 'ACGT' for b in 'ACGT' for c in 'ACGT']

class CodonChange(BaseModel):
    Codon: str
    Translation: str

class TranslateInput(BaseModel):
    input_path: str
    orf: int = 1  # Open reading frame
    converge: bool = True  # Whether to merge non-unique amino acid sequences
    input_changes: Optional[List[CodonChange]] = None  # Custom codon modifications
    translate_selection: str = "Standard"  # Genetic code selection
    output_format: str = "csv"

def load_translation_mapping(translate_selection: str = "Standard") -> pd.DataFrame:

    try:
        with open(TRANSLATIONS_FILE, 'r') as f:
            translations_data = json.load(f)
        
        # Filter for selected genetic code
        translation_df = pd.DataFrame(translations_data)
        translation_df = translation_df[
            translation_df['whichGroup'] == translate_selection
        ][['Codon', 'Translation']].copy()
        
        if translation_df.empty:
            raise ValueError(f"Genetic code '{translate_selection}' not found")
        
        return translation_df
    
    except FileNotFoundError:
        # Return standard genetic code as fallback
        return get_standard_genetic_code()

def get_standard_genetic_code() -> pd.DataFrame:
    """
    Return standard genetic code mapping.
    
    Returns:
        DataFrame with standard codon to amino acid mappings
    """
    standard_code = {
        # Standard genetic code
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    
    return pd.DataFrame(list(standard_code.items()), columns=['Codon', 'Translation'])

def apply_translation_changes(
    translation_df: pd.DataFrame,
    input_changes: Optional[List[CodonChange]]
) -> pd.DataFrame:
    """
    Apply custom modifications to genetic code translation mapping.
    
    Args:
        translation_df: Base translation mapping
        input_changes: List of custom codon changes
    
    Returns:
        Modified translation DataFrame
    """
    if not input_changes:
        return translation_df
    
    # Convert input_changes to DataFrame
    changes_data = []
    for change in input_changes:
        # Remove special characters
        codon = re.sub(r'[^a-zA-Z0-9]', '', change.Codon).upper()
        translation = re.sub(r'[^a-zA-Z0-9*]', '', change.Translation).upper()
        
        # Validate: codons must be 3 characters, translations must be 1 character
        if len(codon) == 3 and len(translation) == 1:
            changes_data.append({'Codon': codon, 'Translation': translation})
    
    if not changes_data:
        return translation_df
    
    changes_df = pd.DataFrame(changes_data)
    
    # Merge changes with original translation (changes override originals)
    # R equivalent: dplyr::full_join(...) with conditional mutation
    merged_df = translation_df.merge(
        changes_df,
        on='Codon',
        how='outer',
        suffixes=('_original', '_new')
    )
    
    # Use new translation if available, otherwise use original
    merged_df['Translation'] = merged_df['Translation_new'].fillna(
        merged_df['Translation_original']
    )
    
    result_df = merged_df[['Codon', 'Translation']].copy()
    
    return result_df

def translate_sequence(sequence: str, translation_map: Dict[str, str], orf: int = 1) -> str:
    # Convert to uppercase and apply ORF
    sequence = sequence.upper()
    sequence = sequence[orf - 1:]  # orf is 1-based in R
    
    # Remove trailing bases if length not divisible by 3
    remainder = len(sequence) % 3
    if remainder > 0:
        sequence = sequence[:-remainder]
    
    # Translate codon by codon
    amino_acids = []
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        # Use translation map, default to 'X' for ambiguous codons
        amino_acid = translation_map.get(codon, 'X')
        amino_acids.append(amino_acid)
    
    return ''.join(amino_acids)

def translate_sequence_list(sequences: List[str], translation_map: Dict[str, str], orf: int = 1) -> List[str]:
    """Translate many sequences at once; same rules and output as translate_sequence on each."""
    # Uppercase, apply ORF, drop trailing bases that do not complete a codon
    seqs = [s.upper()[orf - 1:] for s in sequences]
    seqs = [s[:len(s) - len(s) % 3] for s in seqs]
    n_codons = np.fromiter((len(s) // 3 for s in seqs), dtype=np.int64, count=len(seqs))
    if n_codons.sum() == 0:
        return [''] * len(seqs)

    # Every codon of every sequence as a row of nucleotide codes
    joined = ''.join(seqs)
    codes = _NT_CODE[np.frombuffer(joined.encode('ascii', errors='replace'), dtype=np.uint8)].reshape(-1, 3)
    acgt = (codes >= 0).all(axis=1)

    # ACGT codons: look up through a 64-entry table; any other codon: the map itself (default 'X')
    table = np.array([translation_map.get(c, 'X') for c in _ACGT_CODONS], dtype=object)
    amino_acids = np.empty(len(codes), dtype=object)
    amino_acids[acgt] = table[codes[acgt] @ np.array([16, 4, 1])]
    for i in np.flatnonzero(~acgt):
        amino_acids[i] = translation_map.get(joined[3 * i:3 * i + 3], 'X')

    # Split the codon translations back into one protein per sequence
    ends = np.cumsum(n_codons)
    starts = ends - n_codons
    return [''.join(amino_acids[a:b]) for a, b in zip(starts, ends)]

@router.post("/translate")
async def translate_sequences(params: TranslateInput):
    if not params.input_path:
        raise HTTPException(status_code=400, detail="input_path is required")
    
    if params.orf < 1:
        raise HTTPException(status_code=400, detail="orf must be >= 1")
    
    filepath = UPLOAD_DIR / params.input_path
    if not filepath.exists():
        raise HTTPException(
            status_code=404,
            detail=f"File not found: {params.input_path}"
        )
    
    try:
        # Step 1: Read input file
        fa_df = read_file(filepath)
        
        if fa_df.empty:
            raise HTTPException(status_code=400, detail="Input file is empty")
        
        # Step 2: Load and modify translation mapping
        translation_df = load_translation_mapping(params.translate_selection)
        
        # Apply custom changes if provided
        translation_df = apply_translation_changes(
            translation_df,
            params.input_changes
        )
        
        # Convert to dictionary for faster lookup
        translation_map = dict(zip(
            translation_df['Codon'],
            translation_df['Translation']
        ))
        
        # Step 3: Translate sequences (all at once; a per-row loop is slow on large files)
        translated_sequences = translate_sequence_list(
            fa_df[ColumnName.SEQUENCES].astype(str).tolist(),
            translation_map,
            params.orf
        )
        
        # Create result DataFrame
        translate_df = fa_df.copy()
        translate_df[ColumnName.SEQUENCES] = translated_sequences

        # U of the input's RPU, carried into the output IDs (None for inputs that do not record it)
        u = scaling_factor_of(translate_df[ColumnName.U]) if ColumnName.U in translate_df.columns else None
        u_suffix = f';{ColumnName.U}={u}' if u is not None else ''

        # Step 4: Optionally converge (merge non-unique amino acid sequences)
        if params.converge:
            # Group by sequences and sum reads/RPU
            converged_df = (
                translate_df
                .groupby(ColumnName.SEQUENCES, as_index=False)
                .agg({
                    ColumnName.READS: 'sum',
                    ColumnName.RPU: 'sum'
                })
            )
            
            # Sort by reads (descending) and assign new rank
            converged_df = converged_df.sort_values(
                ColumnName.READS,
                ascending=False
            ).reset_index(drop=True)
            # Summed RPU picks up float noise (e.g. 0.30000000000000004); keep the usual decimals
            converged_df[ColumnName.RPU] = converged_df[ColumnName.RPU].round(RPU_DECIMALS)

            converged_df[ColumnName.RANK] = range(1, len(converged_df) + 1)
            
            # Count unique nucleotide sequences for each amino acid sequence
            unique_nt_counts = translate_df.groupby(
                ColumnName.SEQUENCES
            ).size().reset_index(name='Unique_Nts')
            
            converged_df = converged_df.merge(
                unique_nt_counts,
                on=ColumnName.SEQUENCES,
                how='left'
            )
            
            # Recreate ID column in semicolon format (consistent with count_service pattern)
            converged_df[ColumnName.ID] = (
                'Rank=' + converged_df[ColumnName.RANK].astype(str) + ';' +
                'Reads=' + converged_df[ColumnName.READS].astype(str) + ';' +
                'RPU=' + format_rpu(converged_df[ColumnName.RPU]) + u_suffix
            )

            output_cols = [ColumnName.ID, ColumnName.RANK, ColumnName.READS, ColumnName.RPU]
            if u is not None:
                converged_df[ColumnName.U] = u
                output_cols.append(ColumnName.U)
            translate_df = converged_df[output_cols + ['Unique_Nts', ColumnName.SEQUENCES]].copy()
        else:
            # When not converging, update ID to maintain consistent format
            translate_df = translate_df.copy()
            translate_df[ColumnName.ID] = (
                'Rank=' + translate_df[ColumnName.RANK].astype(str) + ';' +
                'Reads=' + translate_df[ColumnName.READS].astype(str) + ';' +
                'RPU=' + format_rpu(translate_df[ColumnName.RPU]) + u_suffix
            )
        
        # Step 5: Add sequence length column
        translate_df['Length'] = translate_df[ColumnName.SEQUENCES].str.len()
        
        # Reorder columns to put Length before Sequences
        cols = translate_df.columns.tolist()
        seq_idx = cols.index(ColumnName.SEQUENCES)
        cols.insert(seq_idx, cols.pop(cols.index('Length')))
        translate_df = translate_df[cols]
        
        # Step 6: Save results
        base_name = os.path.splitext(os.path.basename(params.input_path))[0]
        output_format = params.output_format.lower()
        output_path = UPLOAD_DIR / f"{base_name}_translated.{output_format}"
        save_sequences(translate_df, str(output_path), params.output_format)
        
        return {
            "status": "ok",
            "result": os.path.basename(output_path)
        }
    
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Translation failed ({type(e).__name__}): {str(e)}"
        )
        
        
        