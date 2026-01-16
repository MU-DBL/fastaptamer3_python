from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from pathlib import Path
from typing import Optional, Dict, List
import pandas as pd
import numpy as np
import os
import re
import json
from services.file_service import read_file, save_sequences
from services.constants import ColumnName

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))
TRANSLATIONS_FILE = Path(__file__).parent.parent / "data" / "translations.json"

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
        
        # Step 3: Translate sequences
        translated_sequences = []
        
        for idx, row in fa_df.iterrows():
            sequence = str(row[ColumnName.SEQUENCES])
            translated_seq = translate_sequence(
                sequence,
                translation_map,
                params.orf
            )
            translated_sequences.append(translated_seq)
        
        # Create result DataFrame
        translate_df = fa_df.copy()
        translate_df[ColumnName.SEQUENCES] = translated_sequences
        
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
                'RPU=' + converged_df[ColumnName.RPU].astype(str)
            )
            
            translate_df = converged_df[[
                ColumnName.ID,
                ColumnName.RANK,
                ColumnName.READS,
                ColumnName.RPU,
                'Unique_Nts',
                ColumnName.SEQUENCES
            ]].copy()
        else:
            # When not converging, update ID to maintain consistent format
            translate_df = translate_df.copy()
            translate_df[ColumnName.ID] = (
                'Rank=' + translate_df[ColumnName.RANK].astype(str) + ';' +
                'Reads=' + translate_df[ColumnName.READS].astype(str) + ';' +
                'RPU=' + translate_df[ColumnName.RPU].astype(str)
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
            detail=f"Translation failed: {str(e)}"
        )
        
        
        