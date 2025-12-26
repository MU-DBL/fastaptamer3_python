from itertools import combinations
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from pathlib import Path
from typing import Dict, List, Set, Optional
import pandas as pd
import os
from services.file_service import read_file, save_sequences
from services.constants import ColumnName

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))

class DataMergeInput(BaseModel):
    input_paths: List[str]  # List of FASTA file paths
    merge_type: str = "outer"  # "outer", "inner", "left", "right"
    output_format: str = "csv"

class SeqPersistenceInput(BaseModel):
    merged_file_path: str  # Path to merged CSV/TSV file
    
class UpSetInput(BaseModel):
    merged_file_path: str  # Path to merged CSV/TSV file
    fasta_names: Optional[List[str]] = None  # Optional: override auto-detected names

class PersistenceResponse(BaseModel):
    status: str
    data: List[Dict[str, int]]  # [{"freq": 1, "seqCount": 150}, ...]
    total_unique_sequences: int
    max_populations: int
    populations_analyzed: int

class UpSetResponse(BaseModel):
    status: str
    sets: List[str]
    set_sizes: Dict[str, int]
    intersections: List[Dict]  # [{"sets": ["A", "B"], "size": 50, "sequences": [...]}, ...]
    total_unique_sequences: int


@router.post("/data-merge")
async def data_merge(params: DataMergeInput):
    """
    Merge multiple FASTA dataframes by Sequences column.
    
    Args:
        params: Input parameters including file paths and merge type
        
    Returns:
        Dict with status and result file path
    """
    # Validate inputs
    if not params.input_paths or len(params.input_paths) < 2:
        raise HTTPException(
            status_code=400, 
            detail="At least 2 input files required for merging"
        )
    
    if len(params.input_paths) > 26:
        raise HTTPException(
            status_code=400, 
            detail="Maximum 26 files supported (a-z suffixes)"
        )
    
    # Validate merge type
    valid_merge_types = ["outer", "inner", "left", "right"]
    if params.merge_type not in valid_merge_types:
        raise HTTPException(
            status_code=400, 
            detail=f"Invalid merge_type. Must be one of: {', '.join(valid_merge_types)}"
        )
    
    # Load all FASTA files
    fa_df_list = []
    for file_path in params.input_paths:
        input_file = UPLOAD_DIR / file_path
        if not input_file.exists():
            raise HTTPException(
                status_code=404, 
                detail=f"File not found: {file_path}"
            )
        fa_df_list.append(read_file(str(input_file)))
    
    # Rename columns for each population
    population_columns = ["ID", "Rank", "Reads", "RPU", "Cluster", "RankInCluster", "LED"]
    renamed_df_list = []
    
    for i, df in enumerate(fa_df_list):
        # Create suffix (a, b, c, ...)
        suffix = chr(ord('a') + i)
        
        # Create rename mapping for columns that exist
        rename_map = {}
        for col in population_columns:
            if col in df.columns:
                rename_map[col] = f"{col}.{suffix}"
        
        # Rename columns (keep Sequences unchanged)
        df_renamed = df.rename(columns=rename_map)
        renamed_df_list.append(df_renamed)
    
    # Perform the merge
    # Start with first dataframe
    merge_df = renamed_df_list[0]
    
    # Sequentially merge remaining dataframes
    for i in range(1, len(renamed_df_list)):
        merge_df = pd.merge(
            merge_df,
            renamed_df_list[i],
            on=ColumnName.SEQUENCES,
            how=params.merge_type
        )
    
    # Move Sequences column to first position
    cols = merge_df.columns.tolist()
    if ColumnName.SEQUENCES in cols:
        cols.remove(ColumnName.SEQUENCES)
        cols.insert(0, ColumnName.SEQUENCES)
        merge_df = merge_df[cols]
    
    # Save output
    base = Path(params.input_paths[0]).stem
    output = UPLOAD_DIR / f"{base}_merged.{params.output_format}"
    
    save_sequences(merge_df, output, params.output_format)
    
    return {
        "status": "ok",
        "result": output.name
    }

@router.post("/sequence-persistence", response_model=PersistenceResponse)
async def sequence_persistence(params: SeqPersistenceInput):
    """
    Analyze sequence persistence from merged data.
    Counts how many populations each sequence appears in.
    
    Args:
        params: Input with path to merged file
        
    Returns:
        Persistence distribution data for Angular visualization
    """
    # Load merged file
    merged_file = UPLOAD_DIR / params.merged_file_path
    if not merged_file.exists():
        raise HTTPException(
            status_code=404, 
            detail=f"Merged file not found: {params.merged_file_path}"
        )
    
    # Read merged data
    merge_df = read_file(str(merged_file))
    
    if ColumnName.SEQUENCES not in merge_df.columns:
        raise HTTPException(
            status_code=400, 
            detail="Merged file must contain 'Sequences' column"
        )
    
    # Detect populations from column suffixes (e.g., ID.a, ID.b, Reads.a, Reads.b)
    # Get unique suffixes from column names
    suffixes = set()
    for col in merge_df.columns:
        if '.' in col and col != ColumnName.SEQUENCES:
            suffix = col.split('.')[-1]
            suffixes.add(suffix)
    
    if not suffixes:
        raise HTTPException(
            status_code=400,
            detail="No population suffixes found in merged data. Expected columns like 'ID.a', 'Reads.b', etc."
        )
    
    populations = sorted(list(suffixes))
    
    # For each sequence, count how many populations it appears in
    # A sequence "appears" in a population if ANY column for that population is not null
    persistence_counts = []
    
    for idx, row in merge_df.iterrows():
        sequence = row[ColumnName.SEQUENCES]
        
        # Count how many populations have data for this sequence
        pop_count = 0
        for suffix in populations:
            # Check if any column with this suffix has a non-null value
            suffix_cols = [col for col in merge_df.columns if col.endswith(f'.{suffix}')]
            if any(pd.notna(row[col]) for col in suffix_cols):
                pop_count += 1
        
        persistence_counts.append(pop_count)
    
    # Count how many sequences appear in 1 pop, 2 pops, 3 pops, etc.
    persistence_series = pd.Series(persistence_counts).value_counts().sort_index()
    
    # Format for Angular
    data = [
        {"freq": int(freq), "seqCount": int(count)}
        for freq, count in persistence_series.items()
    ]
    
    return PersistenceResponse(
        status="ok",
        data=data
    )

@router.post("/upset-data", response_model=UpSetResponse)
async def upset_data(params: UpSetInput):
    """
    Compute UpSet intersection data from merged file.
    Returns JSON data for Angular to visualize.
    
    Args:
        params: Input with path to merged file and optional set names
        
    Returns:
        UpSet intersection data for Angular visualization
    """
    # Load merged file
    merged_file = UPLOAD_DIR / params.merged_file_path
    if not merged_file.exists():
        raise HTTPException(
            status_code=404,
            detail=f"Merged file not found: {params.merged_file_path}"
        )
    
    # Read merged data
    merge_df = read_file(str(merged_file))
    
    if ColumnName.SEQUENCES not in merge_df.columns:
        raise HTTPException(
            status_code=400,
            detail="Merged file must contain 'Sequences' column"
        )
    
    # Detect populations from column suffixes
    suffixes = set()
    for col in merge_df.columns:
        if '.' in col and col != ColumnName.SEQUENCES:
            suffix = col.split('.')[-1]
            suffixes.add(suffix)
    
    if not suffixes:
        raise HTTPException(
            status_code=400,
            detail="No population suffixes found in merged data"
        )
    
    populations = sorted(list(suffixes))
    
    # Use provided names or default to suffix names
    if params.fasta_names:
        if len(params.fasta_names) != len(populations):
            raise HTTPException(
                status_code=400,
                detail=f"Number of fasta_names ({len(params.fasta_names)}) must match detected populations ({len(populations)})"
            )
        set_names = params.fasta_names
    else:
        set_names = [f"Population_{suffix}" for suffix in populations]
    
    # Build dictionary of sequences for each population
    sequence_dict: Dict[str, Set[str]] = {}
    
    for suffix, name in zip(populations, set_names):
        # Get all sequences where this population has data
        suffix_cols = [col for col in merge_df.columns if col.endswith(f'.{suffix}')]
        
        # Filter rows where at least one column for this population is not null
        mask = merge_df[suffix_cols].notna().any(axis=1)
        sequences = set(merge_df.loc[mask, ColumnName.SEQUENCES].tolist())
        
        sequence_dict[name] = sequences
    
    # Calculate set sizes
    set_sizes = {name: len(seqs) for name, seqs in sequence_dict.items()}
    
    # Calculate all intersections
    intersections = []
    
    # For each possible combination of sets
    for r in range(1, len(sequence_dict) + 1):
        for combo in combinations(sequence_dict.keys(), r):
            combo_list = list(combo)
            
            # Get intersection of selected sets
            intersection_seqs = set.intersection(*[sequence_dict[name] for name in combo_list])
            
            # Exclude sequences that are in OTHER sets (strict intersection)
            for other_name in sequence_dict.keys():
                if other_name not in combo_list:
                    intersection_seqs -= sequence_dict[other_name]
            
            if intersection_seqs:
                intersections.append({
                    "sets": combo_list,
                    "size": len(intersection_seqs),
                    "sequences": sorted(list(intersection_seqs))[:100]  # Limit to first 100
                })
    
    # Sort by size descending
    intersections.sort(key=lambda x: x["size"], reverse=True)
    
    # Calculate total unique sequences
    all_sequences = set()
    for seqs in sequence_dict.values():
        all_sequences.update(seqs)
    
    return UpSetResponse(
        status="ok",
        sets=set_names,
        set_sizes=set_sizes,
        intersections=intersections,
        total_unique_sequences=len(all_sequences)
    )
