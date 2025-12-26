from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from pathlib import Path
from typing import List, Optional, Dict
import pandas as pd
import networkx as nx
from Levenshtein import distance as levenshtein_distance
import os
from services.file_service import read_file
from services.constants import ColumnName

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))

class MutationNetworkInput(BaseModel):
    input_path: str  # FASTA file path
    start_node: str
    end_node: str
    max_cost: int = 1
    output_format: str = "csv"

@router.post("/mutation-network")
async def mutation_network(params: MutationNetworkInput):
    # Validate inputs
    if not params.start_node or not params.end_node:
        raise HTTPException(status_code=400, detail="Start and end nodes required")
    
    if params.max_cost < 1:
        raise HTTPException(status_code=400, detail="Max cost must be at least 1")
    
    # Load FASTA file
    input_file = UPLOAD_DIR / params.input_path
    if not input_file.exists():
        raise HTTPException(status_code=404, detail=f"File not found: {params.input_path}")
    
    fa_df = read_file(str(input_file))
    
    # Get unique sequences and add start/end nodes
    seq_list = list(set(fa_df[ColumnName.SEQUENCES].tolist() + [params.start_node, params.end_node]))
    
    print("Starting distance calculations.")
    
    # Build connections based on edit distance
    connections = []
    n = len(seq_list)
    
    for i in range(n - 1):
        for j in range(i + 1, n):
            edit_dist = levenshtein_distance(seq_list[i], seq_list[j])
            
            if edit_dist <= params.max_cost:
                connections.append({
                    'from_vertex': seq_list[i],
                    'to_vertex': seq_list[j],
                    'cost': edit_dist
                })
    
    if not connections:
        raise HTTPException(
            status_code=404, 
            detail="No connections found within max cost threshold"
        )
    
    print("Starting path calculations.")
    
    # Create graph using NetworkX
    G = nx.Graph()
    for conn in connections:
        G.add_edge(
            conn['from_vertex'], 
            conn['to_vertex'], 
            weight=conn['cost']
        )
    
    # Find shortest path
    try:
        mutation_path = nx.shortest_path(
            G, 
            source=params.start_node, 
            target=params.end_node, 
            weight='weight'
        )
    except nx.NetworkXNoPath:
        # If no path found, return just start and end
        mutation_path = [params.start_node, params.end_node]
    except nx.NodeNotFound as e:
        raise HTTPException(
            status_code=404, 
            detail=f"Node not found in graph: {str(e)}"
        )
    
    # Check if path is trivial (only one node or direct connection failed)
    if len(mutation_path) == 1:
        raise HTTPException(
            status_code=404, 
            detail="No valid mutation path found"
        )
    
    # Build result dataframe
    result_rows = []
    for i in range(1, len(mutation_path)):
        from_seq = mutation_path[i - 1]
        to_seq = mutation_path[i]
        transition_cost = levenshtein_distance(from_seq, to_seq)
        
        result_rows.append({
            'From_Sequence': from_seq,
            'To_Sequence': to_seq,
            'Transition_Cost': transition_cost
        })
    
    result_df = pd.DataFrame(result_rows)
    
    # Save output
    base = Path(params.input_path).stem
    output = UPLOAD_DIR / f"{base}_mutation_network.{params.output_format}"
    
    if params.output_format == "csv":
        result_df.to_csv(output, index=False)
    elif params.output_format == "tsv":
        result_df.to_csv(output, sep='\t', index=False)
    else:
        raise HTTPException(status_code=400, detail=f"Unsupported format: {params.output_format}")
    
    return {
        "status": "ok",
        "result": output.name
    }