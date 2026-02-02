from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from pathlib import Path
from typing import List
import pandas as pd
import networkx as nx
from Levenshtein import distance as levenshtein_distance
import os

from services.file_service import read_file
from services.constants import ColumnName

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))


class MutationNetworkInput(BaseModel):
    input_path: str
    start_node: str
    end_node: str
    max_cost: int = 1
    output_format: str = "csv"


def norm(seq: str) -> str:
    """Normalize DNA / protein sequences"""
    return seq.strip().upper()


@router.post("/mutation-network")
async def mutation_network(params: MutationNetworkInput):

    if params.max_cost < 1:
        raise HTTPException(status_code=400, detail="Max cost must be at least 1")

    input_file = UPLOAD_DIR / params.input_path
    if not input_file.exists():
        raise HTTPException(
            status_code=404,
            detail=f"File not found: {params.input_path}"
        )

    # Load FASTA / CSV
    fa_df = read_file(input_file)

    if ColumnName.SEQUENCES not in fa_df.columns:
        raise HTTPException(
            status_code=400,
            detail=f"Missing column: {ColumnName.SEQUENCES}"
        )

    # Normalize sequences
    fa_df[ColumnName.SEQUENCES] = fa_df[ColumnName.SEQUENCES].map(norm)
    start_node = norm(params.start_node)
    end_node = norm(params.end_node)

    # Build unique sequence list (ensure start/end included)
    seq_list = list(
        set(fa_df[ColumnName.SEQUENCES].tolist() + [start_node, end_node])
    )

    # Build connections
    connections = []
    n = len(seq_list)

    for i in range(n - 1):
        for j in range(i + 1, n):
            d = levenshtein_distance(seq_list[i], seq_list[j])
            if d <= params.max_cost:
                connections.append(
                    (seq_list[i], seq_list[j], d)
                )

    if not connections:
        raise HTTPException(
            status_code=404,
            detail="No connections found within max cost threshold"
        )

    # Create graph
    G = nx.Graph()

    G.add_nodes_from(seq_list)

    for u, v, w in connections:
        G.add_edge(u, v, weight=w)

    # Check connectivity before shortest path
    if start_node not in G or end_node not in G:
        raise HTTPException(
            status_code=404,
            detail="Start or end node has no valid connections"
        )

    try:
        mutation_path = nx.shortest_path(
            G,
            source=start_node,
            target=end_node,
            weight="weight"
        )
    except nx.NetworkXNoPath:
        raise HTTPException(
            status_code=404,
            detail="No mutation path found between start and end nodes"
        )

    # Build output
    rows = []
    for i in range(1, len(mutation_path)):
        a = mutation_path[i - 1]
        b = mutation_path[i]
        rows.append({
            "From_Sequence": a,
            "To_Sequence": b,
            "Transition_Cost": levenshtein_distance(a, b)
        })

    result_df = pd.DataFrame(rows)

    base = Path(params.input_path).stem
    output = UPLOAD_DIR / f"{base}_mutation_network.{params.output_format}"

    if params.output_format == "csv":
        result_df.to_csv(output, index=False)
    elif params.output_format == "tsv":
        result_df.to_csv(output, sep="\t", index=False)
    else:
        raise HTTPException(
            status_code=400,
            detail=f"Unsupported format: {params.output_format}"
        )

    return {
        "status": "ok",
        "result": output.name,
        "path_length": len(mutation_path)
    }
