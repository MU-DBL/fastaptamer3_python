import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import tempfile
import os
from pathlib import Path
from services.constants import ColumnName
from pathlib import Path
from typing import List, Dict, Any, Optional
import numpy as np
import subprocess


def perform_msa(
    cluster_df: pd.DataFrame, preserve_columns: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Perform multiple sequence alignment using MUSCLE5.
    """
    # Create SeqRecord objects
    seq_records = [
        SeqRecord(Seq(row[ColumnName.SEQUENCES]), id=row[ColumnName.ID], description="")
        for _, row in cluster_df.iterrows()
    ]

    # Create temporary files
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False
    ) as input_file:
        input_path = input_file.name

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False
    ) as output_file:
        temp_output = output_file.name

    try:
        # Write sequences
        SeqIO.write(seq_records, input_path, "fasta")

        # Run MUSCLE5 using subprocess
        # MUSCLE5 uses different syntax than MUSCLE3
        process = subprocess.run(
            ["muscle", "-super5", input_path, "-output", temp_output],
            capture_output=True,
            text=True,
        )

        # Read aligned sequences
        alignment = AlignIO.read(temp_output, "fasta")
        msa_dict = {record.id: str(record.seq) for record in alignment}

    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"MUSCLE5 alignment failed: {e.stderr}")
    except FileNotFoundError:
        raise RuntimeError("MUSCLE5 not found. Check installation.")
    finally:
        for path in [input_path, temp_output]:
            if os.path.exists(path):
                os.remove(path)

    # Create MSA dataframe
    msa_df = pd.DataFrame(
        {
            ColumnName.ID: list(msa_dict.keys()),
            ColumnName.SEQUENCES: list(msa_dict.values()),
        }
    )

    # Merge with original metadata
    metadata_cols = [col for col in cluster_df.columns if col != ColumnName.SEQUENCES]
    metadata_df = cluster_df[metadata_cols]

    if preserve_columns is not None:
        metadata_cols = [ColumnName.ID] + [
            col for col in preserve_columns if col in cluster_df.columns
        ]
    else:
        metadata_cols = [
            col for col in cluster_df.columns if col != ColumnName.SEQUENCES
        ]

    metadata_df = cluster_df[metadata_cols]
    return msa_df.merge(metadata_df, on=ColumnName.ID, how="left")
