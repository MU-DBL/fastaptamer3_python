import pandas as pd
import numpy as np
from services import file_service
from services.file_service import parse_fasta
from services.constants import ColumnName


def run_recount(input_paths=None, output_path=None, output_format='fasta', scaling_factor=1e6):
    """
    Merge and recount sequences from 2 to 5 FASTA files.
    input_paths: list of file paths (2 to 5).

    We only care about the total Reads per sequence; Rank/RPU/ID are
    recomputed after merging, so we don't need to carry those columns
    through each merge. To avoid suffix collisions on ID/Rank/RPU, we
    restrict all intermediate dataframes to just [SEQUENCES, READS].
    """
    if not input_paths or len(input_paths) < 2:
        raise ValueError("At least 2 input files are required")
    if len(input_paths) > 5:
        raise ValueError("At most 5 input files are allowed")

    # Parse all FASTA files and keep only sequences + reads
    dfs = [parse_fasta(p)[[ColumnName.SEQUENCES, ColumnName.READS]] for p in input_paths]

    # Start from first dataframe
    countMerge = dfs[0]

    # Iteratively outer-merge on sequence and sum reads
    for i in range(1, len(dfs)):
        countMerge = pd.merge(
            countMerge,
            dfs[i],
            on=ColumnName.SEQUENCES,
            how='outer',
            suffixes=('_x', '_y'),
        )
        countMerge[ColumnName.READS] = (
            countMerge[[f"{ColumnName.READS}_x", f"{ColumnName.READS}_y"]]
            .fillna(0)
            .sum(axis=1)
        )
        countMerge = countMerge.drop(
            columns=[f"{ColumnName.READS}_x", f"{ColumnName.READS}_y"]
        )

    # Now compute rank, RPU, length, ID just like before
    countMerge = countMerge.sort_values(ColumnName.READS, ascending=False).reset_index(
        drop=True
    )
    countMerge[ColumnName.RANK] = range(1, len(countMerge) + 1)
    countMerge[ColumnName.RPU] = np.round(
        countMerge[ColumnName.READS]
        / (countMerge[ColumnName.READS].sum() / scaling_factor),
        0,
    ).astype(int)
    countMerge[ColumnName.LENGTH] = countMerge[ColumnName.SEQUENCES].str.len()
    countMerge[ColumnName.ID] = (
        ColumnName.RANK
        + "="
        + countMerge[ColumnName.RANK].astype(str)
        + ";"
        + ColumnName.READS
        + "="
        + countMerge[ColumnName.READS].astype(str)
        + ";"
        + ColumnName.RPU
        + "="
        + countMerge[ColumnName.RPU].astype(str)
    )
    countMerge = countMerge[
        [
            ColumnName.ID,
            ColumnName.RANK,
            ColumnName.READS,
            ColumnName.RPU,
            ColumnName.LENGTH,
            ColumnName.SEQUENCES,
        ]
    ]
    file_service.save_sequences(countMerge, output_path, output_format)
    return output_path


