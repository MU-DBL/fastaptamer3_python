import pandas as pd
import numpy as np
from services import file_service
from services.file_service import parse_fasta
from constants import ColumnName

def run_recount(input_path_1=None, input_path_2=None, output_path=None, output_format='fasta',scaling_factor=1e6):

    df1 = parse_fasta(input_path_1)
    df2 = parse_fasta(input_path_2)

    countMerge = pd.merge(df1, df2, on=ColumnName.SEQUENCES, how='outer', suffixes=('_x', '_y'))
    countMerge[ColumnName.READS] = countMerge[[f'{ColumnName.READS}_x', f'{ColumnName.READS}_y']].fillna(0).sum(axis=1)
    countMerge = countMerge.drop(columns=[f'{ColumnName.READS}_x', f'{ColumnName.READS}_y'])
    
    # Sort by Reads in descending order
    countMerge = countMerge.sort_values(ColumnName.READS, ascending=False).reset_index(drop=True)
    
    # Add Rank (1-indexed to match R's rowid_to_column)
    countMerge[ColumnName.RANK] = range(1, len(countMerge) + 1)
    countMerge[ColumnName.RPU] = np.round(countMerge[ColumnName.READS] / (countMerge[ColumnName.READS].sum() / scaling_factor), 0).astype(int)
    countMerge[ColumnName.LENGTH] = countMerge[ColumnName.SEQUENCES].str.len()

    countMerge[ColumnName.ID] = ColumnName.RANK + '=' + countMerge[ColumnName.RANK].astype(str) + ';' + ColumnName.READS + '=' + countMerge[ColumnName.READS].astype(str) + ';' + ColumnName.RPU + '=' + countMerge[ColumnName.RPU].astype(str)
    countMerge = countMerge[[ColumnName.ID, ColumnName.RANK, ColumnName.READS, ColumnName.RPU, ColumnName.LENGTH, ColumnName.SEQUENCES]]
    file_service.save_sequences(countMerge, output_path, output_format)
    return output_path  


