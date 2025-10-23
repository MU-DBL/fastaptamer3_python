import pandas as pd
import numpy as np
from services import file_service
from services.file_service import parse_fasta_with_reads

def run_recount(input_path_1=None, input_path_2=None, output_path=None, output_format='fasta',scaling_factor=1e6):

    df1 = parse_fasta_with_reads(input_path_1)
    df2 = parse_fasta_with_reads(input_path_2)

    countMerge = pd.merge(df1, df2, on='Sequence', how='outer', suffixes=('_x', '_y'))
    countMerge['Read'] = countMerge[['Read_x', 'Read_y']].fillna(0).sum(axis=1)
    countMerge = countMerge.drop(columns=['Read_x', 'Read_y'])
    
    # Sort by Reads in descending order
    countMerge = countMerge.sort_values('Read', ascending=False).reset_index(drop=True)
    
    # Add Rank (1-indexed to match R's rowid_to_column)
    countMerge['Rank'] = range(1, len(countMerge) + 1)
    countMerge['RPU'] = np.round(countMerge['Read'] / (countMerge['Read'].sum() / scaling_factor), 0).astype(int)
    countMerge['Length'] = countMerge['Sequence'].str.len()
    
    countMerge['ID'] = 'rank=' + countMerge['Rank'].astype(str) + ';read=' + countMerge['Read'].astype(str) + ';RPU=' + countMerge['RPU'].astype(str)
    countMerge = countMerge[['ID', 'Rank', 'Read', 'RPU', 'Length', 'Sequence']]
    file_service.save_sequences(countMerge, output_path, output_format)
    return output_path


