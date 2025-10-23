from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def save_sequences(seq_df, output_path, output_format='fasta'):
    if output_format in ['fastq', 'fasta']:
        records = []
        for idx, row in seq_df.iterrows():
            seq_record = SeqRecord(
                Seq(row['Sequence']),
                id=row['ID'].split()[0],
                description=row['ID'] 
            )
            
            # Add quality scores for FASTQ output
            if output_format == 'fastq' and 'Quality' in seq_df.columns:
                # Convert Phred scores back to quality string
                seq_record.letter_annotations["phred_quality"] = row['Quality']
            
            records.append(seq_record)
        
        # Write to file
        SeqIO.write(records, output_path, output_format)
