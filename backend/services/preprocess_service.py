import asyncio
import shutil
import time
import numpy as np
import pandas as pd
from Bio import SeqIO
import regex
import os
from routers.progress import send_progress
from services import file_service
import subprocess
import state
from numba import jit

@jit(nopython=True)
def calculate_avg_error_fast(quality_scores):
    """JIT-compiled function for fast error calculation"""
    total_error = 0.0
    for q in quality_scores:
        total_error += 10 ** (-q / 10.0)
    return total_error / len(quality_scores)


def run_preprocess_slow(input_path=None, const5p="", const3p="", 
                  length_range=(10, 100), max_error=0.005,
                  output_path=None, output_format='fasta'):
    """
    Preprocess FASTA/FASTQ files with trimming and filtering
    """
    # Get file extension
    if input_path.endswith(('.fq', '.fastq', '.fq.gz', '.fastq.gz')):
        file_format = 'fastq'
    elif input_path.endswith(('.fa', '.fasta', '.fa.gz', '.fasta.gz')):
        file_format = 'fasta'
    else:
        return None
    
    # Read sequences
    records = list(SeqIO.parse(input_path, file_format))
    
    # Create DataFrame
    seq_df = pd.DataFrame({
        'ID': [rec.description for rec in records],
        'Sequence': [str(rec.seq) for rec in records]
    })
    
    # Add quality scores if FASTQ
    if file_format == 'fastq':
        seq_df['Quality'] = [rec.letter_annotations['phred_quality'] for rec in records]
    
    # Trim 5' constant region
    if const5p != "":
        seq_df = trim_constant_region(seq_df, const5p, file_format)
    
    # Trim 3' constant region
    if const3p != "":
        seq_df = trim_constant_region(seq_df, const3p, file_format)
    
    # Length filtering
    seq_df = seq_df[
        (seq_df['Sequence'].str.len() >= length_range[0]) &
        (seq_df['Sequence'].str.len() <= length_range[1])
    ]
    
    # Quality filtering
    if 'Quality' in seq_df.columns:
        seq_df['Avg_Error'] = seq_df['Quality'].apply(
            lambda q: sum([10**(-score/10) for score in q]) / len(q)
        )
        seq_df = seq_df[seq_df['Avg_Error'] <= max_error]
    
    file_service.save_sequences(seq_df, output_path, output_format)
    
    return output_path


def trim_constant_region(seq_df, const_region, file_format):
    """Trim constant region with fuzzy matching"""
    # max.distance = 0.1
    max_distance = int(len(const_region) * 0.1) + 1
    
    trimmed_sequences = []
    trimmed_qualities = []
    
    for idx, row in seq_df.iterrows():
        sequence = row['Sequence']
        
        # Fuzzy search for pattern
        match = regex.search(
            f"({const_region}){{e<={max_distance}}}", 
            sequence
        )
        
        if match:
            # Remove matched region
            start, end = match.span()
            new_seq = sequence[:start] + sequence[end:]
            trimmed_sequences.append(new_seq)
            
            # Handle quality scores
            if file_format == 'fastq':
                quality = row['Quality']
                new_qual = quality[:start] + quality[end:]
                trimmed_qualities.append(new_qual)
        else:
            trimmed_sequences.append(sequence)
            if file_format == 'fastq':
                trimmed_qualities.append(row['Quality'])
    
    seq_df['Sequence'] = trimmed_sequences
    if file_format == 'fastq':
        seq_df['Quality'] = trimmed_qualities
    
    return seq_df



async def run_preprocess(job_id, input_path, const5p="", const3p="",
                                  trim5_fixed=0, trim3_fixed=0,
                                  length_range=(10, 100), max_error=0.005,
                                  adapter_error_rate=0.1,
                                  output_path=None, output_format='fasta'):
    """
    Preprocessing with real-time progress updates
    """
    try:
        start_time = time.time()
        temp_trimmed = "/tmp/trimmed.fastq"
        
        await send_progress(job_id, 'init', 'Starting preprocessing...', 0)
        
        if input_path.endswith(('.fq', '.fastq', '.fq.gz', '.fastq.gz')):
            file_format = 'fastq'
        else:
            file_format = 'fasta'
        
        # Cutadapt command
        cmd = ['cutadapt']
        
        if trim5_fixed > 0:
            cmd.extend(['-u', str(trim5_fixed)])
        if trim3_fixed > 0:
            cmd.extend(['-u', str(-trim3_fixed)])
        if const5p:
            cmd.extend(['-g', f'{const5p}'])
        if const3p:
            cmd.extend(['-a', f'{const3p}'])
        
        cmd.extend([
            '-j', '0',  # use all available cores
            '-m', str(length_range[0]),
            '-M', str(length_range[1]),
            '-e', str(adapter_error_rate),
            '-o', temp_trimmed,
            input_path
        ])

        await send_progress(job_id, 'trimming', 'Running adapter trimming...', 10)
        cutadapt_start = time.time()
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        state.current_process = proc

        # Poll until cutadapt finishes, sending heartbeat updates so SSE stays live
        while proc.poll() is None:
            await asyncio.sleep(5)
            elapsed = int(time.time() - cutadapt_start)
            await send_progress(job_id, 'trimming', f'Adapter trimming in progress... ({elapsed}s)', min(10 + elapsed, 35))

        stdout = proc.stdout.read() if proc.stdout else ''
        stderr = proc.stderr.read() if proc.stderr else ''
        state.current_process = None
        cutadapt_time = time.time() - cutadapt_start

        if proc.returncode != 0:
            error_msg = stderr or stdout
            await send_progress(job_id, 'error', f'Cutadapt failed: {error_msg}', None)
            raise RuntimeError(f"Cutadapt failed: {error_msg}") 
        
        await send_progress(job_id, 'trimming', f'Adapter trimming completed', 40, {'time': cutadapt_time, 'output': stdout})
        
        # Fast quality filtering
        if file_format == 'fastq':
            await send_progress(job_id, 'QC', 'Starting quality filtering...', 50)
            qc_start = time.time()

            # Parse total read count from cutadapt stdout to avoid loading all records
            before_qc = 0
            for line in stdout.splitlines():
                if 'Reads written (passing filters)' in line:
                    # e.g. "Reads written (passing filters):     1,234,567 (72.1%)"
                    before_qc = int(line.split(':')[1].strip().split()[0].replace(',', ''))
                    break

            await send_progress(job_id, 'QC', f'Filtering {before_qc:,} sequences', 55)

            # Stream records one by one — no full list in memory
            after_qc = 0
            processed = 0
            batch_size = max(1, before_qc // 10) if before_qc else 10000

            with open(output_path, 'w') as out_handle:
                for rec in SeqIO.parse(temp_trimmed, 'fastq'):
                    qual_array = np.array(rec.letter_annotations['phred_quality'], dtype=np.float64)
                    avg_error = calculate_avg_error_fast(qual_array)

                    if avg_error <= max_error:
                        SeqIO.write([rec], out_handle, output_format)
                        after_qc += 1

                    processed += 1
                    if before_qc and processed % batch_size == 0:
                        progress = 55 + int(30 * processed / before_qc)
                        await send_progress(job_id, 'QC', f'Filtered {processed:,}/{before_qc:,} sequences', progress)

            qc_time = time.time() - qc_start
            pass_rate = 100 * after_qc / before_qc if before_qc > 0 else 0

            await send_progress(job_id, 'qc',
                        f'Quality filtering completed: {after_qc:,}/{before_qc:,} passed ({pass_rate:.1f}%)', 95, {'time': qc_time, 'passed': after_qc, 'total': before_qc})
            
        else:
            await send_progress(job_id, 'convert', 'Converting file format...', 85)
            conversion_start = time.time()
            if output_format != file_format:
                records = SeqIO.parse(temp_trimmed, file_format)
                SeqIO.write(records, output_path, output_format)
            else:
                shutil.move(temp_trimmed, output_path)
            conversion_time = time.time() - conversion_start
            await send_progress(job_id, 'convert', 'Conversion completed', 95, {'time': conversion_time})
        
        if os.path.exists(temp_trimmed):
            os.remove(temp_trimmed)
        
        total_time = time.time() - start_time
        await send_progress(job_id, 'complete', 'Preprocessing completed successfully!', 100, {'total_time': total_time, "output_path":os.path.basename(output_path)})
        
    except Exception as e:
        await send_progress(job_id, 'error', str(e), None)
        raise


