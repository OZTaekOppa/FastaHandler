# This script to find and trim anchor sequecnes from input fasta.
# Usage: python3 find_anchro_trim.py --input-fa ref_seq.fa --anchor1-fa anchor1_seq.fa --anchor2-fa anchor2_seq.fa --out ./ --t 1 --mem 10

#!/usr/bin/env python3
import argparse
import os
import sys
import gzip
import bz2
import zipfile
import logging
import multiprocessing as mp
import parasail
from Bio import SeqIO

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Script to find anchors in sequences and trim input sequence accordingly.')
    parser.add_argument('--input-fa', required=True,
                        help='Indicate the input file and path of fasta file')
    parser.add_argument('--anchor1-fa', required=True,
                        help='Indicate the first anchor file and path of fasta file')
    parser.add_argument('--anchor2-fa', required=True,
                        help='Indicate the second anchor file and path of fasta file')
    parser.add_argument('--out', required=True,
                        help='Indicate the output directory path')
    parser.add_argument('--t', type=int, default=1,
                        help='Number of CPUs with only numbers (integers)')
    parser.add_argument('--mem', type=int, default=10,
                        help='Number of memory in GB with only numbers (integers)')
    args = parser.parse_args()

    if args.t < 1 or not isinstance(args.t, int):
        sys.exit("Not an integer number. Please provide only readable integer numbers to call the proper CPUs and Memory.")

    if args.mem < 1 or not isinstance(args.mem, int):
        sys.exit("Not an integer number. Please provide only readable integer numbers to call the proper CPUs and Memory.")

    return args

def setup_logging(output_dir):
    log_file = os.path.join(output_dir, 'find_anchor_trim.log')
    logging.basicConfig(filename=log_file,
                        level=logging.INFO,
                        format='%(asctime)s %(levelname)s:%(message)s')

def decompress_file(filepath):
    if filepath.endswith('.gz'):
        decompressed_path = filepath[:-3]
        with gzip.open(filepath, 'rt', encoding='utf-8') as f_in, open(decompressed_path, 'w', encoding='utf-8') as f_out:
            f_out.write(f_in.read())
        return decompressed_path
    elif filepath.endswith('.bz2'):
        decompressed_path = filepath[:-4]
        with bz2.open(filepath, 'rt', encoding='utf-8') as f_in, open(decompressed_path, 'w', encoding='utf-8') as f_out:
            f_out.write(f_in.read())
        return decompressed_path
    elif filepath.endswith('.zip'):
        with zipfile.ZipFile(filepath, 'r') as zip_ref:
            zip_ref.extractall(os.path.dirname(filepath))
            extracted_files = zip_ref.namelist()
            if len(extracted_files) == 1:
                return os.path.join(os.path.dirname(filepath), extracted_files[0])
            else:
                sys.exit("Zip file contains multiple files. Please provide a zip file with a single file.")
    else:
        return filepath  # No decompression needed

def is_valid_file(filename):
    extensions = ['.fasta', '.fa', '.fna', '.gz', '.zip', '.bz2']
    if any(filename.endswith(ext) for ext in extensions):
        return True
    else:
        return False

def read_fasta_file(filepath):
    try:
        sequences = []
        for record in SeqIO.parse(filepath, 'fasta'):
            seq = str(record.seq)
            sequences.append(seq)
        if len(sequences) == 0:
            sys.exit(f"No sequences found in {filepath}")
        # Concatenate sequences into a single sequence
        full_sequence = ''.join(sequences).upper()
        return full_sequence
    except Exception as e:
        sys.exit(f"Failed to read fasta file {filepath}: {str(e)}")

def write_fasta_file(filepath, header, sequence):
    with open(filepath, 'w') as f:
        f.write(f">{header}\n")
        f.write(f"{sequence}\n")

def smith_waterman_alignment(input_seq, anchor_seq):
    # Perform Smith-Waterman local alignment using parasail
    # Use the appropriate substitution matrix and gap penalties
    matrix = parasail.matrix_create("ACGTN", 2, -1)
    # Use the function that returns alignment statistics
    result = parasail.sw_stats(anchor_seq, input_seq, 10, 1, matrix)
    return result

def find_alignments(args):
    input_anchor, input_seq, anchor_seq = args
    # Perform alignment
    result = smith_waterman_alignment(input_seq, anchor_seq)
    # Extract alignment information
    matches = result.matches
    length = result.length
    mismatches = length - matches
    seq_identity = (matches / length) * 100 if length > 0 else 0
    idx_start = result.end_ref - length + 1  # Convert to 1-based index
    idx_end = result.end_ref + 1  # Inclusive
    return {
        'Input_anchor': input_anchor,
        'Idx_start': idx_start,
        'Idx_end': idx_end,
        'Mism_num': mismatches,
        'Seq_idty': seq_identity
    }

def main():
    args = parse_arguments()
    output_dir = args.out
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    setup_logging(output_dir)
    logging.info("Starting processing.")

    # Handle input-fa
    if not is_valid_file(args.input_fa):
        sys.exit('Not a proper file format. Please provide a readable file format to call the "find_anchor_trim module."')
    input_fa_file = decompress_file(args.input_fa)

    # Handle anchor1-fa
    if not is_valid_file(args.anchor1_fa):
        sys.exit('Not a proper file format. Please provide a readable file format to call the "find_anchor_trim module."')
    anchor1_fa_file = decompress_file(args.anchor1_fa)

    # Handle anchor2-fa
    if not is_valid_file(args.anchor2_fa):
        sys.exit('Not a proper file format. Please provide a readable file format to call the "find_anchor_trim module."')
    anchor2_fa_file = decompress_file(args.anchor2_fa)

    # Read sequences
    input_seq = read_fasta_file(input_fa_file)
    anchor1_seq = read_fasta_file(anchor1_fa_file)
    anchor2_seq = read_fasta_file(anchor2_fa_file)

    # Write input_index.txt
    input_index_path = os.path.join(output_dir, 'input_index.txt')
    input_length = len(input_seq)
    with open(input_index_path, 'w') as f:
        f.write(f"Input Index: 1-{input_length}\n")
    logging.info(f"input_index.txt saved to {input_index_path}")

    # Prepare alignment tasks
    tasks = [
        ('Anchor1', input_seq, anchor1_seq),
        ('Anchor2', input_seq, anchor2_seq)
    ]

    # Perform alignments using multiprocessing
    pool = mp.Pool(args.t)
    results = pool.map(find_alignments, tasks)
    pool.close()
    pool.join()

    # Filter results based on mismatch number and sequence identity
    filtered_results = []
    for res in results:
        mism_num = res['Mism_num']
        seq_idty = res['Seq_idty']
        if mism_num <= 50 and seq_idty >= 90:
            res['Seq_idty'] = f"{seq_idty:.0f}%"
            filtered_results.append(res)

    # Write anchor_match_index.txt
    anchor_match_index_path = os.path.join(output_dir, 'anchor_match_index.txt')
    with open(anchor_match_index_path, 'w') as f:
        f.write("Input_anchor\tIdx_start\tIdx_end\tMism_num\tSeq_idty\n")
        for res in filtered_results:
            f.write(f"{res['Input_anchor']}\t{res['Idx_start']}\t{res['Idx_end']}\t{res['Mism_num']}\t{res['Seq_idty']}\n")
    logging.info(f"anchor_match_index.txt saved to {anchor_match_index_path}")

    # Find best matches for trimming
    best_matches = {}
    for res in filtered_results:
        input_anchor = res['Input_anchor']
        if input_anchor not in best_matches:
            best_matches[input_anchor] = res
        else:
            current_best = best_matches[input_anchor]
            if res['Mism_num'] < current_best['Mism_num']:
                best_matches[input_anchor] = res
            elif res['Mism_num'] == current_best['Mism_num'] and float(res['Seq_idty'].strip('%')) > float(current_best['Seq_idty'].strip('%')):
                best_matches[input_anchor] = res

    # Trim input sequence based on best matches
    # For Anchor1, trim from start to idx_end
    # For Anchor2, trim from idx_start to end
    trim_start = 0
    trim_end = len(input_seq)
    if 'Anchor1' in best_matches:
        trim_start = int(best_matches['Anchor1']['Idx_end'])
    if 'Anchor2' in best_matches:
        trim_end = int(best_matches['Anchor2']['Idx_start']) - 1  # -1 because idx_start is inclusive
    trimmed_seq = input_seq[trim_start:trim_end]

    # Write trimmed_input.fa
    trimmed_input_path = os.path.join(output_dir, 'trimmed_input.fa')
    # Get header from input file
    for record in SeqIO.parse(input_fa_file, 'fasta'):
        header = record.id
        break
    write_fasta_file(trimmed_input_path, header, trimmed_seq)
    logging.info(f"trimmed_input.fa saved to {trimmed_input_path}")

    logging.info("Processing completed successfully.")

if __name__ == "__main__":
    main()
