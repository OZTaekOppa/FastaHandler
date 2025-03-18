# This script to find unique and similariy sequecnes against its own input sequecne
# Created by Hyungtaek Jung
# Usage: python3 ./size_uniq_ptrn_search.py --input-fa input.fa.gz --target-fa target.fa.gz --char-size 100 --out output.txt (default both direction search)

#!/usr/bin/env python3
import argparse
import os
import sys
import gzip
import bz2
import zipfile
import multiprocessing as mp
import logging
from functools import partial
from Bio.Align import PairwiseAligner

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Script to find mismatched numbers and sequence identity in fasta sequences.')
    parser.add_argument('--input-fa', required=True,
                        help='Indicate the input file and path of fasta file')
    parser.add_argument('--target-fa', required=True,
                        help='Indicate the target file and path of fasta file')
    parser.add_argument('--char-size', required=True, type=int,
                        help='Indicate the character or sequence size to search pattern')
    parser.add_argument('--both-srch', action='store_true',
                        help='Indicate the --char-size direction to both ways (ascending and descending positions)')
    parser.add_argument('--out', required=True,
                        help='Indicate the out file and path')
    parser.add_argument('--t', type=int, default=1,
                        help='Number of CPUs with only numbers (integers)')
    parser.add_argument('--mem', type=int, default=10,
                        help='Number of memory with only numbers (integers)')
    args = parser.parse_args()

    # Default to both search if no direction is specified
    if not args.both_srch:
        args.both_srch = True

    if args.t < 1 or not isinstance(args.t, int):
        sys.exit("Not an integer number. Please provide only readable integer numbers to call the proper CPUs and Memory.")

    if args.mem < 1 or not isinstance(args.mem, int):
        sys.exit("Not an integer number. Please provide only readable integer numbers to call the proper CPUs and Memory.")

    return args

def setup_logging(output_dir):
    log_file = os.path.join(output_dir, 'size_uniq_ptrn_search.log')
    logging.basicConfig(filename=log_file,
                        level=logging.INFO,
                        format='%(asctime)s %(levelname)s:%(message)s')

def is_fasta_file(filename):
    extensions = ['.fa', '.fasta', '.fna']
    compressed_extensions = ['.fa.gz', '.fasta.gz', '.fna.gz',
                             '.fa.bz2', '.fasta.bz2', '.fna.bz2',
                             '.fa.zip', '.fasta.zip', '.fna.zip']
    if any(filename.endswith(ext) for ext in extensions + compressed_extensions):
        return True
    else:
        return False

def decompress_file(filepath):
    if filepath.endswith('.gz'):
        decompressed_path = filepath[:-3]
        with gzip.open(filepath, 'rt') as f_in, open(decompressed_path, 'w') as f_out:
            f_out.write(f_in.read())
        return decompressed_path
    elif filepath.endswith('.bz2'):
        decompressed_path = filepath[:-4]
        with bz2.open(filepath, 'rt') as f_in, open(decompressed_path, 'w') as f_out:
            f_out.write(f_in.read())
        return decompressed_path
    elif filepath.endswith('.zip'):
        with zipfile.ZipFile(filepath, 'r') as zip_ref:
            zip_ref.extractall(os.path.dirname(filepath))
            return os.path.join(os.path.dirname(filepath), zip_ref.namelist()[0])
    else:
        return filepath

def read_fasta(filepath):
    if not os.path.exists(filepath):
        sys.exit(f"File not found: {filepath}")

    seq = ''
    with open(filepath, 'r') as file:
        for line in file:
            if line.startswith('>'):
                continue  # Skip header lines
            else:
                seq += line.strip()
    return seq.lower()

def find_match(input_seq, target_seq):
    idx = input_seq.find(target_seq)
    if idx != -1:
        start_idx = idx
        end_idx = idx + len(target_seq) - 1
        mismatches = sum(1 for a, b in zip(input_seq[start_idx:end_idx+1], target_seq) if a != b)
    else:
        # If not found, set indices to -1
        start_idx = -1
        end_idx = -1
        mismatches = len(target_seq)
    return start_idx, end_idx, mismatches

def write_input_index(input_seq_len, output_dir):
    idx_info_path = os.path.join(output_dir, 'Input_Index.txt')
    with open(idx_info_path, 'w') as f:
        # Adjust indices to start from 1
        f.write(f"Input Index: 1-{input_seq_len}\n")
    logging.info(f"Input Index information saved to {idx_info_path}")

def write_match_index(match_start, match_end, mismatches, output_dir):
    match_idx_path = os.path.join(output_dir, 'Input_Match_Index.txt')
    with open(match_idx_path, 'w') as f:
        if match_start != -1:
            # Adjust indices to start from 1
            f.write(f"Input Match Index: {match_start+1}-{match_end+1}\n")
        else:
            f.write("Input Match Index: Not Found\n")
        f.write(f"Input Mismatch Number: {mismatches}\n")
    logging.info(f"Input Match Index information saved to {match_idx_path}")

def generate_chunks(input_seq_len, match_start, match_end, char_size):
    chunks = []
    # Forward search
    idx = match_end + 1 if match_end != -1 else 0
    while idx < input_seq_len:
        idx_end = min(idx + char_size - 1, input_seq_len - 1)
        if idx_end >= input_seq_len:
            break
        chunks.append((idx, idx_end))
        idx += char_size
    # Backward search
    idx = match_start - char_size if match_start != -1 else input_seq_len - char_size
    while idx >= 0:
        idx_end = idx + char_size - 1
        if idx_end >= input_seq_len:
            idx -= char_size
            continue
        chunks.append((idx, idx_end))
        idx -= char_size
    # Exclude matched region
    if match_start != -1:
        chunks = [chunk for chunk in chunks if not (chunk[0] >= match_start and chunk[1] <= match_end)]
    return sorted(chunks, key=lambda x: x[0])

def compute_min_mismatches_and_identity(chunk_info):
    chunk_seq, input_seq, chunk_pos = chunk_info
    min_mismatches = len(chunk_seq)
    max_identity = 0.0
    # Exclude the chunk's own position
    excluded_start = chunk_pos[0]
    excluded_end = chunk_pos[1]
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = 0
    aligner.extend_gap_score = 0
    for i in range(len(input_seq) - len(chunk_seq) + 1):
        if i >= excluded_start and i <= excluded_end:
            continue  # Skip the chunk's own position
        window_seq = input_seq[i:i+len(chunk_seq)]
        # Perform Smith-Waterman alignment
        alignments = aligner.align(chunk_seq, window_seq)
        if alignments:
            # Take the alignment with the highest score
            best_alignment = alignments[0]
            # Calculate sequence identity
            aligned_chunk_seq = best_alignment.aligned[0]
            aligned_window_seq = best_alignment.aligned[1]
            matches = sum(1 for pos1, pos2 in zip(aligned_chunk_seq, aligned_window_seq) if chunk_seq[pos1[0]:pos1[1]] == input_seq[pos2[0]:pos2[1]])
            total = sum(pos1[1] - pos1[0] for pos1 in aligned_chunk_seq)
            identity = (matches / total) * 100 if total > 0 else 0.0
            mismatches = total - matches
            if mismatches < min_mismatches:
                min_mismatches = mismatches
            if identity > max_identity:
                max_identity = identity
            if min_mismatches == 0 and max_identity == 100.0:
                break  # Best possible match found
    return min_mismatches, f"{max_identity:.2f}%"

def process_chunks(chunks, input_seq, t):
    results = []
    pool = mp.Pool(processes=t)
    chunk_infos = [ (input_seq[start:end+1], input_seq, (start, end)) for start, end in chunks ]
    # Map function with multiple outputs
    mismatches_identities = pool.map(compute_min_mismatches_and_identity, chunk_infos)
    pool.close()
    pool.join()
    for (start, end), (mism, identity) in zip(chunks, mismatches_identities):
        results.append((start, end, mism, identity))
    return results

def write_size_unique_search(results, output_dir):
    output_path = os.path.join(output_dir, 'size_unique_search.txt')
    with open(output_path, 'w') as f:
        f.write("Idx_start\tIdx_end\tMism_num\tSeq_idty\n")
        for start, end, mism, identity in results:
            # Adjust indices to start from 1
            f.write(f"{start+1}\t{end+1}\t{mism}\t{identity}\n")
    logging.info(f"Unique pattern search results saved to {output_path}")

def main():
    args = parse_arguments()
    output_dir = args.out
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    setup_logging(output_dir)

    # Handle input-fa
    if not is_fasta_file(args.input_fa):
        sys.exit('Not a proper file format. Please provide a readable file format to call the "size_uniq_ptrn_search module."')
    input_fa = decompress_file(args.input_fa)

    # Handle target-fa
    if not is_fasta_file(args.target_fa):
        sys.exit('Not a proper file format. Please provide a readable file format to call the "size_uniq_ptrn_search module."')
    target_fa = decompress_file(args.target_fa)

    # Read sequences
    input_seq = read_fasta(input_fa)
    target_seq = read_fasta(target_fa)

    if len(input_seq) == 0 or len(target_seq) == 0:
        sys.exit("No sequences found in the input files.")

    input_seq_len = len(input_seq)

    # Step 3: Write Input Index
    write_input_index(input_seq_len, output_dir)

    # Step 4: Find match
    match_start, match_end, mismatches = find_match(input_seq, target_seq)
    write_match_index(match_start, match_end, mismatches, output_dir)
    logging.info(f"Input Match Number: {mismatches}")

    # Step 5: Generate chunks
    chunks = generate_chunks(input_seq_len, match_start, match_end, args.char_size)

    # Step 6: Process chunks
    results = process_chunks(chunks, input_seq, args.t)

    # Write results
    write_size_unique_search(results, output_dir)

    logging.info("Processing completed successfully.")

if __name__ == "__main__":
    main()
