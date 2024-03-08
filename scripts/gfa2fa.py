# FASTAhandler gfa to fasta singleline created by Hyungtaek Jung
# Example command: python3 gfa2fa.py --input-gfa test.gfa --out test_out.fa (w/ optional for --t cpu and --mem memory)

#!/usr/bin/env python3
import argparse
import gzip
import logging
import multiprocessing
from itertools import islice

def parse_args():
    parser = argparse.ArgumentParser(description="Convert GFA file to FASTA format.")
    parser.add_argument("--input-gfa", required=True, help="Input GFA file (compressed or uncompressed).")
    parser.add_argument("--out", required=True, help="Output FASTA file path.")
    parser.add_argument("--t", type=int, default=1, help="Number of CPUs. Default is 1.")
    parser.add_argument("--mem", type=int, default=10, help="Amount of memory in gigabytes. Default is 10.")
    return parser.parse_args()

def setup_logging():
    logging.basicConfig(filename='gfa2fa.log', level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')

def open_file(file_path, mode='rt'):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, mode)
    else:
        return open(file_path, mode)

def convert_lines_to_fasta(lines):
    fasta_lines = []
    for line in lines:
        if line.startswith('S'):
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                fasta_lines.append(f">{parts[1]}\n{parts[2]}\n")
    return fasta_lines

def chunked_iterable(iterable, size):
    it = iter(iterable)
    while True:
        chunk = list(islice(it, size))
        if not chunk:
            return
        yield chunk

def gfa_to_fasta(input_gfa, output_fasta, num_cpus, chunk_size=1000):
    with open_file(input_gfa) as f_in, open(output_fasta, 'w') as f_out:
        # Create a pool of workers
        with multiprocessing.Pool(num_cpus) as pool:
            # Process the file in chunks to reduce memory usage
            for chunk in chunked_iterable(f_in, chunk_size):
                fasta_chunks = pool.map(convert_lines_to_fasta, [chunk])
                for fasta_lines in fasta_chunks:
                    f_out.writelines(fasta_lines)

def main():
    args = parse_args()
    setup_logging()

    logging.info("Starting GFA to FASTA conversion.")
    logging.info(f"Using {args.t} CPUs and {args.mem}GB of memory.")

    gfa_to_fasta(args.input_gfa, args.out, args.t)

    logging.info("GFA to FASTA conversion completed successfully.")
    logging.info(f"FASTA file saved to {args.out}")

if __name__ == "__main__":
    main()
