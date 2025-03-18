# This script to make sequence partition based on the name of CL in the fasta header.
# Example usage: python3 chr_pansn_ext.py --input-fa FAfile --output Outfolder (--t and --mem optional)

#!/usr/bin/env python3
import os
import gzip
import argparse
from multiprocessing import Pool, cpu_count
from Bio import SeqIO
import logging

# Set up logging
logging.basicConfig(filename='chr_pansn_ext.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

def parse_args():
    parser = argparse.ArgumentParser(description='Separate fasta sequences into groups based on CL: prefix in headers.')
    parser.add_argument('--input-fa', required=True, help='Input fasta file path (*.fasta, *.fa, *.fasta.gz, *.fa.gz)')
    parser.add_argument('--output', required=True, help='Output directory for separated fasta files')
    parser.add_argument('--t', type=int, default=1, help='Number of CPUs to use (default: 1)')
    parser.add_argument('--mem', type=int, default=10, help='Memory in GB to use (default: 10)')
    return parser.parse_args()

def validate_file_extension(file_path):
    valid_extensions = ['.fasta', '.fa', '.fasta.gz', '.fa.gz']
    if not any(file_path.endswith(ext) for ext in valid_extensions):
        raise ValueError("Not a proper file format. Please provide a readable file format to call the chr_pansn_ext module.")

def read_fasta(file_path):
    if file_path.endswith('.gz'):
        with gzip.open(file_path, 'rt') as f:
            return list(SeqIO.parse(f, 'fasta'))
    else:
        with open(file_path, 'r') as f:
            return list(SeqIO.parse(f, 'fasta'))

def group_sequences(sequences):
    groups = {}
    for record in sequences:
        header = record.description
        if 'CL:' in header:
            cl_prefix = header.split('CL:')[1].split('_')[0]
            if cl_prefix not in groups:
                groups[cl_prefix] = []
            groups[cl_prefix].append(record)
        else:
            logging.warning(f"No 'CL:' found in header: {header}")
    return groups

def write_fasta_groups(groups, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for group, records in groups.items():
        output_file = os.path.join(output_dir, f"{group}.fasta")
        with open(output_file, 'w') as f:
            SeqIO.write(records, f, 'fasta')
        logging.info(f"Wrote {len(records)} sequences to {output_file}")

def process_fasta(args):
    input_fa, output_dir = args
    sequences = read_fasta(input_fa)
    groups = group_sequences(sequences)
    write_fasta_groups(groups, output_dir)

def main():
    args = parse_args()

    try:
        validate_file_extension(args.input_fa)
    except ValueError as e:
        logging.error(e)
        print(e)
        return

    cpu_count_to_use = min(args.t, cpu_count())
    memory_to_use = args.mem
    logging.info(f"Using {cpu_count_to_use} CPUs and {memory_to_use} GB of memory")

    pool = Pool(processes=cpu_count_to_use)
    try:
        pool.map(process_fasta, [(args.input_fa, args.output)])
        pool.close()
        pool.join()
    except Exception as e:
        logging.error(f"Error during multiprocessing: {e}")
        print(f"Error during multiprocessing: {e}")
        return

if __name__ == "__main__":
    main()
