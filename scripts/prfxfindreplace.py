# FastaHandler rename only selected ID pattern like "sed" command created by Hyungtaek Jung
# A multi-fasta (multiline) file to a single-line fasta with replacing ID/Prefix name files only in the input id file
# Example command: python3 prfxfindreplace.py --input-seq test.fa --find-ptrn hifiasm --replace Assembly --out test_out.fa (w/ optional for --t cpu and --mem memory)

#!/usr/bin/env python3

import argparse
import os
import re
import gzip
import logging
from Bio import SeqIO
import multiprocessing

# Set up logging
logging.basicConfig(filename='prfxfindreplace.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

def unzip_file(input_file, output_dir):
    # Check compressed file and extract it
    if input_file.lower().endswith(('.gz', '.tar.gz', '.zip', '.bz2')):
        with gzip.open(input_file, 'rb') as f_in:
            with open(os.path.join(output_dir, os.path.basename(input_file)[:-3]), 'wb') as f_out:
                f_out.write(f_in.read())
        return os.path.join(output_dir, os.path.basename(input_file)[:-3])
    return input_file

def process_sequence(record, find_pattern, replace_pattern):
    new_id = re.sub(find_pattern, replace_pattern, record.id)
    new_description = re.sub(find_pattern, replace_pattern, record.description)
    record.id = new_id
    record.description = new_description
    return record

def replace_pattern_in_fasta(input_fasta, find_pattern, replace_pattern, output_fasta, threads):
    input_fasta = unzip_file(input_fasta, os.path.dirname(input_fasta))
    temp_output_fasta = output_fasta + ".temp"
    with open(input_fasta, "r") as input_handle, open(temp_output_fasta, "w") as temp_output_handle:
        records = SeqIO.parse(input_handle, "fasta")
        pool = multiprocessing.Pool(processes=threads)
        modified_records = pool.starmap(process_sequence, [(record, find_pattern, replace_pattern) for record in records])
        SeqIO.write(modified_records, temp_output_handle, "fasta")
    # Convert multiline sequences to single-line format
    multi2singleline(temp_output_fasta, output_fasta)
    os.remove(temp_output_fasta)
    logging.info(f"Modified fasta file saved as {output_fasta}")

def multi2singleline(input_file, output_file):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        header = None
        sequence = ""
        for line in f_in:
            if line.startswith(">"):
                if header is not None:
                    f_out.write(header + "\n")
                    f_out.write(sequence + "\n")
                header = line.strip()
                sequence = ""
            else:
                sequence += line.strip()
        if header is not None and sequence:
            f_out.write(header + "\n")
            f_out.write(sequence + "\n")

def main():
    parser = argparse.ArgumentParser(description="Replace a pattern in fasta headers with a user-defined string.")
    parser.add_argument("--input-seq", required=True, help="Input fasta file path.")
    parser.add_argument("--find-ptrn", required=True, help="Pattern to find in fasta headers.")
    parser.add_argument("--replace", required=True, help="String to replace the pattern with.")
    parser.add_argument("--out", required=True, help="Output fasta file path.")
    parser.add_argument("--t", type=int, default=1, help="Number of threads for multiprocessing.")
    parser.add_argument("--mem", type=int, default=10, help="Amount of memory in gigabytes.")
    args = parser.parse_args()

    replace_pattern_in_fasta(args.input_seq, args.find_ptrn, args.replace, args.out, args.t)

if __name__ == "__main__":
    main()
