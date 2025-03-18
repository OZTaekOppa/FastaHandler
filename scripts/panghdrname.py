# This script to rename the fasta header text file for pangenome, then as an input for prfxselrename.py.
# Make sure the original Fasta headers do not have any space. The space must be replaced with "_".

#!/usr/bin/env python3
import os
import gzip
import zipfile
import tarfile
import bz2
import argparse
import logging
import multiprocessing

# Set up logging
logging.basicConfig(filename='pang_hdr_name.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

def parse_args():
    parser = argparse.ArgumentParser(description='Process fasta header information and add new column based on conditions.')
    parser.add_argument('--input', required=True, help='Input file path (.txt, .fasta, .fa, .fna, .gz, .bz2, .zip, .tar.gz)')
    parser.add_argument('--output', required=True, help='Output file path')
    parser.add_argument('--t', type=int, default=1, help='Number of CPUs to use (default: 1)')
    parser.add_argument('--mem', type=int, default=10, help='Memory in GB to use (default: 10)')
    return parser.parse_args()

def validate_file_extension(file_path):
    valid_extensions = ['.txt', '.fasta', '.fa', '.fna', '.gz', '.bz2', '.zip', '.tar.gz']
    if not any(file_path.endswith(ext) for ext in valid_extensions):
        raise ValueError("Not a proper file format. Please provide a readable file format to call the pang_hdr_name module.")

def decompress_file(file_path):
    if file_path.endswith('.gz'):
        with gzip.open(file_path, 'rt') as f:
            return f.read()
    elif file_path.endswith('.bz2'):
        with bz2.open(file_path, 'rt') as f:
            return f.read()
    elif file_path.endswith('.zip'):
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            extracted_file = zip_ref.namelist()[0]
            with zip_ref.open(extracted_file, 'r') as f:
                return f.read().decode('utf-8')
    elif file_path.endswith('.tar.gz'):
        with tarfile.open(file_path, 'r:gz') as tar:
            tar.extractall('.')
            extracted_file = tar.getnames()[0]
            with open(extracted_file, 'r') as f:
                return f.read()
    else:
        with open(file_path, 'r') as f:
            return f.read()

def process_line(line):
    parts = line.split()
    if len(parts) < 3:
        return None  # Handle invalid lines gracefully

    # Step 2.1: Grep only characters in the first column
    first_col_char_only = ''.join([i for i in parts[0] if not i.isdigit()])

    # Step 2.2: Extract OriSeqID, check for h1 or h2
    ori_seq_id = parts[1].split('=')[1]
    if ori_seq_id.startswith("h1"):
        ori_seq_tag = "#1#"
    elif ori_seq_id.startswith("h2"):
        ori_seq_tag = "#2#"
    else:
        return None  # If OriSeqID doesn't start with h1 or h2, skip this line

    # Step 2.3: Extract the full OriSeqID after "OriSeqID="
    full_ori_seq = parts[1].split('=')[1]

    # Step 2.4: Combine first column character, ori_seq_tag, and full_ori_seq (no spaces)
    new_column = f"{first_col_char_only}{ori_seq_tag}{full_ori_seq}"

    # Step 2.5: Replace spaces in the first three columns with underscores and create two columns
    original_columns_with_underscores = "_".join(parts[:3])

    return f"{original_columns_with_underscores}\t{new_column}"

def process_file_content(content):
    lines = content.splitlines()
    processed_lines = [process_line(line) for line in lines if line.strip()]
    return [line for line in processed_lines if line]

def write_output_file(processed_lines, output_file):
    with open(output_file, 'w') as f:
        for line in processed_lines:
            f.write(f"{line}\n")

def process_file(input_file, output_file):
    try:
        content = decompress_file(input_file)
        processed_lines = process_file_content(content)
        write_output_file(processed_lines, output_file)
        logging.info(f"Processing complete. Output saved to {output_file}")
    except Exception as e:
        logging.error(f"Error processing file {input_file}: {e}")
        print(f"Error processing file {input_file}: {e}")

def main():
    args = parse_args()

    # Validate file extension
    try:
        validate_file_extension(args.input)
    except ValueError as e:
        logging.error(e)
        print(e)
        return

    # Setup multiprocessing
    cpu_count = min(args.t, multiprocessing.cpu_count())
    memory_limit = args.mem

    logging.info(f"Using {cpu_count} CPUs and {memory_limit} GB memory for processing.")

    # Process file
    process_file(args.input, args.output)

if __name__ == "__main__":
    main()
