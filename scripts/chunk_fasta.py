# -*- coding: utf-8 -*-
import argparse

def split_fasta(input_file, output_prefix, num_chunks):
    try:
        with open(input_file, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print("File '{}' not found.".format(input_file))
        return

    # Find each entry in the Fasta file (lines starting with '>')
    entries = []
    current_entry = ""
    for line in lines:
        if line.startswith('>'):
            if current_entry:
                entries.append(current_entry)
            current_entry = line
        else:
            current_entry += line
    if current_entry:
        entries.append(current_entry)

    # Calculate the chunk size for splitting
    chunk_size = len(entries) // num_chunks

    # Create file chunks
    for i in range(num_chunks):
        start_idx = i * chunk_size
        end_idx = (i + 1) * chunk_size
        chunk_entries = entries[start_idx:end_idx]

        output_file = "{}_{}.fasta".format(output_prefix, i+1)
        with open(output_file, 'w') as f:
            f.write(''.join(chunk_entries))

        print("{} file has been created.".format(output_file))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Split a Fasta file into multiple chunks.')
    parser.add_argument('-i', '--input', type=str, help='Input Fasta file path', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output file name prefix', required=True)
    parser.add_argument('-n', '--num_chunks', type=int, help='Desired number of chunks', required=True)
    args = parser.parse_args()

    split_fasta(args.input, args.output, args.num_chunks)
