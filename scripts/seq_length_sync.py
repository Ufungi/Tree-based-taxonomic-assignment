import argparse

def read_fasta(file_path):
    with open(file_path, 'r') as file:
        sequences = {}
        sequence_id = None
        sequence_data = []
        
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence_id:
                    sequences[sequence_id] = ''.join(sequence_data)
                sequence_id = line
                sequence_data = []
            else:
                sequence_data.append(line)
        
        if sequence_id:
            sequences[sequence_id] = ''.join(sequence_data)
    
    return sequences

def write_fasta(sequences, file_path):
    with open(file_path, 'w') as file:
        for sequence_id, sequence in sequences.items():
            file.write(f"{sequence_id}\n")
            file.write(f"{sequence}\n")

def pad_sequences(sequences):
    max_length = max(len(seq) for seq in sequences.values())
    padded_sequences = {id_: seq + '-' * (max_length - len(seq)) for id_, seq in sequences.items()}
    return padded_sequences

def main(input_file, output_file):
    sequences = read_fasta(input_file)
    padded_sequences = pad_sequences(sequences)
    write_fasta(padded_sequences, output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pad sequences in a FASTA file to the maximum sequence length.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    
    args = parser.parse_args()
    main(args.input, args.output)
