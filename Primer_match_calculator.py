# Author: Wentao Chen
# Importing necessary libraries
import argparse
import pandas as pd
from tqdm import tqdm

# Function to parse the target file containing primers, each primer on a new line
def parse_target_file(target_file_path):
    with open(target_file_path, 'r') as file:
        sequences = [line.strip() for line in file if line.strip()]
    return sequences

# Function to parse a FASTA file, returns a dictionary with headers as keys and sequences as values
def parse_fasta(fasta_file):
    sequences = {}
    header = None
    seq = ""
    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    sequences[header] = seq
                header = line
                seq = ""
            else:
                seq += line
        if header:
            sequences[header] = seq
    return sequences

# Function to count matching sequences, returns a dictionary with reference headers as keys and match counts as values
def count_matching_sequences(target_sequences, reference_sequences):
    counts = {}
    for ref_header, ref_seq in tqdm(reference_sequences.items(), desc="Counting matches", unit="sequence"):
        match_counts = {}
        for target_seq in target_sequences:
            match_counts[target_seq] = ref_seq.count(target_seq)
        counts[ref_header] = match_counts
    return counts

# Setting up command line argument parser
parser = argparse.ArgumentParser(description='Count the occurrences of target sequences (primers) in a reference FASTA file.')
# Argument for the target file path
parser.add_argument('target_file', help='Path to the file containing primers to search for. Each primer should be on its own line.')
# Argument for the reference FASTA file path
parser.add_argument('reference_file', help='Path to the FASTA file to be searched.')
# Optional argument for the output file path
parser.add_argument('--output', help='Path to save the results as an xlsx file. Default is "output.xlsx".', default="output.xlsx")
args = parser.parse_args()

target_sequences = parse_target_file(args.target_file)
reference_sequences = parse_fasta(args.reference_file)

matching_sequences_counts = count_matching_sequences(target_sequences, reference_sequences)

# Creating a DataFrame from the count dictionary, preserving the order of target sequences
df = pd.DataFrame.from_dict(matching_sequences_counts, orient="index", columns=target_sequences)

print(df)

# Saving the results to an Excel file
df.to_excel(args.output)

print(f"Results saved to {args.output}")
