import csv

# Constants
INPUT_FILE = "../data/preliminaryTests.tsv"
OUTPUT_FILE = "../data/preliminaryTests.fasta"

with open(INPUT_FILE, 'r') as input_tsv, open(OUTPUT_FILE, 'w') as output_fasta:
    input_reader = csv.DictReader(input_tsv, delimiter='\t')
    for row in input_reader:

        # Extract data from each row
        symbol = row["Symbol"]
        compound = row["Compound"]
        sequence_str = row["ProteinSequence"]

        # Add sequence to output fasta
        output_fasta.write(f">{symbol}_{compound}\n")
        output_fasta.write(f"{sequence_str}\n")
