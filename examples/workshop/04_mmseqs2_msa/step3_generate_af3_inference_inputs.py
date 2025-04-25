import csv
from pathlib import Path
from af3cli import InputBuilder, ProteinSequence, SMILigand

# Constants
INPUT_FILE = "../data/preliminaryTests.tsv"
MSA_INPUT_DIR = Path("colabfold_msa_outputs")
OUTPUT_DIR = Path("af3_inference_inputs")
OUTPUT_DIR.mkdir(exist_ok=True)

# Read the TSV file
with open(INPUT_FILE, newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        # Extract data from each row
        symbol = row["Symbol"]
        compound = row["Compound"]
        sequence_str = row["ProteinSequence"]
        smiles_str = row["SMILES"]

        # Construct job name and output filename
        job_name = f"{symbol}_{compound}"
        output_filename = OUTPUT_DIR / f"{job_name}.json"
        msa_filename = MSA_INPUT_DIR / f"{job_name}/{job_name}.a3m"

        # Build AF3 input with MSA from colabfold
        msa = MSA(unpaired=msa_filename, unpaired_is_path=True)
        sequence = ProteinSequence(seq_str=sequence_str, msa=msa)
        ligand = SMILigand(ligand_value=smiles_str)
        input_builder = InputBuilder()
        input_builder.set_name(job_name)
        input_builder.add_sequence(sequence)
        input_builder.add_ligand(ligand)
        input_builder.set_version(3)
        input_builder.set_seeds([1])
        internal_input = input_builder.build()

        # Write to file
        internal_input.write(output_filename)
