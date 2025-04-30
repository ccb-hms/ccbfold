---
title: "AlphaFold3 Workshop"
---
# AlphaFold3 Workshop

AlphaFold3 is blah blah blah and way better than AlphaFold2 because of blah blah blah

### Setup

We start by creating a reproducible python environment with all the necessary packages. To do so, we first install [uv](https://docs.astral.sh/uv/), an extremely fast Python project manager written in Rust. We then use `uv` to restore the virtual environment for `ccbfold`, a python package created for the development of the workshop:

```bash
# install uv package manager
curl -LsSf https://astral.sh/uv/install.sh | sh

# add to PATH
source $HOME/.local/bin/env

# download ccbfold repo
git clone https://github.com/ccb-hms/ccbfold.git
cd ccbfold

# restore venv and activate
uv sync
source .venv/bin/activate

# create and move to directory for workshop
mkdir ../af3_workshop && cd $_
```

### Exercise 1: Hello Alphafold3

Let's start off with a simple example to get started. This uses the amino acid sequence for a DNA binding protein and gene editing meganuclease as well as the associated DNA sequence. For this complex, the crystal structure has been experimentally determined using x-ray diffraction and deposited to [PDB entry 7RCE](https://www.rcsb.org/structure/7RCE). Let's create the input file that we will need for AlphaFold3:

```bash
mkdir 01_basic
cd 01_basic
nano pdb_7rce_input.json

# to paste: Ctrl (Cmd on Mac) + Shift + v
# to save: Ctrl (Cmd on Mac) + o then Enter
# to exit nano: Ctrl (Cmd on Mac) + x
```

Save the following into the newly created `pdb_7rce_input.json`:

```json
{
    "name": "Protein-DNA-Ion: PDB 7RCE",
    "modelSeeds": [1],
    "sequences": [
      {
        "protein": {
          "id": "A",
          "sequence": "MASSRRESINPWILTGFADAEGSFGLSILNRNRGTARYHTRLSFTIMLHNKDKSILENIQSTWKVGSILNNGDHYVSLVVYRFEDLKVIIDHFEKYPLITQKLGDYKLFKQAFSVMENKEHLKENGIKELVRIKAKMNWGLNDELKKAFPENISKERPLINKNIPNFKWLAGFTSGDGSFFVRLRKSNVNARVRVQLVFEISQHIRDKNLMNSLITYLGCGHIYEGNKSERSWLQFRVEKFSDINDKIIPVFQENTLIGVKLEDFEDWCKVAKLIEEKKHLTESGLDEIKKIKLNMNKGR"
        }
      },
      {
        "dna": {
          "id": "B",
          "sequence": "GGGGGCATGCAGATCCCACAGGCGCG"
        }
      },
      {
        "ligand": {
          "id": ["C", "D", "E"],
          "ccdCodes": ["CA"]
        }
      },
      {
        "ligand": {
          "id": "F",
          "ccdCodes": ["NA"]
        }
      }
    ],
    "dialect": "alphafold3",
    "version": 2
}
```


We also need a sbatch script in order to submit our input to the Slurm workload manager on O2:

```bash
nano pdb_7rce_submit.sh
```

Paste and save the following:

```bash
#!/bin/bash
#SBATCH -c 20                    # Request 20 cores
#SBATCH --mem=64G                # Memory total in GiB
#SBATCH --partition=gpu_quad     # Partition to run in
#SBATCH -o pdb_7rce_%j.out       # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e pdb_7rce_%j.err       # File to which STDERR will be written, including job ID (%j)
#SBATCH --gres=gpu:l40s:1        # GPU requested
#SBATCH -t 0-01:00               # Runtime in D-HH:MM format

module load alphafold/3.0.1

run_alphafold.py \
   --json_path=pdb_7rce_input.json \
   --output_dir=pdb_7rce_output
```

We are now ready to submit our first AlphaFold3 job to the cluster:

```bash
sbatch pdb_7rce_submit.sh
```

### Exercise 2: Running a Batch of Jobs

If we have a lot of structures that we want to predict, seting up the above manually can be quite tedious and error prone. To make it a lot easier, let's take advantage of the python package `af3cli` to generate our AlphaFold3 input JSON files programatically. The input sequences and ligands will be read in from the following TSV: 

|   Symbol   |   ProteinSequence                            |   Compound  |   SMILES                                        |
|------------|----------------------------------------------|-------------|-------------------------------------------------|
|   RREB1    | MTSSSPAGLEGSDLSSINTMMSAVMSVGKVTENGGSPQ... |   E1879     |   Cn1c(CNC=C(C#N)C(=O)Nc2ccccc2)cc(=O)n(C)c1=O  |
|   SEPTIN5  | MSTGLRYKSKLATPEDKQDIDKQYVGFATLPNQVHRKSV... |   E2964     |   CN(C)c1nc(ncc1NC(=O)CCl)N2CCOCC2              |
|   ATP6V1A  | MDFSKLPKILDEDKESTFGYVHGVSGPVVTACDMAGAAM... |   CL117     |   NC(=O)C=1C=CC=CC1NC(=O)CCl                    |
|   ALDH1A2  | MTSSKIEMPGEVKADPAALMASLHLLPSPTPNLEIKYTKIF... |   CL102     |   COC=1C=CC=C(CN(C)C(=O)CCl)C1                  |
|   ACAT1    | MAVLAALLRSGARSRSPLLRRLVQEIRYVERSYVSKPTLKE... |   E1387     |   FC(F)Sc1ccc(NC(=O)CCl)cc1                     |
|   ACAT1    | MAVLAALLRSGARSRSPLLRRLVQEIRYVERSYVSKPTLKE... |   E2992     |   FC(F)(F)CC(=O)N1CCN(CC1)C(=O)CCl              |
|   ACAT1    | MAVLAALLRSGARSRSPLLRRLVQEIRYVERSYVSKPTLKE... |   CL56      |   CC(=O)NC=1C=CC(NC(=O)CCl)=CC1                 |
|   ACAT1    | MAVLAALLRSGARSRSPLLRRLVQEIRYVERSYVSKPTLKE... |   E1384     |   CSc1cccc(NC(=O)CCl)c1                         |
|   ACAT1    | MAVLAALLRSGARSRSPLLRRLVQEIRYVERSYVSKPTLKE... |   E1383     |   FC(F)Oc1ccc(NC(=O)CCl)cc1                     |
|   ATP6V1A  | MDFSKLPKILDEDKESTFGYVHGVSGPVVTACDMAGAAM... |   E2912     |   Fc1ccc2N(C(=O)CCl)C3(CCC3)C(=O)Nc2c1          |

We want to predit the structure of each `ProteinSequence` with the associated small molecule compound, whose chemical formula is provided as a `SMILES` string. 

Let's start by setting up a directory for this exercise and creating our `af3cli` script:

```bash
# go up one directory if still in 01_basic/
cd ..
mkdir 02_batch_of_jobs
cd 02_batch_of_jobs
nano step1_generate_af3_inputs.py
```

paste the following into this file:

```python
import csv
from pathlib import Path
from af3cli import InputBuilder, ProteinSequence, SMILigand

# Constants
INPUT_FILE = "../data/preliminaryTests.tsv"
OUTPUT_DIR = Path("af3_inputs")
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
        filename = OUTPUT_DIR / f"{job_name}.json"

        # Build AF3 input
        sequence = ProteinSequence(seq_str=sequence_str)
        ligand = SMILigand(ligand_value=smiles_str)
        input_builder = InputBuilder()
        input_builder.set_name(job_name)
        input_builder.add_sequence(sequence)
        input_builder.add_ligand(ligand)
        input_builder.set_version(2)
        input_builder.set_seeds([1])
        internal_input = input_builder.build()

        # Write to file
        internal_input.write(filename)

```

Before we can run the above, we first need to download the data. Since the data will also be used for the next exercise, lets put it somewhere less specific:

```bash
mkdir ../data
curl -L -o ../data/preliminaryTests.tsv "https://raw.githubusercontent.com/ccb-hms/ccbfold/refs/heads/main/examples/workshop/data/preliminaryTests.tsv"
```

We can now generate our 10 AlphaFold3 inputs:

```python
# will save .json files to af3_inputs/
python step1_generate_af3_inputs.py
```

To submit our inputs, we can use a single sbatch script with the `#SBATCH --array=0-9` declaration:

```bash
nano step2_submit_batch.sh
```

Save the following to this file:

```bash
#!/bin/bash
#SBATCH -c 20                    # Request 20 cores
#SBATCH --mem=64G                # Memory total in GiB
#SBATCH --partition=gpu_quad     # Partition to run in
#SBATCH -o logs/af3_job_%A_%a.out       # STDOUT file
#SBATCH -e logs/af3_job_%A_%a.err       # STDERR file
#SBATCH --gres=gpu:l40s:1        # GPU requested
#SBATCH -t 0-01:00               # Runtime in D-HH:MM
#SBATCH --array=0-9              # Job array indices (for 10 .json files)

module load alphafold/3.0.1

# Get the input json file
INPUT_DIR="af3_inputs"
INPUT_FILES=(${INPUT_DIR}/*.json)
INPUT_FILE=${INPUT_FILES[$SLURM_ARRAY_TASK_ID]}

# Derive job name and output directory from file name
JOB_NAME=$(basename "$INPUT_FILE" .json)
OUTPUT_DIR="af3_outputs"

# Make sure output directory exists
mkdir -p "$OUTPUT_DIR"

run_alphafold.py \
   --json_path="$INPUT_FILE" \
   --output_dir="$OUTPUT_DIR"
```

We submit our sbatch script to the Slurm scheduler just as before:

```bash
sbatch step2_submit_batch.sh
```

### Exercise 3: Separate MSA and Inference

For a typical AlphaFold3 run, the multiple sequence alignment (MSA) accounts for the vast majority of the total run-time (~5/6th of total). The MSA step is performed using `jackhmmer`, which runs entirely using CPU. As a result, we can save substantial GPU time by running MSA and inference as separate steps. To setup, we use the same starting input files as for exercise 2:

```bash
cd ..
mkdir 03_separate_msa
cd 03_separate_msa

# copy over af3_inputs from exercise 2
cp -R ../02_batch_of_jobs/af3_inputs .
```

Our initial submission script is also identical except for the following:
- added `--norun_data_pipeline` flag so that AlphaFold3 only runs MSA
- save our MSA results to `af3_msa_outputs`
- we change our partition and remove our request for a GPU

 Save the following to `step2_submit_msa_batch.sh`:

```bash
#!/bin/bash
#SBATCH -c 20                    # Request 20 cores
#SBATCH --mem=64G                # Memory total in GiB
#SBATCH --partition=short        # Partition to run in
#SBATCH -o logs/af3_msa_job_%A_%a.out       # STDOUT file
#SBATCH -e logs/af3_msa_job_%A_%a.err       # STDERR file
#SBATCH --gres=gpu:l40s:1        # GPU requested
#SBATCH -t 0-01:00               # Runtime in D-HH:MM
#SBATCH --array=0-9              # Job array indices (for 10 .json files)

module load alphafold/3.0.1

# Get the input json file
INPUT_DIR="af3_inputs"
INPUT_FILES=(${INPUT_DIR}/*.json)
INPUT_FILE=${INPUT_FILES[$SLURM_ARRAY_TASK_ID]}

# Derive job name and output directory from file name
JOB_NAME=$(basename "$INPUT_FILE" .json)
OUTPUT_DIR="af3_msa_outputs"
JOB_OUTPUT_DIR="${OUTPUT_DIR}/${JOB_NAME}"

# Make sure output directory exists
mkdir -p "$OUTPUT_DIR"

# Run AlphaFold3 MSA only
run_alphafold.py \
   --json_path="$INPUT_FILE" \
   --output_dir="$JOB_OUTPUT_DIR" \
   --norun_inference
```

We submit our MSA just as before:

```bash
sbatch step2_submit_msa_batch.sh
```

Once the above completes we can then submit our inference as well. Save the following to `step3_submit_inference_batch.sh`:

```bash
#!/bin/bash
#SBATCH -c 2                     # Request 4 cores
#SBATCH --mem=16G                # Memory total in GiB
#SBATCH --partition=gpu_quad     # Partition to run in
#SBATCH -o logs/af3_inference_job_%A_%a.out       # STDOUT file
#SBATCH -e logs/af3_inference_job_%A_%a.err       # STDERR file
#SBATCH --gres=gpu:l40s:1        # GPU requested
#SBATCH -t 0-01:00               # Runtime in D-HH:MM
#SBATCH --array=0-9              # Job array indices (for 10 .json files)

module load alphafold/3.0.1

# Get the input json file with MSA created in step 2
INPUT_DIR="af3_msa_outputs"
JOB_NAMES=($(ls "${INPUT_DIR}"))
JOB_NAME=${JOB_NAMES[$SLURM_ARRAY_TASK_ID]}
INPUT_FILE="${INPUT_DIR}/${JOB_NAME}/${JOB_NAME}_data.json"

# Derive job name and output directory from file name
OUTPUT_DIR="af3_inference_outputs"
JOB_OUTPUT_DIR="${OUTPUT_DIR}/${JOB_NAME}"

# Make sure output directory exists
mkdir -p "$OUTPUT_DIR"

# Run AlphaFold3 Inference only
run_alphafold.py \
   --json_path="$INPUT_FILE" \
   --output_dir="$JOB_OUTPUT_DIR" \
   --norun_data_pipeline
```

Note the following changes:
- added `--norun_data_pipeline` flag so AlphaFold3 only runs inference step
- our `INPUT_FILE` for each submission is the output of the MSA step instead of the original JSON
- we save our final results to `af3_inference_outputs`
- we don't request much CPU or RAM, as inference is mostly GPU-bound

We can submit our jobs as we always do. We can also add a dependency requiring that our previous job (the MSA) finished first:

```bash
# TODO: add dependency
sbatch step3_submit_inference_batch.sh
```


### Exercise 4: Perform MSA with GPU mmseqs2

