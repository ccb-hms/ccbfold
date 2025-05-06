---
title: "AlphaFold3 Workshop"
---
# AlphaFold3 on O2 Workshop

Welcome to the AlphaFold3 Workshop!

AlphaFold3 significantly advances protein structure prediction beyond AlphaFold2, offering greater accuracy and broader capabilities.

### Setup

To begin, we'll set up a reproducible Python environment with all required packages. We’ll use [uv](https://docs.astral.sh/uv/), a fast Python project manager written in Rust, to restore the virtual environment for `ccbfold`—a custom Python package developed for this workshop.

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
mkdir ../af3_workshop
cd ../af3_workshop
```

### Exercise 1: Hello Alphafold3


To get started, we'll use a simple example involving a DNA-binding protein and gene-editing meganuclease, along with its the DNA sequence that it binds to. A crystal structure for this complex is available in [PDB entry 7RCE](https://www.rcsb.org/structure/7RCE).

First, create the necessary input file for AlphaFold3:

```bash
mkdir 01_basic
cd 01_basic
nano pdb_7rce_input.json

# to paste: Ctrl (Cmd on Mac) + Shift + v
# to save: Ctrl (Cmd on Mac) + o then Enter
# to exit nano: Ctrl (Cmd on Mac) + x
#
# to check file contents: cat pdb_7rce_input.json
```

Paste the following content into the file:

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


Next, create a Slurm submission script to run AlphaFold3:

```bash
nano pdb_7rce_submit.sh
```

Paste in the following:

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

Now you're ready to submit your first AlphaFold3 prediction job:

```bash
# NOTE: may need to specify account, e.g.:
# sbatch -A zaklab pdb_7rce_submit.sh
sbatch pdb_7rce_submit.sh
```

### Exercise 2: Running a Batch of Jobs

If you want to generate structures for many different protein-ligand pairs, setting up input files one at a time quickly becomes tedious. Instead, we can use a Python utility called `af3cli` to programmatically generate input JSONs.

All of the sequences and corresponding small molecules are stored in a tab-separated file like this:

|   Symbol   |   ProteinSequence                            |   Compound  |   SMILES                                        |
|------------|----------------------------------------------|-------------|-------------------------------------------------|
|   RREB1    | MTSSSPAGLEGSDLSSINTMMSAVMSVGKVTENGGSPQ... |   E1879     |   Cn1c(CNC=C(C#N)C(=O)Nc2ccccc2)cc(=O)n(C)c1=O  |
|   SEPTIN5  | MSTGLRYKSKLATPEDKQDIDKQYVGFATLPNQVHRKSV... |   E2964     |   CN(C)c1nc(ncc1NC(=O)CCl)N2CCOCC2              |
|   ATP6V1A  | MDFSKLPKILDEDKESTFGYVHGVSGPVVTACDMAGAAM... |   CL117     |   NC(=O)C=1C=CC=CC1NC(=O)CCl                    |
|   ...  | ... |   ...     |   ...                  |

Our job is to predict the structure of each protein sequence in complex with its associated small molecule, using the provided SMILES string as the ligand.

Start by setting up a directory and writing the script that converts the above TSV into AlphaFold3-compatible inputs:

```bash
# go up one directory if still in 01_basic/
cd ..
mkdir 02_batch_of_jobs
cd 02_batch_of_jobs
nano step1_generate_af3_inputs.py
```

Paste the following script:

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
        input_builder = InputBuilder()
        sequence = ProteinSequence(seq_str=sequence_str)
        ligand = SMILigand(ligand_value=smiles_str)
        input_builder.set_name(job_name)
        input_builder.add_sequence(sequence)
        input_builder.add_ligand(ligand)
        input_builder.set_version(2)
        input_builder.set_seeds([1])
        internal_input = input_builder.build()

        # Write to file
        internal_input.write(filename)

```

Before running this script, download the input dataset:

```bash
mkdir ../data
curl -L -o ../data/preliminaryTests.tsv "https://raw.githubusercontent.com/ccb-hms/ccbfold/refs/heads/main/examples/workshop/data/preliminaryTests.tsv"
```

Now generate the input JSON files:

```python
# will save .json files to af3_inputs/
python step1_generate_af3_inputs.py
```

To submit all of the jobs to the cluster at once, we’ll take advantage of Slurm’s job arrays:


```bash
nano step2_submit_batch.sh
```

Paste the following script:

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

# Make sure output directory exists
OUTPUT_DIR="af3_outputs"
mkdir -p "$OUTPUT_DIR"

run_alphafold.py \
   --json_path="$INPUT_FILE" \
   --output_dir="$OUTPUT_DIR"
```

And submit the batch job array like so:

```bash
sbatch step2_submit_batch.sh
```

### Exercise 3: Separate MSA and Inference

When running AlphaFold3, the majority of execution time — often more than 80% — is spent creating the multiple sequence alignment (MSA), a step that only uses CPUs. To better utilize GPU resources, we can split the run into two phases: first run MSA, and then run inference separately.

Set up the directory and copy over inputs from the previous exercise:

```bash
cd ..
mkdir 03_separate_msa
cd 03_separate_msa

# copy over af3_inputs from exercise 2
cp -R ../02_batch_of_jobs/af3_inputs .
```

Create a Slurm script to run only the MSA step:

```bash
nano step2_submit_msa_batch.sh
```

Save the following:

```bash
#!/bin/bash
#SBATCH -c 20                    # Request 20 cores
#SBATCH --mem=64G                # Memory total in GiB
#SBATCH --partition=short        # Partition to run in
#SBATCH -o logs/af3_msa_job_%A_%a.out       # STDOUT file
#SBATCH -e logs/af3_msa_job_%A_%a.err       # STDERR file
#SBATCH -t 0-01:00               # Runtime in D-HH:MM
#SBATCH --array=0-9              # Job array indices (for 10 .json files)

module load alphafold/3.0.1

# Get the input json file
INPUT_DIR="af3_inputs"
INPUT_FILES=(${INPUT_DIR}/*.json)
INPUT_FILE=${INPUT_FILES[$SLURM_ARRAY_TASK_ID]}


# Make sure output directory exists
OUTPUT_DIR="af3_msa_outputs"
mkdir -p "$OUTPUT_DIR"

# Run AlphaFold3 MSA only
run_alphafold.py \
   --json_path="$INPUT_FILE" \
   --output_dir="$OUTPUT_DIR" \
   --norun_inference
```

Note that the submission script has the following changes:
- adds `--norun_data_pipeline` flag so that AlphaFold3 only runs MSA
- we save the MSA results to `af3_msa_outputs`
- we no longer need a GPU so we use the `short` partition and remove the `--gres` flag

Go ahead and submit your MSA job:

```bash
sbatch step2_submit_msa_batch.sh
#Submitted batch job <JOB_ID_HERE>
```

Let's save the job ID from the MSA step so we can make the inference step wait until it has successfully completed:

```bash
MSA_JOB_ID=<JOB_ID_HERE>
```

Now for the inference step, we need another batch script that picks up where the MSA step left off:

```bash
nano step3_submit_inference_batch.sh
```

Paste this:

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

You can now submit the inference step, making sure to run it only after the MSA jobs complete:

```bash
sbatch --dependency=afterok:$MSA_JOB_ID step3_submit_inference_batch.sh
```

### Exercise 4: Run MSAs on GPU with mmseqs2

Previously, we split MSA and inference execution to reduce wasted GPU time. In this exercise, we’ll take it a step further and use `mmseqs2`, which includes GPU acceleration for MSA generation, allowing throughput improvements of ~10x as compared to `jackhmmer`.

We’ll start by generating a `.fasta` file with all protein sequences:

```bash
cd ..
mkdir 04_mmseqs2_msa
cd 04_mmseqs2_msa
nano step1_generate_colabfold_msa_input.py
```

Paste the following python script into this file and save it:

```python
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

```

Run it to create the `preliminaryTests.fasta` file:

```bash
python step1_generate_colabfold_msa_input.py
```

Now create the submission script to run `mmseqs2` on GPU:

```bash
nano step2_submit_colabfold_msa.sh
```

Paste this:

```bash
#!/bin/bash
#SBATCH -c 20                    # Request 20 cores
#SBATCH --mem=64G                # Memory total in GiB
#SBATCH --partition=gpu_quad     # Partition to run in
#SBATCH -o logs/colabfold_msa_job_%A_%a.out       # STDOUT file
#SBATCH -e logs/colabfold_msa_job_%A_%a.err       # STDERR file
#SBATCH --gres=gpu:l40s:1        # GPU requested
#SBATCH -t 0-01:00               # Runtime in D-HH:MM

module use /n/shared_db/tmp/module
module load colabfold
module load mmseqs2

OUTPUT_DIR="colabfold_msa_outputs"

# Shared database directory for colabfold MSAs
DATABASE_DIR="/n/shared_db/tmp/colabfold"

# Run MSAs with mmseqs2 on GPU
colabfold_search \
   --gpu 1 \
   ../data/preliminaryTests.fasta \
   $DATABASE_DIR \
   $OUTPUT_DIR
```

Submit the job:

```bash
sbatch step2_submit_colabfold_msa.sh
# Submitted batch job <JOB_ID_HERE>

# save job id for later
MSA_JOB_ID=<JOB_ID_HERE>
```

Then, generate AlphaFold3 input files using the MSAs from `mmseqs2`:

```bash
nano step3_generate_af3_inference_inputs.py
```

Paste and save:

```python
import csv
from pathlib import Path
from af3cli import InputBuilder, ProteinSequence, SMILigand, MSA

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
        msa_filename = MSA_INPUT_DIR / f"{job_name}.a3m"

        # Build AF3 input with MSA from colabfold
        msa = MSA(unpaired=str(msa_filename), unpaired_is_path=True)
        sequence = ProteinSequence(seq_str=sequence_str, msa=msa)
        ligand = SMILigand(ligand_value=smiles_str)
        input_builder = InputBuilder()
        input_builder.set_name(job_name)
        input_builder.add_sequence(sequence)
        input_builder.add_ligand(ligand)
        input_builder.set_version(2)
        input_builder.set_seeds([1])
        internal_input = input_builder.build()

        # Write to file
        internal_input.write(output_filename)
```

Run the script to generate the input JSONs:

```bash
python step3_generate_af3_inference_inputs.py
```

Now prepare the batch inference submission script:

```bash
nano step4_submit_af3_inference_batch.sh
```

Paste the following:

```bash
#!/bin/bash
#SBATCH -c 20                    # Request 20 cores
#SBATCH --mem=64G                # Memory total in GiB
#SBATCH --partition=gpu_quad     # Partition to run in
#SBATCH -o logs/af3_inference_job_%A_%a.out       # STDOUT file
#SBATCH -e logs/af3_inference_job_%A_%a.err       # STDERR file
#SBATCH --gres=gpu:l40s:1        # GPU requested
#SBATCH -t 0-01:00               # Runtime in D-HH:MM
#SBATCH --array=0-9              # Job array indices (for 10 .json files)

module load alphafold/3.0.1

# Get the input json file
INPUT_DIR="af3_inference_inputs"
INPUT_FILES=(${INPUT_DIR}/*.json)
INPUT_FILE=${INPUT_FILES[$SLURM_ARRAY_TASK_ID]}

# Derive job name and output directory from file name
JOB_NAME=$(basename "$INPUT_FILE" .json)
OUTPUT_DIR="af3_inference_outputs"

# Make sure output directory exists
mkdir -p "$OUTPUT_DIR"

run_alphafold.py \
   --json_path="$INPUT_FILE" \
   --output_dir="$OUTPUT_DIR"
```

Again, we submit ensuring the MSA step finishes first:

```bash
sbatch --dependency=afterok:$MSA_JOB_ID step4_submit_af3_inference_batch.sh
```


### Exercise 4: Perform MSA with GPU mmseqs2

## Exercise 5: Reading, seeing, and understanding the results

At the end of the inference step, you should have the following files in the results folder.
```
.
├── TERMS_OF_USE.md
├── protein-dna-ion_pdb_7rce_confidences.json
├── protein-dna-ion_pdb_7rce_data.json
├── protein-dna-ion_pdb_7rce_model.cif
├── protein-dna-ion_pdb_7rce_summary_confidences.json
├── ranking_scores.csv
├── seed-2530573219_sample-0
├── seed-2530573219_sample-1
├── seed-2530573219_sample-2
├── seed-2530573219_sample-3
└── seed-2530573219_sample-4
```
These result files should contain the <"name"> from the input json. The actual structure is contained in a macromolecular Crystallographic Information File (.mmCIF or .cif), a standard text file format used for representing crystallographic information. Since 2014, this is the default file format for the [Protein Data Bank](https://www.wwpdb.org/). The top-scoring model and confidence is in the root directory. Sub-directories contain the results for each sample and seed the rank of each is described in `ranking_scores.csv`.

### Visualizing the top-ranking predicted model structure

The mmcif file format is incredibly powerful. There is not enough time to go into the details of this format but know it is well supported by visualization applications such as Jmol, Chimera, and OpenRasMol and structure determination systems such as CCP4 and Phenix.

We will use a simple WebGL WebGL protein viewer widget to quickly render the structure


```python
import nglview as nv
filepath = "result/protein-dna-ion_pdb_7rce/protein-dna-ion_pdb_7rce_model.cif"
res_7rce = nv.show_file(filepath, ext='cif')
res_7rce
```


    



    NGLWidget()


As a demonstration of structure analysis possible in python, we will download the cif file of the corresponding PDB crystallographic entry and superpose the experimentally derived model on the predicted model. First we will retreive the 7RCE model with biopython.


```python
from Bio.PDB import PDBList

# Create an instance of the PDBList class
pdb_list = PDBList()

# Specify the PDB ID of the structure you want to download
pdb_id = "7rce" 
# Download the MMCIF file using the retrieve_pdb_file method
pdb_filename = pdb_list.retrieve_pdb_file(pdb_id, pdir="data/PDB_files", file_format="mmCif")

# Print the name of the downloaded file
print(pdb_filename)

```

    Structure exists: 'data/PDB_files/7rce.cif' 
    data/PDB_files/7rce.cif


Using the superpose approximation from NGviewer js library


```python
from Bio.PDB import MMCIFParser
parser = MMCIFParser()
#parse mmcif files
structure1 = parser.get_structure("prediction","result/protein-dna-ion_pdb_7rce/protein-dna-ion_pdb_7rce_model.cif")
structure2 = parser.get_structure("source","data/PDB_files/7rce.cif")

#construct the superpose view
view = nv.NGLWidget()
ref = nv.BiopythonStructure(structure1)
pred = nv.BiopythonStructure(structure2)
view.add_structure(ref)
view.update_cartoon(color = "coral", component=0)
view.add_structure(pred)
view.update_cartoon(color = "teal", component=1)
view._js("""
this.superpose(0, 1, true, ":A", ":A") 
""")
view
```

    /home/daedalus/projects/ccbfold/.venv/lib/python3.12/site-packages/Bio/PDB/StructureBuilder.py:100: PDBConstructionWarning: WARNING: Chain A is discontinuous at line 3113.
      warnings.warn(
    /home/daedalus/projects/ccbfold/.venv/lib/python3.12/site-packages/Bio/PDB/StructureBuilder.py:100: PDBConstructionWarning: WARNING: Chain B is discontinuous at line 3123.
      warnings.warn(
    /home/daedalus/projects/ccbfold/.venv/lib/python3.12/site-packages/Bio/PDB/StructureBuilder.py:100: PDBConstructionWarning: WARNING: Chain C is discontinuous at line 3124.
      warnings.warn(



    NGLWidget()


Other GUI options for protein structures include [pymol](https://biogrids.org/software/titles/pymol-open-source/) or [mol*](https://molstar.org/). The former is a commercial python package with a comprehensive software stack for creating molecule structures suitable for publications. The latter is a web-centric molecule viewer with integrations in the RCSB PDB, EMBL-EBI protein databases, and alphafold server. 

### Quality evaluation of Alphafold3 predicted structures using confidence scores

In order to hold practical significance, predicted protein structures must be accompanied by reliable measures of confidence. The output of AlphaFold3 returns confidence measures that can be used to identify regions of the predicted structure, and relative positions of different regions, that are more or less reliable. On a practial basis, the output of Alphafold3 contains the main metrics in the following files:

- `*_confidences.json`
    - the PAE and pLDDT
- `*_summary_confidences_.json`
    - the pTM and iPTM

#### PAE and pLDDT : Understanding predictive confidence 

These metrics broadly determine how the model trusts its atomic placement on a local(per atom/residue) and global(relative to domains) scale. These can help determine the plausibility of a structure by giving an estimate of how well the prediction will agree with an experimental structure. 


##### plDDT
pLDDT measures confidence in the local structure, estimating how well the prediction would agree with an experimental structure. It is based on the local distance difference test Cα (lDDT-Cα), which is a score that does not rely on superposition but assesses the correctness of the local distances (Mariani et al., 2013).

This per-atom confidence estimate is on a 0-100 scale where a higher value indicates higher confidence. pLDDT aims to predict a modified LDDT score that only considers distances to polymers. For proteins this is similar to the lDDT-Cα metric but with more granularity as it can vary per atom not just per residue. For ligand atoms, the modified LDDT considers the errors only between the ligand atom and polymers, not other ligand atoms. For DNA/RNA a wider radius of 30 Å is used for the modified LDDT instead of 15 Å.

The [mol* model viewer](https://molstar.org/) can color the 3D model by this scale allowing assessment of particular parts of the model.

![plddt scale](05_visualization/plddt_scale.png)


```python
from pathlib import Path 
from ipymolstar import PDBeMolstar
fpth= Path('/result/protein-dna-ion_pdb_7rce/protein-dna-ion_pdb_7rce_model.cif')
custom_data = {
    'data': fpth.read_bytes(),
    'format': 'cif',
    'binary': False,
    }
view = PDBeMolstar(
    custom_data=custom_data,
    alphafold_view=True
)
view
```

Another visualization method is to plot the confidence score for each atom. Elevated, stable values indicate well-defined domains, while dips highlight regions that might be disordered or flexible.


```python
#import functions and libraries that will aid in parsing and visualizing
import json
from eval_plots import plot_plddt
#load the json file containing the pLDDT and PAE data
F= open('/result/protein-dna-ion_pdb_7rce/protein-dna-ion_pdb_7rce_confidences.json')
data = json.load(F)
#plot the pLDDT values
plot_plddt(data['atom_plddts'], data['atom_chain_ids'])
```


    
![png](05_visualization/step1_visualize_af3_outputs_files/step1_visualize_af3_outputs_15_0.png)

##### PAE
PAE is defined as the expected positional error at residue X, measured in Ångströms (Å), if the predicted and actual structures were aligned on residue Y. Therefore, PAE is effectively a measure of how confident the model is that the domains are well packed and that the relative placement of the domains in the predicted structure is correct. Lower PAE values between residues belonging to different domains indicate well-defined relative positions and orientations—suggesting stable domain-domain relationships. Conversely, higher PAE values for such residue pairs suggest uncertainty in how these domains are arranged relative to each other, indicating that their positions and orientations should not be taken as definitive.

We can plot a heat map to quickly assesses how confidently the model predicts the relative positions and orientations of residues across the protein. A dashed line (`show_chains=True`) draws the boundary between different protein chains, or different types of molecules/atoms (ie, RNA/DNA/ligand).


```python
from eval_plots import plot_pae
#plot the PAE values
plot_pae(data['pae'], data['token_chain_ids'], show_chains=True)
```


    
![png](05_visualization/step1_visualize_af3_outputs_files/step1_visualize_af3_outputs_17_0.png)
    

---

You're now equipped with a variety of workflows for running AlphaFold3:

- A single end-to-end prediction of a protein-DNA complex
- Automatically generating inputs and running batch predictions
- Separately running MSA and inference to optimize compute usage
- Accelerating MSA using GPU-based mmseqs2

Feel free to extend these approaches to your own data and explore more complex modeling tasks. Happy folding!

```mermaid
flowchart TD
    A:::start@{label: "Start: Do you need to predict structures?"} --> B@{shape: "diamond", label: "How many structures?"}
    B --> C:::reqs@{shape: "rect", label: "#8226; 30 or fewer per day<br/>#8226; <5000 tokens per &nbsp;&nbsp;&nbsp;structure<br/>#8226; biologically common &nbsp;&nbsp;&nbsp;ligands / PTMs"} --> D:::server@{shape: "rect", label: "Use AlphaFold Server"}
    B --> E:::reqs@{shape: "rect", label: "#8226; A handful<br/>#8226; not matching server &nbsp;&nbsp;&nbsp;criteria"} --> F:::o2@{shape: "rect", label: "Use O2 AlphaFold3 module<br/>without modification"}
    B --> G:::reqs@{shape: "st-rect", label: "#8226; More than a handful<br/>#8226; exact AlphaFold3 spec &nbsp;&nbsp;&nbsp;desired"} --> H:::splitMsa@{shape: "rect", label: "Split MSA and Inference<br/>Run MSAs on CPU-only node"}
    B --> I:::reqs@{shape: "st-rect", label: "#8226; Many structures<br/>#8226; high throughput desired"} --> J:::mmseqs2@{shape: "rect", label: "Run MSAs with mmseqs2<br/>then use AlphaFold3 module for inference"}
    classDef reqs fill:#fff,stroke:#C9C9C9,text-align:left;
    classDef server fill:#DAE8FC,stroke:#6C8EBF;
    classDef o2 fill:#F8CECC,stroke:#B85450;
    classDef splitMsa fill:#FFF2CC,stroke:#D6B656;
    classDef mmseqs2 fill:#E1D5E7,stroke:#9673A6;
    classDef start fill: #fff, stroke:#C9C9C9; 
```
