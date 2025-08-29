# EF-G Mutations and Gentamicin Resistance in E. coli
*A reproducible Python workflow for structural analysis and mutation visualisation*

This repository contains part of my PhD research investigating the structural impact of mutations in the **elongation factor G (EF-G)** on bacterial susceptibility to gentamicin, using computational protein analysis.

## Background

A novel in vitro method was used to investigate the development of gentamicin resistance in Escherichia coli, Pseudomonas aeruginosa, and Klebsiella pneumoniae. This approach created a spatial gradient of increasing antibiotic concentrations, allowing bacteria to grow and adapt across zones of different drug exposure. As resistant subpopulations emerged, they could be isolated for further analysis.

DNA from gentamicin‑resistant colonies was sequenced, and the resulting data were aligned and annotated to detect genetic variants. Many resistant mutants carried substitutions in the fusA gene, which encodes elongation factor G (EF‑G) — a ribosomal GTPase essential for protein synthesis and a known target in antibiotic resistance studies.

🧪 Objectives
This project showcases a Python-based workflow for protein structure visualization and mutation analysis, developed to explore structural variations in EF-G and their potential relevance to antibiotic resistance. The emphasis is on:

- Generating EF-G variant structures using AlphaFold 3 via the AlphaFold server ([AlphaFold server link](https://alphafoldserver.com/)), based on curated FASTA sequences.

- Highlighting specific amino acid substitutions within the predicted models to assess their spatial positioning and potential structural implications

- Superimposing EF-G variants onto the 70S ribosome complexed with gentamicin to explore mutation proximity and generate functional hypotheses.

- Visualising and comparing wild-type and mutant structures using Biopython and NGL Viewer, with annotated mutation sites.

- Structuring the codebase for modularity and reproducibility, including Conda environment setup and clear function separation. 

- Documenting the workflow with markdown-based scientific reporting, combining code summaries, visual outputs, and interpretation  scaffolding in the domain IV of EF-G.



🧾 Disclaimer
This repository was created independently by Razan to demonstrate technical expertise in computational biology, Python-based structural analysis, and reproducible workflow design. While the underlying biological study originates from her PhD research, the coding workflow, visualisation pipeline, and documentation presented here were developed separately for portfolio and educational purposes. Scientific questions and interpretations were guided by the original supervisory team.


## Methods

🧬 1. Structural Analysis
- For EF-G structural context, crystal structures of the Thermus thermophilus 70S ribosome in complex with EF-G were retrieved from the Protein Data Bank:

- Pre-translocational state: `PDB ID: 4WPO`
- Post-translocational state: `PDB ID: 4V5F`

- EF-G chains were programmatically extracted using Biopython, and structural alignment was performed using MDAnalysis and NumPy to calculate RMSD (33.2 Å), reflecting large-scale conformational shifts during translocation.

- Mutation site Pro659 was highlighted within the aligned models using NGL Viewer, demonstrating its location in a dynamic region of EF-G potentially relevant to ribosomal movement.

- For comparative analysis, EF-G (post-translocational, 4V5F) was superimposed with the ribosome bound to gentamicin (PDB ID: 4V53), yielding an RMSD of 0.0005 Å, confirming structural consistency and enabling proximity-based visualisation of mutation effects.

🧪 2. Protein Modeling
- EF-G amino acid sequences were extracted from sequencing data and submitted to AlphaFold 3 via the AlphaFold server for structure prediction.

Predicted models showed confidence scores >88%, with most residues scoring >90 except in flexible loop regions.

- Wild-type and mutant EF-G structures were visualized and superimposed using Biopython and NGL Viewer within a Jupyter Notebook environment.

- RMSD between wild-type and mutant structures:
   0.52 Å for E. coli 36099
   0.46 Å for E. coli MG1655 
These minor structural shifts were localised to domains II and IV, and visualized in proximity to gentamicin within the ribosomal complex.

🧠 Findings
- Pro659 is located in EF-G Domain V, a region involved in ribosomal translocation. Structural alignment of pre- and post-translocational EF-G models revealed domain-level shifts, suggesting that the Pro659Leu substitution may influence conformational dynamics relevant to gentamicin susceptibility.

- Phe593Leu, located in Domain IV of EF-G, was visualised in close proximity to gentamicin within the ribosome-bound complex. Structural comparison between wild-type and mutant models showed a transition in Domain IV, supporting a potential role in resistance. This observation aligns with prior literature (Quiroga et al., 2018) identifying Phe593Leu as a resistance-associated mutation.

- The workflow provides a reproducible framework for visualising mutation positions, comparing structural variants, and generating hypothesis-driven insights into antimicrobial resistance mechanisms using Python-based tools.

## Repository Structure

efg_gentamicin_resistance/
├── data/                        # PDB files, sequencing data, and other raw inputs
├── models/                      # Superimposed and AlphaFold-generated structures
├── notebooks/                   # Jupyter Notebooks for analysis
├── scripts/                     # Python scripts for structural extraction and modelling
├── figures/                     # Plots and rendered structural images
├── efg_gentamicin_report.ipynb  # Summary notebook
├── environment.yml              # Conda environment for reproducibility
├── LICENSE                      # MIT License
└── README.md                    # Project overview

## Tools and Resources

- **Python** for protein extraction, superimposition, and RMSD calculation
- **Biopython** for sequence parsing and structure manipulation 
- **AlphaFold 3** for protein structure prediction ([AlphaFold server link](https://alphafoldserver.com/))  
- **NGL Viewer** for 3D visualization of protein structures
- **Jupyter Notebook** for workflow integration and scientific reporting 
- **PDB Database** for reference crystal structures

## Environment Setup

To ensure reproducibility, this project uses a Conda environment. You can recreate it using the provided `environment.yml` file:

```bash
conda env create -f environment.yml
conda activate ef-g-resistance
```
## Installation

To get started, clone the repository and install dependencies:

```bash
git clone https://github.com/RazanAb88/efg_gentamicin_resistance.git
cd efg_gentamicin_resistance
conda env create -f environment.yml
conda activate ef-g-resistance
```

## How to Run

1. Launch Jupyter Notebook
```bash
jupyter notebook
```
   Open the notebook from the notebooks/ folder (e.g., finding1_3d_alignment.ipynb).

2. Run the analysis
   Ensure the data/ folder contains the PDB or mmCIF files (4WPO, 4V5F, and AlphaFold models).
   Execute each cell in order to:
    -Extract EF-G from ribosome structures
    -Perform superimposition
    -Calculate RMSD
    -Visualise in 3D with NGL Viewer

3. View results
   Structural overlays are saved in the models/ folder.
   Plots and rendered images are stored in figures/.

