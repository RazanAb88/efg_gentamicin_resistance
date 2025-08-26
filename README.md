# EF-G Mutations and Gentamicin Resistance in E. coli

This repository contains part of my PhD research investigating the structural impact of mutations in the **elongation factor G (EF-G)** on bacterial susceptibility to gentamicin, using computational protein analysis.

## Background

Antibiotic Resistance Growth Plates (ARGP) were used to investigate antibiotic resistance development  in *E. coli* and *P. aeruginosa* *K. pneumoniae*. The ARGP provides a spatial gradient of increasing antibiotic concentrations that initiate mutants from strains previously reported to be motile in drug-free media.  

DNA from gentamicin-resistant bacteria was sequenced at the **MicrobesNG facility, University of Birmingham**. Sequencing data were aligned and annotated to detect variants. Many gentamicin-resistant mutants carried mutations in the **fusA** gene, which encodes **EF-G**.

## Objective
This project focuses on understanding the potential role of EF-G mutations on antibiotic susceptibility. Specifically, we studied:

- A clinical strain of E. coli (E. coli 36099) and its mutant resistant to 10×MIC of gentamicin, which carried a Pro659Leu mutation in the domain V of EF-G.

- A reference strain of E. coli (E. coli MG1655) and its mutant resistant to 10×MIC of gentamicin, which had a Phe593Leu mutation in the domain IV of EF-G.

## Methods

### 1. Structural Analysis
-  For E. coli 36099: Crystal structures of **Thermus thermophilus 70S ribosome in complex with EF-G** were used:
  - Pre-translocational state: **PDB ID: 4WPO**
  - Post-translocational state: **PDB ID: 4V5F**
- EF-G was extracted from these ribosome models using Python, superimposed, and RMSD calculated (**33.2 Å**). This high RMSD reflects large-scale conformational shifts inherent to EF-G’s function during translocation, not modeling error.
- Key observation: **Pro659** is located in a region of EF-G critical for movement during the transition between states. The **Pro659Leu mutation** may affect this motion.

- For E. coli MG1655: EF-G in the post-translocational state (PDB ID: 4V5F) was superimposed with the bacterial ribosome in complex with gentamicin (PDB ID: 4V53).
  - RMSD = 0.0005 Å, indicating an almost perfect superimposition and confirming structural consistency between the ribosome-bound gentamicin and EF-G.
### 2. Protein Modeling
- The FASTA sequence of EF-G was extracted from sequencing data and modelled using **AlphaFold 3** with confidence > 88% for the models. Most residues scored >90, with lower confidence in flexible loop regions.
- Models for both the parent strain and 10×MIC-resistant mutant were visualised in **Jupyter Notebook** and superimposed.
- RMSD between wild type and mutant: **0.52 Å** for E. coli 36099 & **0.46 Å** for E. coli MG1655 , indicating minor structural alterations in **domains II and IV**, which may impact gentamicin binding in the 70S ribosome.

## Findings

- Pro659 is positioned in the domain V in a functionally important region of EF-G. Structural comparison reveals domain changes that could influence gentamicin susceptibility.

- Phe593Leu in EF-G Domain IV is located close to the gentamicin binding site. This mutation may interfere with gentamicin binding. A transition in Domain IV was observed between parent and mutant models. This finding is consistent with Quiroga et al. (2018), who reported Phe593Leu as a causative factor for gentamicin resistance in E. coli MG1655.

- Together, the studies provide a structural framework for understanding EF-G–mediated antimicrobial resistance

## Repository Structure

efg_gentamicin_resistance/
├── data/ # PDB files, sequencing data, and other raw inputs
├── models/ # Superimposed and AlphaFold-generated structures
├── notebooks/ # Jupyter Notebooks for analysis
├── scripts/ # Python scripts for structural extraction and modelling
├──figures/ # Plots and rendered structural images
├── efg_gentamicin_report.ipynb
├── environment.yml
├── LICENSE # MIT License
└── README.md # Project overview

## Tools and Resources

- **Python** for protein extraction, superimposition, and RMSD calculation
- **AlphaFold 3** for protein structure prediction ([AlphaFold server link](https://alphafoldserver.com/))  
- **NGL Viewer** for 3D visualization of protein structures
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
   jupyter notebook
   Open the notebook from the notebooks/ folder (e.g., finding1_3d_alignment.ipynb).

2. Run the analysis
   Ensure the data/ folder contains the PDB or mmCIF files (4WPO, 4V5F, and AlphaFold models).
- Execute each cell in order to:

    -Extract EF-G from ribosome structures
    -Perform superimposition
    -Calculate RMSD
    -Visualise in 3D with NGL Viewer

3. View results

   Structural overlays are saved in the models/ folder.

   Plots and rendered images are stored in figures/.

