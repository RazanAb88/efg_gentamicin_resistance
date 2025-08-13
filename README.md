# EF-G Mutations and Gentamicin Resistance in E. coli

This repository contains part of my PhD research investigating the structural impact of mutations in the **elongation factor G (EF-G)** on bacterial susceptibility to gentamicin, using computational protein analysis.

## Background

Antibiotic Resistance Growth Plates (ARGP) were used to study genotypic mutations underlying antibiotic resistance in *E. coli* and *P. aeruginosa*. The ARGP provides a spatial gradient of increasing antibiotic concentrations to investigate resistance development in strains previously reported to be motile in drug-free media.  

DNA from gentamicin-resistant bacteria was sequenced at the **MicrobesNG facility, University of Birmingham**. Sequencing data were aligned and annotated to detect variants. Many gentamicin-resistant mutants carried mutations in the **fusA** gene, which encodes **EF-G**.

## Objective

This project focuses on understanding the potential role of EF-G mutations on antibiotic susceptibility. Specifically, we studied a clinical strain of *E. coli* (**E. coli 36099**) and its mutants resistant to 10×MIC of gentamicin, which carried a **Pro659Leu mutation**.

## Methods

### 1. Structural Analysis
- Crystal structures of **Thermus thermophilus 70S ribosome in complex with EF-G** were used:
  - Pre-translocational state: **PDB ID: 4WPO**
  - Post-translocational state: **PDB ID: 4V5F**
- EF-G was extracted from these ribosome models using Python, superimposed, and RMSD calculated (**33.2 Å**).  
  > Note: The large RMSD reflects conformational changes between the pre- and post-translocational states.
- Key observation: **Pro659** is located in a region of EF-G critical for movement during the transition between states. The **Pro659Leu mutation** may affect this motion.

### 2. Protein Modeling
- The FASTA sequence of EF-G was extracted from sequencing data and modeled using **AlphaFold 3**.
- Models for both the parent strain and 10×MIC-resistant mutant were visualized in **Jupyter Notebook** and superimposed.
- RMSD between wild type and mutant: **0.52 Å**, indicating minor structural alterations in **domains II and IV**, which may impact gentamicin binding in the 70S ribosome.

## Findings

- Pro659 is positioned in a functionally important region of EF-G.
- Structural comparison between wild type and mutant EF-G reveals subtle domain changes that could influence gentamicin susceptibility.
- Provides a structural framework for understanding EF-G–mediated antimicrobial resistance.

## Repository Structure

efg_gentamicin_resistance/
├── data/ # PDB files, sequencing data, and other raw inputs
├── models/ # Superimposed and AlphaFold-generated structures
├── notebooks/ # Jupyter Notebooks for analysis
├── scripts/ # Python scripts for structural extraction and modeling
├── LICENSE # MIT License
└── README.md # Project overview

## Tools and Resources

- **Python** for protein extraction, superimposition, and RMSD calculation
- **AlphaFold 3** for protein structure prediction
- **NGL Viewer** for 3D visualization of protein structures
- **PDB Database** for reference crystal structures

## How to Run

1. **Clone the repository**
   ```bash
   git clone https://github.com/RazanAb88/efg_gentamicin_resistance.git
   cd efg_gentamicin_resistance
2. **Install dependencies**
   Make sure you have Python 3.9+ and install required packages:
    pip install -r requirements.txt

3. Launch Jupyter Notebook
   jupyter notebook
   Open the notebook from the notebooks/ folder (e.g., finding1_3d_alignment.ipynb).

4. Run the analysis
   Ensure the data/ folder contains the PDB or mmCIF files (4WPO, 4V5F, and AlphaFold models).
- Execute each cell in order to:

    -Extract EF-G from ribosome structures
    -Perform superimposition
    -Calculate RMSD
    -Visualise in 3D with NGL Viewer

5. View results

   Structural overlays are saved in the models/ folder.

   Plots and rendered images are stored in figures/.
