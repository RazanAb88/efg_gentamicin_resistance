<p align="center"> 
  <img src="figures/workflow2_gentamicin_proximity2.png" alt="Structural overlay of mutated E. coli EF-G Phe593Leu mutant" width="600"/>
</p>

<p align="center"><em>Structural overlay highlighting the Phe683Leu mutation in EF-G — visualising potential structural proximity to gentamicin binding site potentially linked to gentamicin resistance.</em></p>

# EF-G Mutations and Gentamicin Resistance in E. coli
*A reproducible Python workflow for structural analysis and mutation visualisation*

This repository contains part of my PhD research investigating the structural impact of mutations in the **elongation factor G (EF-G)** on bacterial susceptibility to gentamicin, using computational protein analysis.

## Background


A novel in vitro method was used to investigate the development of gentamicin resistance in Escherichia coli, Pseudomonas aeruginosa, and Klebsiella pneumoniae. This approach created a spatial gradient of increasing antibiotic concentrations, allowing bacteria to grow and adapt across zones of different drug exposure. As resistant subpopulations emerged, they could be isolated for further analysis.

DNA from gentamicin‑resistant colonies was sequenced, and the resulting data were aligned and annotated to detect genetic variants. Many resistant mutants carried substitutions in the fusA gene, which encodes elongation factor G (EF‑G) — a ribosomal GTPase essential for protein synthesis and a known target in antibiotic resistance studies.


## Objectives

This project showcases a Python-based workflow for protein structure visualization and mutation analysis:

- Generating EF-G variant structures using **AlphaFold 3**  
- Highlighting specific amino acid substitutions in predicted models  
- Superimposing EF-G variants onto the 70S ribosome with gentamicin  
- Visualising wild-type and mutant structures using **Biopython** and **NGL Viewer**  
- Structuring the codebase for modularity and reproducibility  


## Machine Learning Component

This workflow integrates AlphaFold 3, a deep-learning–based protein structure prediction model, as the machine learning component of the analysis. EF-G amino acid sequences derived from whole-genome sequencing data were submitted to the AlphaFold server, and the resulting structures were used as ML-derived inputs for downstream structural analysis.

Rather than treating AlphaFold outputs as ground truth, per-residue pLDDT confidence scores were explicitly extracted and analysed to guide interpretation. Confidence metrics were used to:

- Validate structural reliability at mutation sites
- Distinguish well-resolved domains from flexible or disordered regions
- Support cautious, evidence-based interpretation during structural superposition and visualisation

Predicted EF-G models consistently showed mean pLDDT scores >88, with lower confidence largely confined to flexible loop regions. This supports robust domain-level inference for mutation localisation while appropriately accounting for structural uncertainty.

AlphaFold-derived structures were subsequently integrated into downstream analyses, including:
- Structural superimposition with crystallographic ribosome–gentamicin complexes
- RMSD quantification between wild-type and mutant EF-G variants
- Spatial mapping of resistance-associated mutations relative to functional ribosomal sites
- Confidence-guided visualisation of mutation effects using Biopython and NGL Viewer

Together, this approach demonstrates the responsible integration of machine learning–based structure prediction into a reproducible computational biology workflow.


## Disclaimer 

This repository was created independently by Razan to demonstrate technical expertise in computational biology, Python-based structural analysis, and reproducible workflow design. While the underlying biological study originates from her PhD research, the coding workflow, visualisation pipeline, and documentation presented here were developed separately for portfolio and educational purposes. Scientific questions and interpretations were guided by the original supervisory team.


## Scientific Workflows

### Workflow 1: AlphaFold EF-G Modeling & Comparison

- EF-G amino acid sequences were submitted to AlphaFold 3 for structure prediction.  
- Predicted models had confidence scores >88%, with most residues >90%.  
- Wild-type and mutant EF-G structures were visualized and superimposed using Biopython and NGL Viewer.  
- **RMSD between wild-type and mutant structures**:  
  - 0.52 Å for *E. coli 36099*  
  - 0.46 Å for *E. coli MG1655*  


**Findings:**  
- **Pro659** is located in EF-G Domain V, potentially influencing ribosomal translocation, since visualising the structural alignment of pre- and post-translocational EF-G models revealed domain-level shifts, suggesting that the Pro659Leu substitution may influence conformational dynamics relevant to gentamicin susceptibility.

- **Phe593Leu** in Domain IV of  EF-G, and the visualising the structural comparison between wild-type and mutant models showed a transition in Domain IV, supporting a potential role in resistance. 

#### AlphaFold Confidence & Mutation Context 
AlphaFold 3–predicted EF-G models were analysed using per-residue pLDDT scores to assess structural confidence at mutation sites.

For E. coli 36099:

- Mean pLDDT: 91.78
- Mutation-site pLDDT: 92.93
- Mutation domain: Domain V

For E. coli MG1655:

- Mutation-site pLDDT: 96.22
- Mutation domain: Domain IV

Both mutations occur in regions of high predicted structural confidence, supporting robust localisation within functional EF-G domains. Domain-specific mutation positioning highlights potential strain-dependent structural mechanisms underlying gentamicin resistance.


### Workflow 2: EF-G + Ribosome Superimposition

- EF-G (post-translocational, 4V5F) was superimposed with the ribosome bound to gentamicin (PDB ID: 4V53).  
- **RMSD**: 0.7652 Å, confirming structural consistency.  
- Enables proximity-based visualisation of mutation effects within ribosome complexes.

**Findings:**
- Phe593Leu was visualised in close proximity to gentamicin within the ribosome-bound complex. This observation aligns with prior literature (Quiroga et al., 2018) identifying Phe593Leu as a resistance-associated mutation.

## Repository Structure

```text
efg_gentamicin_resistance/
├── data/                      # CIF/PDB files, sequencing data
├── models/                    # Superimposed and AlphaFold-generated structures
├── notebooks/                 # Workflow notebooks (Workflow 1, 2 and alphafold_features )
├── src/                       # Python scripts (bio_structures, visualise_structures, ribosome_drug_proximity & alphafold_features)
├── figures/                   # Plots and structural images
├── efg_gentamicin_report.ipynb # Summary notebook
├── environment.yml            # Conda environment
├── conda                      # conda file
├── requirements.txt           # Pip dependencies (optional)
├── LICENSE                    # MIT License
└── README.md                  # Project overview
```

## Tools and Resources

- Python – protein extraction, superimposition, RMSD calculation

- Biopython – sequence parsing & structure manipulation

- AlphaFold 3 – protein structure prediction (AlphaFold server: https://alphafold.ebi.ac.uk/)

- NGL Viewer – 3D protein visualization

- Jupyter Notebook – workflow integration and reporting

- PDB Database – reference crystal structures

## Environment Setup

Recreate the Conda environment:

```bash
conda env create -f environment.yml
conda activate ef-g-resistance
```
Optional Pip installation:

```bash 
pip install -r requirements.txt
```

## How to Run

1- Launch Jupyter Notebook:
``` bash
jupyter notebook
```
2- Open notebooks in notebooks/ folder.

3- Ensure data/ contains CIF or PDB files.

4- Execute the cells to:

- Extract EF-G from ribosome structures

- Perform superimposition

- Calculate RMSD

- Visualize in 3D

5- Results:

- Superimposed structures → models/

- Figures & visualizations → figures/

### Skills Demonstrated

- Protein structure prediction with AlphaFold 3

- Structural alignment and RMSD calculation using Biopython

- Mutation visualization with NGL Viewer

- Reproducible workflow design with Conda and modular scripts

- Notebook-based scientific reporting



