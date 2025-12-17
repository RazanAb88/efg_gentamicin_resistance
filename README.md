<p align="center"> 
  <img src="figures/workflow2_gentamicin_proximity2.png" alt="Structural overlay of mutated E. coli EF-G Phe593Leu mutant" width="600"/>
</p>

<p align="center"><em>Structural overlay highlighting the Phe683Leu mutation in EF-G — visualising potential structural proximity to gentamicin binding site potentially linked to gentamicin resistance.</em></p>

# EF-G Mutations and Gentamicin Resistance in E. coli
*A reproducible Python workflow for structural analysis and mutation visualisation*

This repository contains part of my PhD research investigating the structural impact of mutations in the **elongation factor G (EF-G)** on bacterial susceptibility to gentamicin, using computational protein analysis.

➡️ Full reproducible structural analysis is available in the notebooks/ directory, including RMSD analysis, PCA-based clustering, and supervised feature interpretation.


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
- Quantifying mutation-specific structural perturbation using backbone and side-chain RMSD  
- Applying supervised and unsupervised learning to interpret structural variation  


## Machine Learning Component

Machine learning is incorporated into this project in two complementary ways:

### 1. AlphaFold 3 Structural Prediction (Deep Learning Component)
EF-G amino acid sequences were submitted to the AlphaFold 3 server, and the predicted structures were used as ML-derived inputs for downstream analysis. Per-residue pLDDT confidence scores were extracted to assess structural reliability at mutation sites, ensuring responsible use of predicted models in mechanistic interpretation.

AlphaFold predictions showed consistently high confidence (mean pLDDT >88), with mutation-site confidence >92. These scores support robust domain-level inference while 
appropriately highlighting flexible loop regions that require cautious interpretation.

### 2. Interpretable Learning on Structural Features (Supervised & Unsupervised)
Two lightweight machine-learning analyses were performed:

- **Supervised regression (A1-enhanced):**  
  Used to interpret which structural features (position, domain, side-chain size, backbone RMSD) contribute to side-chain displacement.  
  *Due to n=2 samples, the model is used strictly for interpretability, not prediction.*

- **Unsupervised PCA clustering:**  
  Used to explore whether mutations occupy distinct structural regimes based on   global RMSD, domain and positional context.

Together, these analyses demonstrate *responsible, interpretable integration* of both deep-learning–derived structural models and classical ML methods within a reproducible 
computational biology workflow.



## Disclaimer

This repository was created independently by Razan as a demonstration of technical skills in computational biology, protein structure analysis, and reproducible pipeline 
design. Although the mutations originate from her PhD research, all code, workflows, visualisation methods, and documentation presented here were developed separately for 
portfolio and educational purposes. Scientific interpretation was guided by the original supervisory team, but the computational implementation and presentation are 
the author's own.



## Scientific Workflows

### Workflow 1: AlphaFold EF-G Modeling & Comparison

- EF-G amino acid sequences were submitted to AlphaFold 3 for structure prediction.  
- Predicted models showed high confidence (mean pLDDT >88%; most residues >90%).  
- Wild-type and mutant EF-G structures were visualised and superimposed using Biopython and NGL Viewer.  
- **Global RMSD between wild-type and mutant EF-G structures**:  
  - 0.52 Å for *E. coli 36099*  
  - 0.46 Å for *E. coli MG1655*

**Findings:**  
- **Pro659** is located in EF-G Domain V. Structural alignment of pre- and post-translocational EF-G models revealed domain-level shifts in this region, suggesting that the Pro659Leu mutation may influence conformational dynamics relevant to gentamicin susceptibility.  
- **Phe593Leu**, located in Domain IV, exhibited a local structural shift relative to the wild-type model. This domain plays a functional role in ribosome–antibiotic interactions, consistent with its reported association with aminoglycoside resistance.

#### AlphaFold Confidence & Mutation Context 

Per-residue pLDDT scores were inspected to assess structural reliability at the mutation sites:

**E. coli 36099 (Pro659Leu):**
- Mean pLDDT: 91.78  
- Mutation-site pLDDT: 92.93  
- Domain: V  

**E. coli MG1655 (Phe593Leu):**
- Mutation-site pLDDT: 96.22  
- Domain: IV  

Both mutations occur in regions of high predicted confidence, supporting robust localisation within EF-G’s functional domains. The domain context also highlights potential strain-dependent structural mechanisms contributing to gentamicin resistance.

---

### Workflow 2: EF-G + Ribosome Superimposition

- EF-G (post-translocational; PDB 4V5F) was superimposed onto the ribosome–gentamicin complex (PDB 4V53).  
- **RMSD of superposition**: 0.7652 Å, indicating strong geometric compatibility.  
- This workflow enables proximity-based visualisation of EF-G mutations within a biologically relevant ribosome–antibiotic context.

**Findings:**  
- Phe593Leu was visualised in close proximity to gentamicin within the ribosome-bound complex.  
- This aligns with prior literature (Quiroga et al., 2018) identifying Phe593Leu as a resistance-associated substitution.

---

### Workflow 3: Local Structural Perturbation Analysis (Backbone + Side-Chain RMSD)

This workflow quantifies mutation-specific structural effects by computing **side-chain RMSD** and **local backbone RMSD** after Cα-based superposition.

#### Methods
For each mutation:
- Wild-type and mutant structures were superimposed using common Cα atoms.  
- Side-chain RMSD was computed for all non-backbone atoms at the mutation site.  
- Backbone RMSD was computed from N, Cα, and C atoms at the same residue.  
- This produces a structural signature per mutation:  
  `{ position, domain, strain, n_sc_atoms, sidechain_rmsd, backbone_rmsd }`

#### Findings
- **Phe593Leu** shows **large side-chain displacement** (~1.77 Å) and moderate backbone adjustment.  
- **Pro659Leu** shows **moderate side-chain displacement** (~0.74 Å) and moderate backbone shift, consistent with the loss of Proline’s conformational rigidity.

These findings reveal **distinct local perturbation regimes** between Domain IV and Domain V mutants and provide mechanistic descriptors for downstream ML analyses.

---

### Workflow 4: Supervised Feature-Based Interpretation (A1-Enhanced)

A supervised regression model was used to explore which structural features best explain mutation-induced side-chain displacement. Predictor variables included:

- Residue position  
- EF-G domain  
- Strain background  
- Side-chain atom count (n_sc_atoms)  
- Local backbone RMSD  

Side-chain RMSD was used as the target variable because it directly reflects the local chemical and geometric impact of each mutation.

#### Purpose
With an intentionally small dataset (n = 2), the model is used **exclusively for interpretability**.  
LOOCV results are shown for completeness but are **not** interpreted as predictive metrics.

#### Key Insight
The mutations occupy **distinct structural regimes**:
- Phe593Leu → high-perturbation  
- Pro659Leu → moderate-perturbation  

Backbone RMSD emerges as an indicator of **local geometric accommodation**, justifying its inclusion as a mechanistic feature.

---

### Workflow 5: Unsupervised Structural Clustering

An unsupervised analysis was performed to determine whether the two EF-G variants differ in structural context. Principal Component Analysis (PCA) was applied to:

- Residue position  
- Domain  
- Strain  
- Global RMSD  

#### Findings
Despite the minimal dataset, PCA cleanly separates the two mutations, driven primarily by their **domain context and mutation position**. This workflow establishes a scalable framework for incorporating additional mutants as they become available.


## Repository Structure

```text
efg_gentamicin_resistance/
├── data/                      # CIF/PDB files, sequencing data
├── models/                    # Superimposed and AlphaFold-generated structures
├── notebooks/                 # Workflow notebooks                
│   ├── workflow1_alphafold.ipynb
│   ├── workflow2_ribosome_superimposition.ipynb
│   ├── workflow3_local_rmsd.ipynb
│   ├── unsupervised_structural_clustering.ipynb
│   ├── supervised_feature_analysis.ipynb
│   └──alphafold_features.ipynb
│
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

- **Python** – core scripting for structure handling, RMSD measurement, and automation  
- **Biopython** – structure parsing, chain extraction, and coordinate manipulation  
- **AlphaFold 3 (EBI server)** – deep-learning–based protein structure prediction  
- **NGL Viewer** – interactive 3D molecular visualisation  
- **PDB / RCSB Database** – reference ribosome and EF-G crystal structures  
- **NumPy / Pandas / Matplotlib** – scientific computing and plotting  
- **scikit-learn** – supervised and unsupervised ML analysis  
- **Jupyter Notebook** – integrated workflow development and scientific reporting


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
2- Open notebooks from the notebooks/ directory.

3- Ensure that the data/ folder contains the required CIF or PDB files (wild-type, mutant, and ribosome structures).

4- Execute workflow notebooks to:

- Extract EF-G from ribosomal structures

- Perform structural superimpositions

- Compute global, backbone, and side-chain RMSDs

- Visualise mutations in 3D (NGL Viewer)

- Run supervised and unsupervised analyses

5- Results:

- Superimposed models → models/

- Plots and figures → figures/

- Processed tables → data/processed/ (if exported)

## Skills Demonstrated

- Structural bioinformatics using Biopython  
- Protein structure prediction (AlphaFold 3)  
- RMSD computation (global, backbone, and side-chain)  
- Structural superimposition of EF-G and ribosome complexes  
- Mutation visualisation using NGL Viewer  
- PCA-based structural clustering  
- Interpretable supervised regression for mutation feature analysis  
- Modular Python workflow design (src/ + notebooks/)  
- Reproducible computational environment with Conda  
- Scientific documentation and figure generation in Jupyter



