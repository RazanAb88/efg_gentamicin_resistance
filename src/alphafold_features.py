from Bio.PDB import MMCIFParser
import numpy as np

# -------------------------------
# EF-G domain map (E. coli)
# -------------------------------
DOMAIN_MAP = {
    (1, 300): "Domain I",
    (301, 399): "Domain II",
    (400, 485): "Domain III",
    (486, 609): "Domain IV",
    (610, 698): "Domain V"
}

def residue_to_domain(residue_number):
    """Map residue number to EF-G domain."""
    for (start, end), domain in DOMAIN_MAP.items():
        if start <= residue_number <= end:
            return domain
    return "Unknown"


def parse_alphafold_structure(cif_file, structure_id="af"):
    parser = MMCIFParser(QUIET=True)
    return parser.get_structure(structure_id, cif_file)


def extract_plddt_scores(structure, chain_id="A"):
    """
    Extract per-residue pLDDT scores.
    AlphaFold stores pLDDT in B-factor.
    """
    scores = {}
    model = structure[0]
    for res in model[chain_id]:
        if "CA" in res:
            resnum = res.id[1]
            scores[resnum] = res["CA"].bfactor
    return scores


def summarize_plddt(plddt_scores):
    return {
        "mean_pLDDT": float(np.mean(list(plddt_scores.values()))),
        "min_pLDDT": float(np.min(list(plddt_scores.values()))),
        "max_pLDDT": float(np.max(list(plddt_scores.values())))
    }


def mutation_plddt(plddt_scores, residue_number):
    return plddt_scores.get(residue_number, None)


def mutation_domain(residue_number):
    return residue_to_domain(residue_number)
