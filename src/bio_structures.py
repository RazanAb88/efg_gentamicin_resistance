#!/usr/bin/env python
# coding: utf-8

# === Imports ===
from Bio.PDB import (
    MMCIFParser, MMCIFIO, Superimposer
)


# --- Parsing & Saving --- #

def parse_structure(file_path, structure_id):
    """Parse a CIF file and return structure object."""
    parser = MMCIFParser(QUIET=True)
    return parser.get_structure(structure_id, file_path)


def rename_chains(structure, new_chain_id):
    """Rename all chains in a structure (AlphaFold → chain A only)."""
    for model in structure:
        for chain in model:
            chain.id = new_chain_id
    return structure


# --- AlphaFold Superimposition --- #

def superimpose_structures(file1, file2):
    """
    Superimpose WT and mutant EF-G (AlphaFold single-chain models).
    WT  → Chain A
    Mut → Chain B
    Returns RMSD, aligned model1, aligned model2.
    """
    s1 = parse_structure(file1, "wt")
    s2 = parse_structure(file2, "mut")

    rename_chains(s1, "A")
    rename_chains(s2, "B")

    m1 = s1[0]
    m2 = s2[0]

    atoms1 = [res["CA"] for res in m1["A"] if "CA" in res]
    atoms2 = [res["CA"] for res in m2["B"] if "CA" in res]

    L = min(len(atoms1), len(atoms2))
    atoms1, atoms2 = atoms1[:L], atoms2[:L]

    sup = Superimposer()
    sup.set_atoms(atoms1, atoms2)
    sup.apply(m2.get_atoms())

    return sup.rms, m1, m2


# --- General Model Manipulation --- #
import copy 
def list_chain_ids(model):
    return [chain.id for chain in model.get_chains()]


def rename_chain(model, old_id, new_id):
    for chain in model.get_chains():
        if chain.id == old_id:
            chain.id = new_id
    return model


def save_model(structure, output_path):
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(output_path)



def extract_chains(model, chain_ids):
    """
    Return a copy of the model containing only the specified chains.
    
    Parameters:
        model: Bio.PDB.Model
        chain_ids: list of chain IDs to keep
    
    Returns:
        new_model: Bio.PDB.Model containing only requested chains
    """
    new_model = copy.deepcopy(model)
    for chain in list(new_model):
        if chain.id not in chain_ids:
            new_model.detach_child(chain.id)
    return new_model


