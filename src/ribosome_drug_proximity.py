#!/usr/bin/env python
# coding: utf-8

"""
Workflow 2: Ribosome + EF-G + S12 alignment + visualization

Superimpose EF-G (from 4V5F) onto the 4V53 ribosome using S12 chain alignment,
then embed EF-G into the 4V53 ribosome, save the structure, and optionally visualize it.
"""

import copy
from Bio.PDB import MMCIFParser, MMCIFIO, Superimposer

from bio_structures import save_model, rename_chains, extract_chains, parse_structure
from visualise_structures import (
    visualize_individual_model,
    visualize_superimposed_models,
    load_and_display_cif
)

def extract_ca_atoms(model, chain_id):
    """Return a list of C-alpha atoms from a specific chain."""
    return [res["CA"] for res in model[chain_id] if "CA" in res]

def run_workflow2(
    cif_4v5f="4v5f.cif",
    cif_4v53="4v53.cif",
    chain_efg="CY",
    chain_s12_4v5f="CL",
    chain_s12_4v53="CL",
    output_file="4v53_with_EFG_superimposed.cif",
    visualize=False
):
    """Run Workflow 2: Ribosome alignment + EF-G placement, optionally visualize."""
    parser = MMCIFParser(QUIET=True)

    # --- Load structures ---
    s_4v5f = parser.get_structure("4v5f", cif_4v5f)
    s_4v53 = parser.get_structure("4v53", cif_4v53)

    # --- Extract chains ---
    efg_s12_model = extract_chains(copy.deepcopy(s_4v5f[0]), [chain_efg, chain_s12_4v5f])
    s12_4v53_model = extract_chains(copy.deepcopy(s_4v53[0]), [chain_s12_4v53])

    # --- Prepare CA atoms for superposition ---
    atoms_4v5f = extract_ca_atoms(efg_s12_model, chain_s12_4v5f)
    atoms_4v53 = extract_ca_atoms(s12_4v53_model, chain_s12_4v53)

    # Equalize number of residues
    L = min(len(atoms_4v5f), len(atoms_4v53))
    atoms_4v5f, atoms_4v53 = atoms_4v5f[:L], atoms_4v53[:L]

    print(f"Aligning S12 C-alpha atoms: {L} residues")

    # --- Superimpose ---
    sup = Superimposer()
    sup.set_atoms(atoms_4v53, atoms_4v5f)
    sup.apply(efg_s12_model.get_atoms())
    print(f"RMSD: {sup.rms:.4f} Å")

    # --- Merge EF-G into full ribosome ---
    merged_ribosome = copy.deepcopy(s_4v53[0])
    merged_ribosome.add(efg_s12_model[chain_efg])

    # --- Save final structure ---
    io = MMCIFIO()
    io.set_structure(merged_ribosome)
    io.save(output_file)
    print(f"Saved file: {output_file}")

    # --- Optional visualization ---
    if visualize:
        # Load and rename chains to avoid duplicates
        structure = parse_structure(output_file, "aligned")
        model = structure[0]
        seen = {}
        for chain in list(model):
            cid = chain.id
            if cid in seen:
                seen[cid] += 1
                chain.id = f"{cid}_{seen[cid]}"
            else:
                seen[cid] = 1
        renamed_file = "alignedS12+EF-G_renamed.cif"
        save_model(structure, renamed_file)

        # Visualize
        view = load_and_display_cif(renamed_file, color="lightgreen")

        # Highlight chains
        chain_colors = {chain_efg: "hotpink", chain_s12_4v5f: "blue", f"{chain_s12_4v5f}_2": "silver"}
        for chain, color in chain_colors.items():
            view.add_representation("cartoon", selection=f":{chain}", color=color)

        # Highlight PHE:582 in EF-G chain
        phe_selection = f"582:{chain_efg}"
        view.add_representation("ball+stick", selection=phe_selection, color="red")
        view.add_representation(
            "label",
            sele=phe_selection,
            labelType="format",
            labelFormat="%(resname)s%(resno)s",
            color="red",
            xOffset=1,
            fixedSize=True
        )
        print(f"Added PHE:582 in chain {chain_efg} as red stick")

        # Highlight gentamicin (LLL)
        gentamicin_selections = ["LLL", "[LLL]", "hetero and not water", "ligand"]
        for selection in gentamicin_selections:
            try:
                view.add_representation("ball+stick", selection=selection, color="red", radius=0.5)
                print(f"Added gentamicin representation with selection: {selection}")
                break
            except Exception:
                continue
        view.add_label(text="Gentamicin", selection="LLL", color="black", size=1)

        from IPython.display import display
        display(view)


if __name__ == "__main__":
    # Example usage
    run_workflow2(
        cif_4v5f="data/4v5f.cif",
        cif_4v53="data/4v53.cif",
        output_file="models/4v53_with_EFG_superimposed.cif",
        visualize=True  # Set True to display the structure
    )
