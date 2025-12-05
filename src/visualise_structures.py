#!/usr/bin/env python
# coding: utf-8


# visualize_structures.py

import nglview as nv
from IPython.display import display



def load_and_display_cif(file_path, color='pink', highlight=None):
    """
    Load and display a CIF with a single color and optional residue highlight.
    """
    view = nv.show_file(file_path)
    view.clear_representations()
    view.add_representation("cartoon", color=color)

    if highlight:
        selection = f"resi {highlight}"           
        view.add_representation("ball+stick", selection=highlight, color="red")
        view.add_representation(
            "label",
            sele=highlight,
            labelType="format",
            labelFormat="%(resname)s%(resno)s",
            color="red",
            fixedSize=True
        )
    return view


def visualize_individual_model(file_path, color='pink', highlight=None):
    """
    Visualize a single CIF/PDB model with optional residue highlight.
    """
    view = nv.show_file(file_path)
    view.clear_representations()
    view.add_representation("cartoon", color=color)

    if highlight:
        # Highlight specific residue(s)
        view.add_representation("ball+stick", selection=highlight, color="red")
        view.add_representation(
            "label",
            sele=highlight,
            labelType="format",
            labelFormat="%(resname)s%(resno)s",
            color="red",
            fixedSize=True
        )
    return view



def visualize_superimposed_models(file_list, colors):
    """
    Display multiple aligned CIF/PDB files as separate components.
    file_list: list of file paths
    colors: list of colors corresponding to each component
    """
    view = nv.NGLWidget()

    for i, (file, color) in enumerate(zip(file_list, colors)):
        comp = view.add_component(file)
        view.clear_representations(component=i)
        view.add_representation("cartoon", color=color, component=i)

    return view


