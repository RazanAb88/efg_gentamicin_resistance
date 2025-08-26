#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Structure parsing and manipulation
from Bio.PDB import (
    MMCIFParser, MMCIFIO, Select, PDBList, Superimposer,
    Model, Structure, FastMMCIFParser, PPBuilder
)
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

# Sequence handling
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

# Visualization
import matplotlib.pyplot as plt
import nglview as nv
import ipywidgets as widgets
from IPython.display import display


# In[ ]:


parser = MMCIFParser()
structure = parser.get_structure("4WPO", "4wpo.cif") 


# In[ ]:


mmcif_dict = MMCIF2Dict("4wpo.cif")

# This will list all chain IDs and associated entity IDs
strand_ids = mmcif_dict.get("_entity_poly.pdbx_strand_id", [])
entity_ids = mmcif_dict.get("_entity_poly.entity_id", [])
polymer_names = mmcif_dict.get("_entity.pdbx_description", [])

for entity, strands, name in zip(entity_ids, strand_ids, polymer_names):
    print(f"Entity ID: {entity}, Chains: {strands}, Molecule: {name}")


# In[ ]:


# Load the structure
parser = MMCIFParser()
structure = parser.get_structure("4WPO", "4wpo.cif")

# Custom Select class to keep only EF-G chains
class EF_G_ChainSelect(Select):
    def __init__(self, chain_ids):
        self.chain_ids = chain_ids

    def accept_chain(self, chain):
        return chain.id in self.chain_ids

# Set up the writer
io = MMCIFIO()
io.set_structure(structure)

# Save only EF-G chains (GB and ND)
io.save("4WPO_EF-G.cif", EF_G_ChainSelect(["BZ"]))


# In[ ]:


widgets.IntSlider()
view = nv.show_file("4WPO_EF-G.cif")
slider = widgets.IntSlider(min=0, max=100, step=1, value=50, description='My Slider')
view.clear_representations()
view.add_representation('cartoon', color='pink')  
view.add_representation('stick', selection='540 and PRO', color='black')
display(view)
display(slider)


# In[ ]:


# Parse the structure
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("4WPO", "4WPO_EF-G.cif")

ppb = PPBuilder()
records = []

# Extract sequences for each chain
for model in structure:
    for chain in model:
        seq = ""
        for pp in ppb.build_peptides(chain):
            seq += str(pp.get_sequence())

        if seq:
            record = SeqRecord(Seq(seq), id=f"Chain_{chain.id}", description="")
            records.append(record)

# Print the sequences in FASTA format
for record in records:
    print(f">{record.id}\n{record.seq}")


# In[ ]:


parser = MMCIFParser()
structure = parser.get_structure("4v5f", "4v5f.cif") 


# In[ ]:


mmcif_dict = MMCIF2Dict("4v5f.cif")

# This will list all chain IDs and associated entity IDs
strand_ids = mmcif_dict.get("_entity_poly.pdbx_strand_id", [])
entity_ids = mmcif_dict.get("_entity_poly.entity_id", [])
polymer_names = mmcif_dict.get("_entity.pdbx_description", [])

for entity, strands, name in zip(entity_ids, strand_ids, polymer_names):
    print(f"Entity ID: {entity}, Chains: {strands}, Molecule: {name}")


# In[ ]:


# Load the structure
parser = MMCIFParser()
structure = parser.get_structure("4v5f", "4v5f.cif")

# Custom Select class to keep only EF-G chains
class EF_G_ChainSelect(Select):
    def __init__(self, chain_ids):
        self.chain_ids = chain_ids

    def accept_chain(self, chain):
        return chain.id in self.chain_ids

# Set up the writer
io = MMCIFIO()
io.set_structure(structure)

# Save only EF-G chains 
io.save("4v5f_EF-G.cif", EF_G_ChainSelect(["CY"]))


# In[ ]:


## for finding proline number in 4v5f structure:
# Parse the structure
parser = MMCIFParser()
structure = parser.get_structure('4V5F', '4v5f_EF-G.cif')

# Iterate over models, chains, and residues to list amino acids and their residue numbers
for model in structure:
    for chain in model:
        print(f"Chain ID: {chain.id}")
        for residue in chain:
            # residue.id is a tuple like (' ', resseq, insertion_code)
            resseq = residue.id[1]  # residue sequence number
            resname = residue.get_resname()  # 3-letter amino acid code
            print(f"Residue {resseq}: {resname}")
        print()  # Blank line between chains


# In[ ]:


view = nv.show_file("4v5f_EF-G.cif")
slider = widgets.IntSlider(min=0, max=100, step=1, value=50, description='My Slider')
view.clear_representations()
view.add_representation('cartoon', color='silver')
selection_string = "648:CY"
view.add_representation('ball+stick', selection=selection_string)
view.add_representation('label', sele=selection_string, labelType='format', labelFormat='%(resname)s%(resno)s', color='red', xOffset=1, fixedSize=True)




display(view)
display(slider)


# In[ ]:


# Parse structures
parser = MMCIFParser()
structure1 = parser.get_structure('s1', '4wpo_EF-G.cif')
structure2 = parser.get_structure('s2', '4v5f_EF-G.cif')

model1 = structure1[0]
model2 = structure2[0]

chain1_id = 'BZ'
chain2_id = 'CY'

atoms1 = [res['CA'] for res in model1[chain1_id] if 'CA' in res]
atoms2 = [res['CA'] for res in model2[chain2_id] if 'CA' in res]

min_len = min(len(atoms1), len(atoms2))
atoms1 = atoms1[:min_len]
atoms2 = atoms2[:min_len]

# Superimpose
sup = Superimposer()
sup.set_atoms(atoms1, atoms2)
sup.apply(model2.get_atoms())
print("RMSD after superimposition:", sup.rms)

# Create new structure with 2 models
combined_structure = Structure.Structure("combined")

# Add first model as model 0
combined_structure.add(model1)

# Rename second model to 1 and add (Bio.PDB requires unique model ids)
model2.id = 1
combined_structure.add(model2)

# Save combined multi-model CIF
io = MMCIFIO()
io.set_structure(combined_structure)
io.save("combined_aligned.cif")


# In[ ]:


view = nv.show_file("combined_aligned.cif")

slider = widgets.IntSlider(min=0, max=100, step=1, value=50, description='My Slider')

view.clear_representations()

# Cartoon for chain CY in silver
view.add_representation('cartoon', selection=':CY', color='silver')

# Cartoon for chain BZ in pink
view.add_representation('cartoon', selection=':BZ', color='pink')

# Highlight residue 648 in chain CY as red ball+stick + label
selection_string = "648:CY"
view.add_representation('ball+stick', selection=selection_string, color='red')
view.add_representation('label', sele=selection_string, labelType='format', labelFormat='%(resname)s%(resno)s', color='red', xOffset=1, fixedSize=True)

display(view)
display(slider)



# In[ ]:


# Load CIF file from AlphaFold 3 server
view = nv.show_file("e.coli36099ef-g.cif")
view.clear_representations()
view.add_representation('cartoon', color='pink')
selection_string = "659:A.CA"
view.add_representation('ball+stick', selection=selection_string, color='red')
view.add_representation('label', sele=selection_string, labelType='format', labelFormat='%(resname)s%(resno)s', color='red', xOffset=1, fixedSize=True)

display(view)


# In[ ]:


# Load CIF file from AlphaFold 3 server
view = nv.show_file("e.coli3609910xmicef-g.cif")
view.clear_representations()
view.add_representation('cartoon', color='siver')
selection_string = "659:A.CA"
view.add_representation('ball+stick', selection=selection_string, color='red')
view.add_representation('label', sele=selection_string, labelType='format', labelFormat='%(resname)s%(resno)s', color='red', xOffset=1, fixedSize=True)

display(view)


# In[ ]:


# Parse structures
parser = MMCIFParser()
structure1 = parser.get_structure('s1', 'e.coli36099ef-g.cif')
structure2 = parser.get_structure('s2', 'e.coli3609910xmicef-g.cif')

model1 = structure1[0]
model2 = structure2[0]

chain1_id = 'A'
chain2_id = 'A'

atoms1 = [res['CA'] for res in model1[chain1_id] if 'CA' in res]
atoms2 = [res['CA'] for res in model2[chain2_id] if 'CA' in res]

min_len = min(len(atoms1), len(atoms2))
atoms1 = atoms1[:min_len]
atoms2 = atoms2[:min_len]

# Superimpose
sup = Superimposer()
sup.set_atoms(atoms1, atoms2)
sup.apply(model2.get_atoms())
print("RMSD after superimposition:", sup.rms)

# Create new structure with 2 models
combined_structure = Structure.Structure("combined")

# Add first model as model 0
combined_structure.add(model1)

# Rename second model to 1 and add (Bio.PDB requires unique model ids)
model2.id = 1
combined_structure.add(model2)

# Save combined multi-model CIF
io = MMCIFIO()
io.set_structure(combined_structure)
io.save("combined_alignedAF3model.cif")


# In[ ]:


view = nv.show_file("combined_alignedAF3model.cif")
view.clear_representations()

# Try different model selection syntaxes
try:
    view.add_representation('cartoon', selection='/0', color='pink')      # Model 0 with slash
    view.add_representation('cartoon', selection='/1', color='silver')    # Model 1 with slash
    print("Using model /0 and /1 - SUCCESS")
except Exception as e:
    print(f"Model selection failed: {e}")
    
    # Fallback: try chain selection
    try:
        view.add_representation('cartoon', selection=':A', color='pink')
        view.add_representation('cartoon', selection=':B', color='silver')
        print("Using chain A and B - SUCCESS")
    except Exception as e2:
        print(f"Chain selection failed: {e2}")
        
        # Final fallback: show everything with automatic coloring
        view.add_representation('cartoon', color='chainid')  # Colors by chain automatically
        print("Using automatic chain coloring - SUCCESS")

# Highlight residue
selection_string = "659:A"
try:
    view.add_representation('ball+stick', selection=selection_string, color='red')
    view.add_representation('label', 
                           selection=selection_string, 
                           labelType='format', 
                           labelFormat='%(resname)s%(resno)s', 
                           color='red', 
                           xOffset=1, 
                           fixedSize=True)
except Exception as e:
    print(f"Residue highlighting failed: {e}")

# Force display
display(view)

# Slider
slider = widgets.IntSlider(min=0, max=100, step=1, value=50, description='My Slider')
display(slider)


# In[ ]:


# Load CIF file from AlphaFold 3 server
view = nv.show_file("e.coli36099ef-g.cif")
view.clear_representations()
view.add_representation('cartoon', color='pink')
selection_string = "593:A.CA"
view.add_representation('ball+stick', selection=selection_string, color='red')
view.add_representation('label', sele=selection_string, labelType='format', labelFormat='%(resname)s%(resno)s', color='red', xOffset=1, fixedSize=True)

display(view)


# In[ ]:


# Load CIF file from AlphaFold 3 server
view = nv.show_file("ef-ge.colimg165510xmic.cif")
view.clear_representations()
view.add_representation('cartoon', color='siver')
selection_string = "593:A.CA"
view.add_representation('ball+stick', selection=selection_string, color='red')
view.add_representation('label', sele=selection_string, labelType='format', labelFormat='%(resname)s%(resno)s', color='red', xOffset=1, fixedSize=True)

display(view)


# In[ ]:


# Parse structures
parser = MMCIFParser()
structure1 = parser.get_structure('s1', 'e.coli36099ef-g.cif')
structure2 = parser.get_structure('s2', 'ef-ge.colimg165510xmic.cif')

model1 = structure1[0]
model2 = structure2[0]

chain1_id = 'A'
chain2_id = 'A'

atoms1 = [res['CA'] for res in model1[chain1_id] if 'CA' in res]
atoms2 = [res['CA'] for res in model2[chain2_id] if 'CA' in res]

min_len = min(len(atoms1), len(atoms2))
atoms1 = atoms1[:min_len]
atoms2 = atoms2[:min_len]

# Superimpose
sup = Superimposer()
sup.set_atoms(atoms1, atoms2)
sup.apply(model2.get_atoms())
print("RMSD after superimposition:", sup.rms)

# Create new structure with 2 models
combined_structure = Structure.Structure("combined")

# Add first model as model 0
combined_structure.add(model1)

# Rename second model to 1 and add (Bio.PDB requires unique model ids)
model2.id = 1
combined_structure.add(model2)

# Save combined multi-model CIF
io = MMCIFIO()
io.set_structure(combined_structure)
io.save("combined_alignedAF3model2.cif")


# In[ ]:


view = nv.show_file("combined_alignedAF3model2.cif")
view.clear_representations()

# Try different model selection syntaxes
try:
    view.add_representation('cartoon', selection='/0', color='pink')      # Model 0 with slash
    view.add_representation('cartoon', selection='/1', color='silver')    # Model 1 with slash
    print("Using model /0 and /1 - SUCCESS")
except Exception as e:
    print(f"Model selection failed: {e}")
    
    # Fallback: try chain selection
    try:
        view.add_representation('cartoon', selection=':A', color='pink')
        view.add_representation('cartoon', selection=':B', color='silver')
        print("Using chain A and B - SUCCESS")
    except Exception as e2:
        print(f"Chain selection failed: {e2}")
        
        # Final fallback: show everything with automatic coloring
        view.add_representation('cartoon', color='chainid')  # Colors by chain automatically
        print("Using automatic chain coloring - SUCCESS")

# Highlight residue
selection_string = "593:A"
try:
    view.add_representation('ball+stick', selection=selection_string, color='red')
    view.add_representation('label', 
                           selection=selection_string, 
                           labelType='format', 
                           labelFormat='%(resname)s%(resno)s', 
                           color='red', 
                           xOffset=1, 
                           fixedSize=True)
except Exception as e:
    print(f"Residue highlighting failed: {e}")

# Force display
display(view)
# Slider
slider = widgets.IntSlider(min=0, max=100, step=1, value=50, description='My Slider')
display(slider)


# In[ ]:


parser = MMCIFParser()
structure = parser.get_structure("4v5f", "4v5f.cif") 


# In[ ]:


# Load the structure
parser = MMCIFParser()
structure = parser.get_structure("4v5f", "4v5f.cif")

# Custom Select class to keep only EF-G_S12 chains
class EF_G_S12_ChainSelect(Select):
    def __init__(self, chain_ids):
        self.chain_ids = chain_ids

    def accept_chain(self, chain):
        return chain.id in self.chain_ids

# Set up the writer
io = MMCIFIO()
io.set_structure(structure)

# Save only EF-G_S12 chains 
io.save("4v5f_EF-G_S12.cif", EF_G_S12_ChainSelect(["CY", "CL"]))


# In[ ]:


## for finding proline number in 4v5f structure:
# Parse the structure
parser = MMCIFParser()
structure = parser.get_structure('4V5F', '4v5f_EF-G_S12.cif')

# Iterate over models, chains, and residues to list amino acids and their residue numbers
for model in structure:
    for chain in model:
        print(f"Chain ID: {chain.id}")
        for residue in chain:
            # residue.id is a tuple like (' ', resseq, insertion_code)
            resseq = residue.id[1]  # residue sequence number
            resname = residue.get_resname()  # 3-letter amino acid code
            print(f"Residue {resseq}: {resname}")
        print()  # Blank line between chains


# In[ ]:


parser = MMCIFParser()
structure = parser.get_structure("4v53", "4v53.cif") 


# In[ ]:


mmcif_dict = MMCIF2Dict("4v53.cif")

# This will list all chain IDs and associated entity IDs
strand_ids = mmcif_dict.get("_entity_poly.pdbx_strand_id", [])
entity_ids = mmcif_dict.get("_entity_poly.entity_id", [])
polymer_names = mmcif_dict.get("_entity.pdbx_description", [])

for entity, strands, name in zip(entity_ids, strand_ids, polymer_names):
    print(f"Entity ID: {entity}, Chains: {strands}, Molecule: {name}")


# In[ ]:


# Load the structure
parser = MMCIFParser()
structure = parser.get_structure("4v53", "4v53.cif")

# Custom Select class to keep only S12 chains
class S12_ChainSelect(Select):
    def __init__(self, chain_ids):
        self.chain_ids = chain_ids

    def accept_chain(self, chain):
        return chain.id in self.chain_ids

# Set up the writer
io = MMCIFIO()
io.set_structure(structure)

# Save only S12 chains 
io.save("4v53_S12.cif", S12_ChainSelect(["CL"]))


# In[ ]:


parser = FastMMCIFParser(QUIET=True)
structure1 = parser.get_structure('s1', '4v5f_EF-G_S12.cif')
structure2 = parser.get_structure('s2', '4v53_S12.cif')


model1 = structure1[0]
model2 = structure2[0]

chain1_id = ["CY", "CL"]  
chain2_id = "CL"

atoms1 = [res["CA"] 
          for cid in chain1_id 
          for res in model1[cid] 
          if "CA" in res]

atoms2 = [res['CA'] for res in model2[chain2_id] if 'CA' in res]

min_len = min(len(atoms1), len(atoms2))
atoms1 = atoms1[:min_len]
atoms2 = atoms2[:min_len]

# Superimpose
sup = Superimposer()
sup.set_atoms(atoms1, atoms2)
sup.apply(model2.get_atoms())
print("RMSD after superimposition:", sup.rms)

# Create new structure with 2 models
combined_structure = Structure.Structure("combined")

# Add first model as model 0
combined_structure.add(model1)

# Rename second model to 1 and add (Bio.PDB requires unique model ids)
model2.id = 1
combined_structure.add(model2)

# Save combined multi-model CIF
io = MMCIFIO()
io.set_structure(combined_structure)
io.save("alignedS12+EF-G.cif")


# In[ ]:


# Rename overlapping chains in model2 to avoid coloring clashes
existing_chains = {chain.id for chain in model1}
for chain in model2:
    if chain.id in existing_chains:
        chain.id = chain.id + "_2"  # e.g., "CL" -> "CL_2"

# Create combined structure
combined_structure = Structure.Structure("combined")
combined_structure.add(model1)
model2.id = 1  # second model
combined_structure.add(model2)

# Save combined CIF
io = MMCIFIO()
io.set_structure(combined_structure)
io.save("alignedS12+EF-G_renamed.cif")


# In[ ]:


view = nv.show_file("alignedS12+EF-G_renamed.cif")
view.clear_representations()

#color by chain names
view.add_representation("cartoon", selection=":CY", color="pink")
view.add_representation("cartoon", selection=":CL", color="blue")
view.add_representation("cartoon", selection=":CL_2", color="silver")

from IPython.display import display
display(view)


# In[ ]:


# Load the previously superimposed structure
parser = MMCIFParser(QUIET=True)
aligned_structure = parser.get_structure('aligned', 'alignedS12+EF-G_renamed.cif')
model2 = aligned_structure[0]  # first model (contains CL_2)

# Load the new structure
new_structure = parser.get_structure('4V53', '4V53.cif')
model_new = new_structure[0]  # first model of 4V53

# Select chains to use for superimposition
chain_old = "CL_2"  # chain from previously superimposed structure
chain_new = "CL"    # chain in 4V53

# Extract CA atoms for alignment
atoms_old = [res["CA"] for res in model2[chain_old] if "CA" in res]
atoms_new = [res["CA"] for res in model_new[chain_new] if "CA" in res]

# Match lengths
min_len = min(len(atoms_old), len(atoms_new))
atoms_old = atoms_old[:min_len]
atoms_new = atoms_new[:min_len]

# Superimpose
sup = Superimposer()
sup.set_atoms(atoms_old, atoms_new)
sup.apply(model_new.get_atoms())
print("RMSD after superimposition:", sup.rms)

# Rename chain in new model if needed to avoid duplicates
if chain_new in [c.id for c in model2]:
    for chain in model_new:
        if chain.id == chain_new:
            chain.id = chain_new + "_new"

# Assign unique model IDs before combining
model2.id = 0
model_new.id = 1

# Combine structures into one
combined_structure = Structure.Structure("super_aligned_2")
combined_structure.add(model2)
combined_structure.add(model_new)

# Save combined structure
io = MMCIFIO()
io.set_structure(combined_structure)
io.save("superimposed_combined.cif")


# In[ ]:


# Load the CIF file directly
view = nv.show_file("superimposed_combined.cif")

# Clear default representations if you want custom colors
view.clear_representations()

# Add cartoon representations for chains
view.add_representation("cartoon", selection=":CY", color="pink")




# Parse the CIF file
parser = MMCIFParser()
structure = parser.get_structure("super", "superimposed_combined.cif")

# Extract chain IDs for model 1
model_index = 0  # Python index
chains_model1 = [chain.id for chain in structure[model_index]]

# Assign colors
colors = [plt.cm.tab20(i % 20) for i in range(len(chains_model1))]

# Load into NGLView directly
view = nv.show_file("superimposed_combined.cif")
view.clear_representations()

# Add cartoon representation for each chain with colors
for chain_id, color in zip(chains_model1, colors):
    view.add_representation("cartoon", selection=f":{chain_id}", color=color)

from IPython.display import display
display(view)


# In[ ]:


# Parse structure and identify all hetero atoms (potential ligands)
parser = MMCIFParser()
structure = parser.get_structure("super", "superimposed_combined.cif")

# Find all hetero atoms/ligands
ligands = []
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.get_id()[0] != ' ':  # hetero atoms
                if residue.get_id()[0] not in ['W', 'H_']:  # exclude water and hydrogens
                    ligands.append((model.id, chain.id, residue.get_resname(), residue.get_id()))

print("Found ligands/hetero atoms:")
for model_id, chain_id, resname, res_id in ligands:
    print(f"Model {model_id}, Chain {chain_id}: {resname} {res_id}")



# Print detailed residue information for debugging
print("\nDetailed residue information:")
for model in structure:
    print(f"Model {model.id}:")
    for chain in model:
        print(f"  Chain {chain.id}:")
        hetero_residues = [res for res in chain if res.get_id()[0] != ' ']
        if hetero_residues:
            for res in hetero_residues:
                print(f"    {res.get_resname()} {res.get_id()}")
        else:
            print("    No hetero residues found")


# In[ ]:


# Parse the CIF file to get chain information
parser = MMCIFParser()
structure = parser.get_structure("super", "superimposed_combined.cif")

# Get all chain IDs from both models
all_chains = []
for model in structure:
    for chain in model:
        all_chains.append(chain.id)

print("Available chains:", all_chains)

# Load the structure in NGLView
view = nv.show_file("superimposed_combined.cif")
view.clear_representations()

# Define colors for each chain (using NGLView compatible color names)
chain_colors = {
   'CY': 'pink', 'CL_new': 'green',
    # Add more chains with different colors as needed
    'AA': 'orange',
    'DB': 'green',
    'CA': 'purple',
    'BB': 'cyan',
    'BA': 'red',
    'BB': 'brown',
    'BL': 'gray',
    'DB': 'magenta'
}


# Add cartoon representations for protein chains
for chain_id in all_chains:
    if chain_id in chain_colors:
        color = chain_colors[chain_id]
    else:
        # Use a default color for unlisted chains
        color = 'green'
    
    view.add_representation("cartoon", 
                           selection=f":{chain_id} and polymer", 
                           color=color,
                           opacity=0.8)
# Highlight PHE:582 in chain CY as red stick
selection_string = "582:CY"
view.add_representation('ball+stick', selection=selection_string, color='red')
view.add_representation('label', sele=selection_string, labelType='format', labelFormat='%(resname)s%(resno)s', color='red', xOffset=1, fixedSize=True)

print("Added PHE:582 in chain CY as green stick")

# Highlight gentamicin (LLL) in red
# Try different possible selections for the ligand
gentamicin_selections = [
    "LLL",  # residue name
    "[LLL]",  # alternative syntax
    "hetero and not water",  # all hetero atoms except water
    "ligand"  # generic ligand selection
]

for selection in gentamicin_selections:
    try:
        view.add_representation("ball+stick", 
                               selection=selection, 
                               color="red",
                               radius=0.5)
        print(f"Added gentamicin representation with selection: {selection}")
        break
    except:
        continue

#Optional: Add surface representation for better visualization
#view.add_representation("surface", 
#                       selection="protein", 
#                        opacity=0.3,
#                        color="white")

# Set background and camera
view.background = "white"
view.camera = "perspective"

# Add labels (optional)
view.add_label(text="Gentamicin", selection="LLL", color="black", size=1)

#RNA visualisation
#view.add_representation("tube", selection="nucleic", color="orange", opacity=0.8)


display(view)

# Print information about the structure
print("\nStructure information:")
for i, model in enumerate(structure):
    print(f"Model {i}:")
    for chain in model:
        residue_count = len([res for res in chain if res.get_id()[0] == ' '])  # protein residues
        hetero_count = len([res for res in chain if res.get_id()[0] != ' '])   # hetero atoms
        print(f"  Chain {chain.id}: {residue_count} protein residues, {hetero_count} hetero atoms")


# In[ ]:




