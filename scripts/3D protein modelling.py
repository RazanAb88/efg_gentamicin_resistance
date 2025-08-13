#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Bio.PDB import MMCIFParser, MMCIFIO


# In[2]:


parser = MMCIFParser()
structure = parser.get_structure("4WPO", "4wpo.cif") 


# In[4]:


from Bio.PDB import PDBList
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

mmcif_dict = MMCIF2Dict("4wpo.cif")

# This will list all chain IDs and associated entity IDs
strand_ids = mmcif_dict.get("_entity_poly.pdbx_strand_id", [])
entity_ids = mmcif_dict.get("_entity_poly.entity_id", [])
polymer_names = mmcif_dict.get("_entity.pdbx_description", [])

for entity, strands, name in zip(entity_ids, strand_ids, polymer_names):
    print(f"Entity ID: {entity}, Chains: {strands}, Molecule: {name}")


# In[5]:


from Bio.PDB import MMCIFParser, MMCIFIO, Select

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


# In[7]:


import nglview as nv
import ipywidgets as widgets
from IPython.display import display
widgets.IntSlider()
view = nv.show_file("4WPO_EF-G.cif")
slider = widgets.IntSlider(min=0, max=100, step=1, value=50, description='My Slider')
view.clear_representations()
view.add_representation('cartoon', color='pink')  
view.add_representation('stick', selection='540 and PRO', color='black')
display(view)
display(slider)


# In[8]:


from Bio.PDB import MMCIFParser, PPBuilder
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

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


# In[10]:


parser = MMCIFParser()
structure = parser.get_structure("4v5f", "4v5f.cif") 


# In[11]:


from Bio.PDB.MMCIF2Dict import MMCIF2Dict

mmcif_dict = MMCIF2Dict("4v5f.cif")

# This will list all chain IDs and associated entity IDs
strand_ids = mmcif_dict.get("_entity_poly.pdbx_strand_id", [])
entity_ids = mmcif_dict.get("_entity_poly.entity_id", [])
polymer_names = mmcif_dict.get("_entity.pdbx_description", [])

for entity, strands, name in zip(entity_ids, strand_ids, polymer_names):
    print(f"Entity ID: {entity}, Chains: {strands}, Molecule: {name}")


# In[12]:


from Bio.PDB import MMCIFParser, MMCIFIO, Select

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


# In[14]:


## for finding proline number in 4v5f structure:
from Bio.PDB import MMCIFParser

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


# In[13]:


import nglview as nv
import ipywidgets as widgets
from IPython.display import display
view = nv.show_file("4v5f_EF-G.cif")
slider = widgets.IntSlider(min=0, max=100, step=1, value=50, description='My Slider')
view.clear_representations()
view.add_representation('cartoon', color='silver')
selection_string = "648:CY"
view.add_representation('ball+stick', selection=selection_string)
view.add_representation('label', sele=selection_string, labelType='format', labelFormat='%(resname)s%(resno)s', color='red', xOffset=1, fixedSize=True)




display(view)
display(slider)


# In[15]:


from Bio.PDB import MMCIFParser, Superimposer, MMCIFIO, Model, Structure

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


# In[17]:


import nglview as nv
import ipywidgets as widgets
from IPython.display import display

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



# In[37]:


import nglview as nv
import ipywidgets as widgets
from IPython.display import display

# Load CIF file from AlphaFold 3 server
view = nv.show_file("e.coli36099ef-g.cif")
view.clear_representations()
view.add_representation('cartoon', color='pink')
selection_string = "659:A.CA"
view.add_representation('ball+stick', selection=selection_string, color='red')
view.add_representation('label', sele=selection_string, labelType='format', labelFormat='%(resname)s%(resno)s', color='red', xOffset=1, fixedSize=True)

display(view)


# In[61]:


import nglview as nv
import ipywidgets as widgets
from IPython.display import display

# Load CIF file from AlphaFold 3 server
view = nv.show_file("e.coli3609910xmicef-g.cif")
view.clear_representations()
view.add_representation('cartoon', color='siver')
selection_string = "659:A.CA"
view.add_representation('ball+stick', selection=selection_string, color='red')
view.add_representation('label', sele=selection_string, labelType='format', labelFormat='%(resname)s%(resno)s', color='red', xOffset=1, fixedSize=True)

display(view)


# In[105]:


from Bio.PDB import MMCIFParser, Superimposer, MMCIFIO, Model, Structure

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


# In[136]:


import nglview as nv
import ipywidgets as widgets
from IPython.display import display

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




