import streamlit as st
from stmol import showmol
import py3Dmol

# demo only

from rdkit import Chem
from rdkit.Chem import AllChem

st.title('RDKit + Py3DMOL')

def makeblock(smi):
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mblock = Chem.MolToMolBlock(mol)
    return mblock

def render_mol(xyz):
    xyzview = py3Dmol.view()#(width=400,height=400)
    xyzview.addModel(xyz,'mol')
    xyzview.setStyle({'stick':{}})
    # xyzview.setStyle({'model': -1}, {"cartoon": {'color': 'spectrum'}})
    xyzview.setBackgroundColor('white')
    xyzview.zoomTo()
    showmol(xyzview,height=500,width=500)

compound_smiles=st.text_input('Input SMILES','COc1ccc2[nH]c([S@@+]([O-])Cc3ncc(C)c(OC)c3C)nc2c1')
blk=makeblock(compound_smiles)
render_mol(blk)

# ---------------------------------------
# import py3Dmol
view = py3Dmol.view(query='pdb:121P')
view.setStyle({'cartoon':{'color':'spectrum'}})
view
