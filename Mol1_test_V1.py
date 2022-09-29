import streamlit as st
from stmol import showmol
import py3Dmol
from urllib.request import urlopen
from urllib.parse import quote

# demo only 

from rdkit import Chem
from rdkit.Chem import AllChem

display_on = False

st.title('Analyser for Chemical Structures (ACS)')
st.write('using py3Dmol, stmol, rdkit, streamlit and data from the National Cancer Institute (https://www.cancer.gov/) ')

def makeblock(smi):
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mblock = Chem.MolToMolBlock(mol)
    return mblock

def CIRconvert(ids):
    try:
        url = 'http://cactus.nci.nih.gov/chemical/structure/' + quote(ids) + '/smiles'
        ans = urlopen(url).read().decode('utf8')
        display_on = True
        return ans
    except:
        display_on = False
        return 'Sorry, this structure could not be found.'
        
    
def render_mol(xyz):
    xyzview = py3Dmol.view()#(width=600,height=400)
    xyzview.addModel(xyz,'mol')
    # xyzview.setStyle({style_choosen: {'radius': 0.1}, 'sphere': {'scale': 0.25}})
    xyzview.setStyle({style_choosen: {}})
    # xyzview.setStyle({'model': -1}, {"cartoon": {'color': 'spectrum'}})
    xyzview.setBackgroundColor(color_b)
    #xyzview.setBackgroundColor('white')
    
    if spin:
        xyzview.spin(True)
    else:
        xyzview.spin(False)
        
    xyzview.zoomTo()
    showmol(xyzview,height=500,width=500)

#compound_smiles = st.text_input('Input SMILES','COc1ccc2[nH]c([S@@+]([O-])Cc3ncc(C)c(OC)c3C)nc2c1')
#compound_input = st.sidebar.selectbox('Input the name of chemical structure: ',['3-Methylheptane', 'Aspirin', 'Diethylsulfate', 'Diethyl sulfate', '50-78-2', 'Adamant'])
compound_input = st.text_input('Input the name of chemical structure:','Aspirin')
compound_smiles = CIRconvert(compound_input)
st.write(compound_smiles)

style_choosen = st.sidebar.selectbox('style',['stick','sphere','cartoon','clicksphere', 'line'])
spin = st.sidebar.checkbox('Spin', value = True)
color_b = st.sidebar.color_picker('Pick Background color', '#000000')
st.sidebar.write('The current color is', color_b)

# st.write(display_on)

if "Sorry" not in compound_smiles:
    blk=makeblock(compound_smiles)
    render_mol(blk)


# ---------------------------------------
# import py3Dmol
#view = py3Dmol.view(query='pdb:121P')
#view.setStyle({'cartoon':{'color':'spectrum'}})
#view 
