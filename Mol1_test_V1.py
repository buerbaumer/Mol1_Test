import streamlit as st
from stmol import showmol
import py3Dmol
from urllib.request import urlopen
from urllib.parse import quote

# demo only 

from rdkit import Chem
from rdkit.Chem import AllChem

display_on = False
compound_rings_calc = 0

# sidebar:
style_choosen = st.sidebar.selectbox('style',['stick','sphere','cartoon','clicksphere', 'line'])
spin = st.sidebar.checkbox('Spin', value = True)
color_b = st.sidebar.color_picker('Pick Background color', '#ffffff')
st.sidebar.write('The current color is', color_b)

# Main page:
st.title('Analyser for Chemical Structures (ACS)')
st.write('using py3Dmol, stmol, rdkit, streamlit and data from the National Cancer Institute (https://www.cancer.gov/) ')
#compound_smiles = st.text_input('Input SMILES','COc1ccc2[nH]c([S@@+]([O-])Cc3ncc(C)c(OC)c3C)nc2c1')
#compound_input = st.sidebar.selectbox('Input the name of chemical structure: ',['3-Methylheptane', 'Aspirin', 'Diethylsulfate', 'Diethyl sulfate', '50-78-2', 'Adamant'])

col1, col2 = st.columns([1,3])
with col1:
    compound_input = st.text_input('Input the name of chemical structure:','Aspirin')

def makeblock(smi):
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mblock = Chem.MolToMolBlock(mol)
    return mblock

def CIRconvert(ids):
    try:
        # url = 'http://cactus.nci.nih.gov/chemical/structure/' + quote(ids) + '/smiles'
        url = 'http://cactus.nci.nih.gov/chemical/structure/' + quote(ids) + '/smiles'
        ans = urlopen(url).read().decode('utf8')
        display_on = True
        # calc compound_rings_calc

        for x in ans:
            if x.isdigit() == True:
            x_int = int(x)
            if x_int > compound_rings_calc:
                compound_rings_calc = x_int

        return ans
    except:
        display_on = False
        return 'Sorry, this structure could not be found.'
 
def CIRconvert_MW(ids):
    try:
        url = 'http://cactus.nci.nih.gov/chemical/structure/' + quote(ids) + '/mw'
        ansRing = urlopen(url).read().decode('utf8')
        display_on = True
        return ansRing
    except:
        display_on = False
        return 'Sorry, this structure could not be found.'
    
def render_mol(xyz):
    with col2:
        xyzview = py3Dmol.view()
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
        showmol(xyzview, height=500, width=500)

compound_smiles = CIRconvert(compound_input)
compound_mw = CIRconvert_MW(compound_input)

with col1:
    #st.write("1")
    st.write('')
    st.markdown('**SMILES:**')
    st.write(compound_smiles)
    st.write('')
    st.markdown('**Molecular Weight:**')
    st.write(compound_mw)
    st.markdown('**Number of rings:**')
    st.write(compound_rings_calc)
    
with col2:
    #st.write("2")

    if "Sorry" not in compound_smiles:
        blk=makeblock(compound_smiles)
        render_mol(blk)


# ---------------------------------------
# import py3Dmol
#view = py3Dmol.view(query='pdb:121P')
#view.setStyle({'cartoon':{'color':'spectrum'}})
#view 
