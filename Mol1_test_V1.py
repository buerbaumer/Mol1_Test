import streamlit as st
from stmol import showmol
import py3Dmol
from urllib.request import urlopen
from urllib.parse import quote

# demo only 

from rdkit import Chem
from rdkit.Chem import AllChem

display_on = False

# sidebar:
style_choosen = st.sidebar.selectbox('style',['stick','balls & stick', 'sphere', 'cartoon', 'clicksphere', 'line'])
add_Hs = st.sidebar.checkbox('add Hs', value = True)
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
    #if add_Hs == 'True':
    #    mol = Chem.AddHs(mol)
    #else:
    #    mol = Chem.RemoveHs(mol)
    mol = Chem.AddHs(mol)   
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    mblock = Chem.MolToMolBlock(mol)
    return mblock

#def calc_rings(comp_str):
#    compound_rings_calc = 0
#    # calc the rings of the compound :-)
#    # note: compound_rings_calc is global
#    for x in comp_str:
#        if x.isdigit() == True:
#            x_int = int(x)
#            if x_int > compound_rings_calc:
#                compound_rings_calc = x_int
#    return compound_rings_calc

#def CIRconvert(ids):
#    try:
#        # url = 'http://cactus.nci.nih.gov/chemical/structure/' + quote(ids) + '/smiles'
#        url = 'http://cactus.nci.nih.gov/chemical/structure/' + quote(ids) + '/smiles'
#        ans = urlopen(url).read().decode('utf8')
#        display_on = True
#        return ans
#    except:
#        display_on = False
#        return 'Sorry, this structure could not be found.'
 
#def CIRconvert_MW(ids):
#    try:
#        url = 'http://cactus.nci.nih.gov/chemical/structure/' + quote(ids) + '/mw'
#        ansRing = urlopen(url).read().decode('utf8')
#        display_on = True
#        return ansRing
#    except:
#        display_on = False
#        return 'Sorry, this structure could not be found.'
    
def render_mol(M3D_input):
    with col2:
        M3Dview = py3Dmol.view()
        M3Dview.addModel(M3D_input,'mol')
                                              
        if style_choosen == 'balls & stick':
            M3Dview.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}})
        else:
            M3Dview.setStyle({style_choosen: {}})
                                              
        M3Dview.setBackgroundColor(color_b)
        #M3Dview.setBackgroundColor('white')
    
        if spin:
            M3Dview.spin(True)
        else:
            M3Dview.spin(False)
        
        M3Dview.zoomTo()
        showmol(M3Dview, height=500, width=500)

#compound_smiles = CIRconvert(compound_input)
#compound_mw = CIRconvert_MW(compound_input)
#compound_rings = calc_rings(compound_smiles)

with col1:
#    #st.write("1")
#    st.write('')
#    st.markdown('**SMILES:**')
#    st.write(compound_smiles)
#    st.write('')
#    st.markdown('**Molecular Weight:**')
#    st.write(compound_mw)
#    st.markdown('**Number of rings:**')
#    st.write(compound_rings)
#    st.write('')
#    st.write('')
    protein_input = st.text_input('Input the name of chemical structure:','121P')
    
with col2:
    #st.write("2")

#    if "Sorry" not in compound_smiles:
#        blk=makeblock(compound_smiles)
#        render_mol(blk)

    #https://3dmol.csb.pitt.edu/viewer.html?pdb=1YCR
    #& select=chain:A
    #& style=cartoon;stick:radius~0.1
    #& surface=opacity:0.8;
    # colorscheme:whiteCarbon
    #& select=chain:B
    #& style=cartoon;line
    #& select=resi:19,23,26;chain:B
    #& style=cartoon;stick
    #& labelres=backgroundOpacity:0.8;fontSize:14
    
    prot_str = protein_input
    M3Dview2 = py3Dmol.view(query='pdb:'+prot_str)
    #M3Dview2.setStyle({style_choosen:{'color':'spectrum'}})
    
    if style_choosen == 'balls & stick':
        M3Dview2.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}})
    else:
        M3Dview2.setStyle({style_choosen: {}})
        
    M3Dview2.setBackgroundColor(color_b)
    if spin:
        M3Dview2.spin(True)
    else:
        M3Dview2.spin(False)
    M3Dview2.zoomTo()
    showmol(M3Dview2,height=500,width=800) 
