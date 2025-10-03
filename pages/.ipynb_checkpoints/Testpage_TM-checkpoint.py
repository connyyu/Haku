import streamlit as st
import py3Dmol
from stmol import showmol
import requests
import re
import pandas as pd

st.set_page_config(page_title="Haku - find TM regions", page_icon="💮")

def fetch_uniprot_sequence(uniprot_ac):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        sequence = response.text.split('\n', 1)[1].replace('\n', '').replace('\r', '')
        return sequence
    else:
        st.error(f"Failed to fetch sequence for UniProt ID: {uniprot_ac}")
        return None

def fetch_uniprot_tm_helices(uniprot_ac):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}.txt"
    response = requests.get(url)
    tm_helices = []
    out_seq = []  
    in_seq = []   
    if response.status_code == 200:
        lines = response.text.splitlines()
        i = 0
        while i < len(lines):
            line = lines[i]
            tm_match = re.search(r"FT\s+TRANSMEM\s+(\d+)\.+(\d+)", line)
            if tm_match:
                start, end = map(int, tm_match.groups())
                tm_helices.append((start, end))
            topo_match = re.search(r"FT\s+TOPO_DOM\s+(\d+)\.+(\d+)", line)
            if topo_match and i + 1 < len(lines):
                start, end = map(int, topo_match.groups())
                next_line = lines[i + 1]
                if "Extracellular" in next_line:
                    out_seq.append((start, end))
                elif "Cytoplasmic" in next_line:
                    in_seq.append((start, end))
            i += 1
    else:
        st.error(f"Failed to fetch TM helices for UniProt ID: {uniprot_ac}")
    return tm_helices, out_seq, in_seq

def convert_tm_helices_to_pred(sequence, tm_helices, out_seq, in_seq):
    pred_uniprot = ["g"] * len(sequence)  
    for start, end in tm_helices:
        for i in range(start - 1, end):
            pred_uniprot[i] = "M"
    for start, end in in_seq:
        for i in range(start - 1, end):
            pred_uniprot[i] = "I"
    for start, end in out_seq:
        for i in range(start - 1, end):
            pred_uniprot[i] = "O"
    return "".join(pred_uniprot)

def fetch_pdb_structure(pdb_code):
    pdb_code_dl = pdb_code.lower()
    url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_code_dl}_updated.cif"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        st.error(f"Failed to fetch structure for PDB code: {pdb_code}")
        return None

def viewpdb(structure, pred, sequence):
    st.markdown('<style>div.block-container{padding-top:3rem;}</style>', unsafe_allow_html=True)
    view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js', height=400, width=400)
    view.addModelsAsFrames(structure)
    view.setBackgroundColor('#eeeeee')
    view.spin(False)

    atom_color = dict()
    for nr, res_type in enumerate(pred):
        if res_type == 'O':
            atom_color[nr] = 'powderblue'
        elif res_type == 'M':
            atom_color[nr] = '#830592'
        elif res_type == 'I':
            atom_color[nr] = 'pink'
        elif res_type == 'g':
            atom_color[nr] = 'lightgrey'
        else:
            atom_color[nr] = '#008c74'

    view.setStyle({'model': -1}, {
        'cartoon': {
            'thickness': 0.5,
            'colorscheme': {'prop': 'resi', 'map': atom_color}
        }
    })
    view.zoomTo()
    showmol(view, height=400, width=600)
    st.caption(f"Visualizing predicted transmembrane regions for UniProt AC: {uniprot_ac}.")

default_unp = 'Q63008'
default_pdb = '7UV0'

with st.form(key='uniprot_form'):
    col1, col2, col3 = st.columns([1, 2, 1])
    with col1:
        uniprot_ac = st.text_input("Enter UniProt AC:", default_unp).strip()
    with col2:
        st.markdown("<h5 style='font-size: 10px;'></h5>", unsafe_allow_html=True)
        fetch_data_button = st.form_submit_button("Fetch data from UniProt")
    with col3:
        st.markdown('')
        st.markdown('')
        st.markdown(f'<a href="https://rest.uniprot.org/uniprotkb/{uniprot_ac}.txt" target="_blank">Open UniProt entry</a>', unsafe_allow_html=True)
        
with st.form(key='pdb_form'):
    col1, col2 = st.columns([1, 3])
    with col1:
        pdb_code = st.text_input("Enter PDB code:", default_pdb).strip()
    with col2:
        st.markdown("<h5 style='font-size: 10px;'></h5>", unsafe_allow_html=True)
        fetch_pdb_button = st.form_submit_button("Show PDB with UniProt annotations")

# Initialize variables
alphafold_structure = None
tm_helices_uniprot = None
tm_helices_pred = None
sequence = None
pdb_structure = None

if 'sequence' not in st.session_state or 'tm_helices_uniprot' not in st.session_state:
    st.session_state.sequence = fetch_uniprot_sequence(default_unp)
    st.session_state.tm_helices_uniprot, st.session_state.out_seq, st.session_state.in_seq = fetch_uniprot_tm_helices(default_unp)

    if st.session_state.sequence and st.session_state.tm_helices_uniprot is not None:
        st.session_state.pred_uniprot = convert_tm_helices_to_pred(
            st.session_state.sequence, 
            st.session_state.tm_helices_uniprot, 
            st.session_state.out_seq, 
            st.session_state.in_seq
        )
    else:
        st.warning("Sequence or TM annotations not found.")
        st.session_state.pred_uniprot = ""

if fetch_pdb_button:
    st.session_state.pdb_code = pdb_code
    pdb_structure = fetch_pdb_structure(pdb_code)
    if pdb_structure:
        st.session_state.pdb_structure = pdb_structure

if fetch_data_button:
    st.session_state.uniprot_ac = uniprot_ac
    sequence = fetch_uniprot_sequence(uniprot_ac)
    st.session_state.sequence = sequence

    tm_helices_uniprot, out_seq, in_seq = fetch_uniprot_tm_helices(uniprot_ac)
    st.session_state.tm_helices_uniprot = tm_helices_uniprot
    st.session_state.out_seq = out_seq
    st.session_state.in_seq = in_seq

    if sequence and tm_helices_uniprot is not None:
        pred_uniprot = convert_tm_helices_to_pred(sequence, tm_helices_uniprot, out_seq, in_seq)
        st.session_state.pred_uniprot = pred_uniprot
    else:
        st.warning("Sequence or TM annotations not found.")
        st.session_state.pred_uniprot = ""

col1, col2 = st.columns([3, 2])

with col1:
    sequence = st.session_state.get("sequence", None)
    tm_helices_uniprot = st.session_state.get("tm_helices_uniprot", None)
    pred_uniprot = st.session_state.get("pred_uniprot", None)

    pdb_structure = st.session_state.get("pdb_structure", None) or fetch_pdb_structure(default_pdb)

    viewpdb(pdb_structure, pred_uniprot, sequence)

with col2:
    output_str = ""
    tm_helices_uniprot = st.session_state.get("tm_helices_uniprot", [])
    
    if tm_helices_uniprot:
        st.markdown("")
        st.markdown("##### UniProt TM Annotations")
        
        for start, end in tm_helices_uniprot:
            output_str += f"FT TRANSMEM  {start}  {end} Helical\n"
        st.markdown(f"```\n{output_str}\n```")
        st.markdown(
    '<div style="display: flex; padding: 8px; width: 100%;">'
    '<span style="background-color: white; padding: 2px 5px; text-align: center;">Key: </span>'
    '<span style="background-color: powderblue; padding: 2px 5px; text-align: center;">Outside</span>'
    '<span style="background-color: pink; padding: 2px 5px; text-align: center;">Inside</span>'
    '</div>', unsafe_allow_html=True)
    else:
        st.warning("No TM helices found in UniProt data.")
        
# Manual annotation input using text area
col1, col2 = st.columns([3, 2])

# Display the structure with manual annotation
with col1:
    st.markdown("")
    st.markdown("##### Structure")
    
    # Fetch the structure with the current predictions (manual annotations)
    pred_manual = st.session_state.get("pred_manual", None)
    sequence = st.session_state.get("sequence", None)
    pdb_structure = st.session_state.get("pdb_structure", None) or fetch_pdb_structure(default_pdb)
    
    if pred_manual and sequence:
        viewpdb(pdb_structure, pred_manual, sequence)
    else:
        sequence = st.session_state.get("sequence", None)
        pred_empty = ["g"] * len(sequence)
        pdb_structure = st.session_state.get("pdb_structure", None) or fetch_pdb_structure(default_pdb)
        viewpdb(pdb_structure, pred_empty, sequence)

# Manual annotation
with col2:
    output_str = ""
    manual_tm_helices = st.session_state.get("manual_tm_helices", "")
    st.markdown("")
    st.markdown("##### Manual Annotation")
    if tm_helices_uniprot:       
        st.markdown(f"```\n{manual_tm_helices}\n```")
        st.markdown(
    '<div style="display: flex; padding: 8px; width: 100%;">'
    '<span style="background-color: white; padding: 2px 5px; text-align: center;">Key: </span>'
    '<span style="background-color: powderblue; padding: 2px 5px; text-align: center;">Outside</span>'
    '<span style="background-color: pink; padding: 2px 5px; text-align: center;">Inside</span>'
    '</div>', unsafe_allow_html=True)

with st.form(key='manual_annotation_form'):
    col1, col2 = st.columns([3, 2])
    
    tm_helices_uniprot = st.session_state.get("tm_helices_uniprot", [])
    if tm_helices_uniprot:
        for start, end in tm_helices_uniprot:
            output_str += f"FT TRANSMEM  {start}  {end} Helical\n"
        input_tm = (f"{output_str}\n")

    with col1:
        manual_annotations_input = st.text_area(
            "Input manual annotation:",
            height=100,
            value=input_tm
        )
    with col2:
        st.markdown("")
        st.markdown("")
        submit_annotation_button = st.form_submit_button("Update annotation")
        st.markdown("<p style='font-size: 14px;'>Double click to load the new annotation.</p>", unsafe_allow_html=True)
        tm_helices_uniprot = st.session_state.get("tm_helices_uniprot", [])
        if tm_helices_uniprot:
            for start, end in tm_helices_uniprot:
                output_str += f"FT TRANSMEM  {start}  {end} Helical\n"
            input_tm = (f"{output_str}\n")

# Parse the input text for manual annotations if the form is submitted
if submit_annotation_button:
    manual_tm_helices = []
    lines = manual_annotations_input.splitlines()
    for line in lines:
        # Match the FT TRANSMEM format to extract start and end positions
        match = re.match(r"FT TRANSMEM\s+(\d+)\s+(\d+)\s+Helical", line)
        if match:
            start, end = map(int, match.groups())
            manual_tm_helices.append((start, end))
    
    # Only modify the pred_manual for visualization, not the tm_helices_uniprot
    pred_manual = ['g'] * len(st.session_state.sequence)
    for start, end in manual_tm_helices:
        for i in range(start - 1, end):
            pred_manual[i] = 'M'

    st.session_state.manual_tm_helices = manual_annotations_input
    st.session_state.pred_manual = pred_manual
    st.session_state.pdb_structure = st.session_state.get("pdb_structure", None) or fetch_pdb_structure(default_pdb)
