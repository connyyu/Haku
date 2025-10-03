import streamlit as st
import py3Dmol
import streamlit.components.v1 as components
import requests
import subprocess
import os
import glob
import sys
from pathlib import Path
import re
import shutil
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import time

# Sidebar, title, parameters
# -----------------------------------------------------------------------------
st.set_page_config(page_title="Haku - Transmembrane annotations", page_icon="💮")

default_unp = 'Q63008'
default_pdb = '7UV0'

# Get the directory where this script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
# Add the scripts subdirectory to Python path
scripts_dir = os.path.join(script_dir, "scripts")
if scripts_dir not in sys.path:
    sys.path.insert(0, scripts_dir)

# Import the pdb_chainID function from the scripts folder
try:
    import pdb_chainID
except ImportError as e:
    st.error(f"Failed to import pdb_chainID module: {e}")

st.session_state.guide = False

with st.sidebar:
    st.header("📝 New query")
    uniprot_ac = st.text_input("Enter UniProt AC:", default_unp).strip()
    fetch_data_button = st.button("Fetch data")
    pdb_code = st.text_input("Enter PDB code:", default_pdb).strip()
    fetch_pdb_button = st.button("Show structure")
    st.sidebar.markdown("[UniProt annotation](#pdb-uniprot)")
    st.sidebar.markdown("[DeepTMHMM prediction](#pdb-tmhmm)")
    st.sidebar.markdown("[DeepTMHMM plot](#tmhmm_plot)")
    guide_tm = st.checkbox("Instructions", value=st.session_state.get("guide_tm", False))
    st.session_state.guide_tm = guide_tm

st.markdown("<a name='top_title'></a>", unsafe_allow_html=True)
st.markdown("#### Visualise transmembrane annotation on a protein structure.")

script_dir = os.path.dirname(os.path.abspath(__file__)) # location of pages directory
output_dir = os.path.join(script_dir, "biolib_results")
demo_dir = os.path.join(script_dir, "demo_results") # Use demo_results for Q63008
fasta_file = os.path.join(output_dir, "sequence.fasta")

# Instructions
# -----------------------------------------------------------------------------
if st.session_state.get('guide_tm'):
    st.info("""
    In this order:
    1. Enter UniProt AC and **Fetch data**
    2. Enter PDB code and **Show structure**
    
    Run a DeepTMHMM prediction using the **Run TM prediction** button (optional).  
    Use the **Show structure with TM prediction** button to refresh the viewer.  
    """)
    
# Define functions
# -----------------------------------------------------------------------------

# Function to fetch protein sequence from UniProt API
def fetch_uniprot_sequence(uniprot_ac):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        sequence = response.text.split('\n', 1)[1].replace('\n', '').replace('\r', '')
        return sequence
    else:
        st.error(f"Failed to fetch sequence for UniProt ID: {uniprot_ac}")
        return None

# Function to fetch transmembrane helices from UniProt flat file
def fetch_uniprot_tm_helices(uniprot_ac):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}.txt"
    response = requests.get(url)
    tm_helices = []
    out_seq = []  # Extracellular
    in_seq = []   # Cytoplasmic
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

# Function to reformat UniProt annotation to display on structures
def convert_tm_helices_to_pred(sequence, tm_helices, out_seq, in_seq):
    pred_uniprot = ["g"] * len(sequence)  

    for start, end in tm_helices:
        for i in range(start - 1, end):  # UniProt indices are 1-based
            pred_uniprot[i] = "M"
    for start, end in in_seq:
        for i in range(start - 1, end):
            pred_uniprot[i] = "I"
    for start, end in out_seq:
        for i in range(start - 1, end):
            pred_uniprot[i] = "O"
    return "".join(pred_uniprot)

# Function to fetch PDB structure
def fetch_pdb_structure(pdb_code):
    pdb_code_dl = pdb_code.lower()
    url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_code_dl}_updated.cif"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        st.error(f"Failed to fetch structure for PDB code: {pdb_code}")
        return None

# Function to clear old DeepTMHMM results
def clear_old_results():
    # Prevent accidental deletion of demo directory
    if output_dir == demo_dir:
        return
    # Ensure the directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        return  # Exit early if the directory was just created

    # Clear old results
    for filename in os.listdir(output_dir):
        file_path = os.path.join(output_dir, filename)
        try:
            if os.path.isfile(file_path):
                os.remove(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print(f"Error deleting {file_path}: {e}")

# Function to extract TM helices from DeepTMHMM results
def extract_tm_helices(gff3_file):
    tm_helices = []
    ss_tag = None

    if os.path.exists(gff3_file):
        with open(gff3_file, "r") as f:
            for line in f:
                if line.startswith("#"):  # Skip comments
                    continue
                columns = line.strip().split("\t")
                if len(columns) >= 4 and columns[1] == "TMhelix":
                    start = int(columns[2])
                    end = int(columns[3])
                    tm_helices.append((start, end))
                    ss_tag = "helix"
                elif len(columns) >= 4 and columns[1] == "Beta sheet":
                    start = int(columns[2])
                    end = int(columns[3])
                    tm_helices.append((start, end))
                    ss_tag = "beta"
    return tm_helices, ss_tag

# Function to run DeepTMHMM prediction
def run_deeptmhmm_biolib(sequence):
    # Create status containers
    status_container = st.empty()
    progress_bar = st.progress(0)
    
    # Step 1: Clear old results
    status_container.info("🧹 Clearing previous results...")
    progress_bar.progress(10)
    time.sleep(0.5)

    output_dir = os.path.join(script_dir, "biolib_results")
    clear_old_results()

    # Step 2: Write FASTA file
    status_container.info("📝 Writing sequence to FASTA file...")
    progress_bar.progress(20)
    time.sleep(0.5)

    # Write the sequence to a temporary FASTA file
    with open(fasta_file, "w") as f:
        f.write(f">sequence\n{sequence}")
    
    # Step 3: Run DeepTMHMM prediction in the output directory
    status_container.info("🔬 Running DeepTMHMM prediction... This may take a few minutes.")
    progress_bar.progress(30)
    
    command = f"biolib run DTU/DeepTMHMM --fasta {fasta_file}"
    
    try:
        result = subprocess.run(command, shell=True, capture_output=True, text=True, cwd=script_dir)
        
        # Step 4: Processing results
        status_container.info("📊 Processing prediction results...")
        progress_bar.progress(80)
        time.sleep(0.5)
        
        if result.returncode == 0:
            gff3_path = os.path.join(output_dir, "TMRs.gff3")
            tm_helices, ss_tag = extract_tm_helices(gff3_path)
            
            # Step 5: Complete
            status_container.success("✅ DeepTMHMM prediction completed successfully!")
            progress_bar.progress(100)
            time.sleep(1)
            
            # Clear status indicators after success
            status_container.empty()
            progress_bar.empty()
            
            return tm_helices, output_dir
        else:
            status_container.error(f"❌ Error running DeepTMHMM: {result.stderr}")
            progress_bar.empty()
            return None, None
    except Exception as e:
        status_container.error(f"❌ Unexpected error: {e}")
        progress_bar.empty()
        return None, None

# Function to read demo DeepTMHMM prediction for default_unp
def read_demo_results():
    gff3_path = os.path.join(demo_dir, "TMRs.gff3")
    tm_helices, ss_tag = extract_tm_helices(gff3_path)
    pred = ""
    output_dir_demo = demo_dir

    file_path = os.path.join(demo_dir, "predicted_topologies.3line")
    if os.path.exists(file_path):
        with open(file_path, "r") as f:
            lines = f.readlines()
            # Extract the second line after the sequence line for TM prediction
            if len(lines) > 2:
                pred = lines[2].strip()
    return tm_helices, output_dir_demo, pred

# Replaces stmol.showmol
def showmol(view, height=500, width=500):
    html = view._make_html()
    components.html(html, height=height, width=width)

# Function to display the structure visualization
def viewpdb(structure, pred, sequence, af2_tag):
    view = py3Dmol.view(height=360, width=400)
    view.addModel(structure, 'mmcif')
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
        elif res_type == 'B':
            atom_color[nr] = '#830592'
        elif res_type == 'P':
            atom_color[nr] = 'green'
        elif res_type == 'S':
            atom_color[nr] = 'orange'
        elif res_type == 'g':
            atom_color[nr] = 'lightgrey'
        else:
            atom_color[nr] = '#008c74'

    try:
        chain_info, _ = pdb_chainID.get_chain_ids(pdb_code, uniprot_ac)
        chain_ids_list = [entry['Chain ID'] for entry in chain_info if 'Chain ID' in entry]
        chain_ids = ', '.join(chain_ids_list)        
    except Exception as e:
        chain_ids = []
        st.error(f"Failed to get chain IDs: {e}")

        if not chain_ids:
            chain_ids = ['A']
    except Exception as e:
        st.error(f"Error fetching chain IDs: {e}")
        chain_ids = ['A']

    if 'chain_ids' not in st.session_state:
        st.session_state.chain_ids = chain_ids
    else:
        st.session_state.chain_ids = chain_ids

    if af2_tag == 1:
        chain_ids = ['A']
    
    for chain_id in chain_ids:
        view.setStyle({'model': -1, 'chain': chain_id}, {
            'cartoon': {
                'thickness': 0.5,
                'colorscheme': {'prop': 'resi', 'map': atom_color}
            }
        })

    view.zoomTo()
    showmol(view, height=360, width=400)
    chain_ids = st.session_state.chain_ids

# Read the TM prediction file (predicted_topologies.3line)
def get_pred_from_file():
    pred = ""
    file_path = os.path.join(output_dir, "predicted_topologies.3line")
    if os.path.exists(file_path):
        with open(file_path, "r") as f:
            lines = f.readlines()
            # Extract the second line after the sequence line for TM prediction
            if len(lines) > 2:
                pred = lines[2].strip()
    return pred
    
col1, col2, col3 = st.columns([1, 1, 2])
with col1:
    st.markdown(f'UniProt entry: <a href="https://rest.uniprot.org/uniprotkb/{uniprot_ac}.txt" target="_blank">{uniprot_ac}</a>', unsafe_allow_html=True)
with col2:
    pdb_code = pdb_code.upper()
    st.markdown(f'PDB entry: <a href="https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_code}" target="_blank">{pdb_code}</a>', unsafe_allow_html=True)

# Initialize variables
alphafold_structure = None
tm_helices_uniprot = None
tm_helices_pred = None
sequence = None
pdb_structure = None
af2_structure = None
af2_pred = None

## Use demo_results for default_unp and default_pdb
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

# Fetch PDB structure
if fetch_pdb_button:
    st.session_state.pdb_code = pdb_code
    pdb_structure = fetch_pdb_structure(pdb_code)
    if pdb_structure:
        st.session_state.pdb_structure = pdb_structure

# Fetch and process UniProt data
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
        
st.markdown("<a name='pdb-uniprot'></a>", unsafe_allow_html=True)
col1, col2 = st.columns([3, 2])

with col1:
    st.markdown("##### &nbsp;&nbsp;&nbsp;Structure")
with col2:
    st.markdown("##### UniProt TM Annotation")
    
col1, col2 = st.columns([3, 2])

# Display PDB structure with UniProt TM annotations (output)
with col1:
    if pdb_code == default_pdb:
        pdb_structure = fetch_pdb_structure(default_pdb)
    else:
        pdb_structure = st.session_state.get("pdb_structure", None)
        
    if sequence == None:
        sequence = fetch_uniprot_sequence(default_unp)
    else:
        sequence = st.session_state.get("sequence", None)
    
    tm_helices_uniprot = st.session_state.get("tm_helices_uniprot", None)
    pred_uniprot = st.session_state.get("pred_uniprot", None)
    af2_tag = 0
    viewpdb(pdb_structure, pred_uniprot, sequence, af2_tag)
    chain_ids = st.session_state.chain_ids
    st.caption(f"PDB:{pdb_code} ({chain_ids}) with UniProt annotations for {uniprot_ac}.")

# Display UniProt TM annotations (output)
with col2:
    output_str = ""
    tm_helices_uniprot = st.session_state.get("tm_helices_uniprot", [])
    
    if tm_helices_uniprot:
        st.markdown("")
        for start, end in tm_helices_uniprot:
            output_str += f"FT TRANSMEM  {start}  {end}\n"
        st.markdown(f"```\n{output_str}\n```")
        st.markdown(
    '<div style="display: flex; padding: 8px; width: 100%;">'
    '<span style="background-color: white; padding: 2px 5px; text-align: center; font-size: 14px;">Key: </span>'
    '<span style="background-color: powderblue; padding: 2px 5px; text-align: center; font-size: 14px;">Outside</span>'
    '<span style="background-color: pink; padding: 2px 5px; text-align: center; font-size: 14px;">Inside</span>'
    '</div>', unsafe_allow_html=True)
    else:
        st.warning("No TM annotation found in UniProt.")
        
# PDB with TM prediction
# -----------------------------------------------------------------------------

# Run DeepTMHMM prediction and display on PDB structure
# Use demo_results for default_unp

if uniprot_ac == default_unp:
    tm_helices_pred, output_dir, pred = read_demo_results()

if 'uniprot_ac' not in st.session_state:
    st.session_state.uniprot_ac = default_unp
st.markdown("")
st.markdown("")
st.markdown("#### Visualise transmembrane *prediction* on a protein structure.")
col1, col2, col3 = st.columns([1, 2, 1])
with col1:
    fetch_pred_button = st.button("Run TM prediction")
if 'pdb_code' not in st.session_state:
    st.session_state.pdb_code = default_pdb
with col2:
    display_pdb_button = st.button("Show structure with TM prediction")

# Display PDB structure (button)
if display_pdb_button:
    st.session_state.pdb_code = pdb_code
    pdb_structure = fetch_pdb_structure(pdb_code)
    if pdb_structure:
        st.session_state.pdb_structure = pdb_structure

if fetch_pred_button:
    # Fetch sequence and structure
    sequence = fetch_uniprot_sequence(uniprot_ac)
    if sequence is None:
        st.stop()
    tm_helices_pred, output_dir = run_deeptmhmm_biolib(sequence)
    if tm_helices_pred is None or output_dir is None:
        st.warning("DeepTMHMM Prediction has failed.")
    pred = get_pred_from_file()

st.markdown("<a name='pdb-tmhmm'></a>", unsafe_allow_html=True)
col1, col2 = st.columns([3, 2])
with col1:
    st.markdown("##### &nbsp;&nbsp;&nbsp;Structure")
with col2:
    st.markdown("##### DeepTMHMM Prediction")

col1, col2 = st.columns([3, 2])

# Display PDB structure with DeepTMHMM prediction (output)
with col1:
    if pdb_code == default_pdb:
        pdb_structure = fetch_pdb_structure(default_pdb)
    else:
        pdb_structure = st.session_state.get("pdb_structure", None)
        
    if sequence == None:
        sequence = fetch_uniprot_sequence(default_unp)
    else:
        sequence = st.session_state.get("sequence", None)
        pred = get_pred_from_file()
        
    viewpdb(pdb_structure, pred, sequence, af2_tag=0)
    chain_ids = st.session_state.chain_ids
    st.caption(f"PDB:{pdb_code} ({chain_ids}) with DeepTMHMM predictions.")

with col2:
    output_str = ""
    if uniprot_ac == default_unp:
        gff3_path = os.path.join(demo_dir, "TMRs.gff3")
    else:
        gff3_path = os.path.join(output_dir, "TMRs.gff3")
    if os.path.exists(gff3_path):
        st.markdown("")
        tm_helices, ss_tag = extract_tm_helices(gff3_path)
        if ss_tag == None:
            st.error(f"No TM prediction available.")
            print (f"Error in tm_helices: {tm_helices}")
        else:
            for start, end in tm_helices:
                if ss_tag == "helix":
                    output_str += f"FT TRANSMEM  {start}  {end} Helical\n"
                elif ss_tag == "beta":
                    output_str += f"FT TRANSMEM  {start}  {end} Beta stranded\n"

        st.markdown(f"```\n{output_str}\n```")
        st.markdown(
    '<div style="display: flex; padding: 8px; width: 100%;">'
    '<span style="background-color: white; padding: 2px 5px; text-align: center; font-size: 14px;">Key: </span>'
    '<span style="background-color: powderblue; padding: 2px 5px; text-align: center; font-size: 14px;">Outside</span>'
    '<span style="background-color: pink; padding: 2px 5px; text-align: center; font-size: 14px;">Inside</span>'
    '</div>', unsafe_allow_html=True)

# DeepTMHMM plot
# -----------------------------------------------------------------------------
col1, col2 = st.columns([3, 2])
with col1:        
    if uniprot_ac == default_unp:
        plot_path = os.path.join(demo_dir, "plot.png")
    else:
        plot_path = os.path.join(output_dir, "plot.png")
    st.markdown("<a name='tmhmm_plot'></a>", unsafe_allow_html=True)
    if os.path.exists(plot_path):
        st.markdown("##### DeepTMHMM Plot")
        img = mpimg.imread(plot_path)
        plt.imshow(img)
        plt.axis('off')
        st.pyplot(plt)
    else:
        st.warning("TM prediction not available.")

# Acknowledgement
# -----------------------------------------------------------------------------

st.markdown(
    """
    ---
    The transmembrane topology prediction was generated using [DeepTMHMM](https://dtu.biolib.com/DeepTMHMM).
    """
)