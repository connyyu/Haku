import csv
import re
import os
import glob
import requests
import sys
from pathlib import Path
import pandas as pd
import streamlit as st
from collections import defaultdict
import py3Dmol
import streamlit.components.v1 as components

# Sidebar, title, parameters
# -----------------------------------------------------------------------------
st.set_page_config(page_title="Haku - Structure viewer", page_icon="💮")

# Initialize session state variables
if 'all_data' not in st.session_state:
    st.session_state.all_data = []
if 'pdb_codes' not in st.session_state:
    st.session_state.pdb_codes = []

default_pdb = '9eii 9eij 9eih'

with st.sidebar:
    st.header("📝 New query")
    pdb_input = st.text_area("Enter PDB codes (comma or space-separated):", default_pdb).strip()
    fetch_button = st.button("Analyse structure")
    if fetch_button:
        st.session_state.guide_viewer = False
    st.sidebar.markdown("[Proteins](#simplified-data)")
    st.sidebar.markdown("[Structures](#structure)")
    st.sidebar.markdown("[PyMOL commands](#pymol-commands)")
    guide_viewer = st.checkbox("Instructions", value=st.session_state.get("guide_viewer", True))
    st.session_state.guide_viewer = guide_viewer

st.markdown("<a name='top_title'></a>", unsafe_allow_html=True)
st.markdown("#### Navigate protein structures in a 3D viewer.")

# Get the directory where this script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
# Add the scripts subdirectory to Python path
scripts_dir = os.path.join(script_dir, "scripts")
if scripts_dir not in sys.path:
    sys.path.insert(0, scripts_dir)
output_dir = os.path.join(script_dir, "chain_info")

# Import the pdb_unpID function from the scripts folder
try:
    import pdb_chainID
except ImportError as e:
    st.error(f"Failed to import pdb_chainID module: {e}")

# Instructions
# -----------------------------------------------------------------------------
if st.session_state.get('guide_viewer'):
    st.info("""
    Input:&emsp;**one or multiple PDB codes**  
      
    Output:
    * all protein components mapped in UniProt
    * resolution and method of each structure
    * **selected protein shown in a 3D viewer**
    * PyMOL command to colour and label each protein in the structures
    
    """)

# Define functions
# -----------------------------------------------------------------------------
def download_cif_file(pdb_code):
    pdb_code_dl = pdb_code.lower()
    url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_code_dl}_updated.cif"
    response = requests.get(url)
    if response.status_code == 200:
        file_path = os.path.join(output_dir, f"{pdb_code_dl}.cif")
        with open(file_path, 'wb') as f:
            f.write(response.content)        
        return f"{pdb_code_dl}.cif"
    else:
        print(f"\033[1mFailed to download CIF file for {pdb_code_dl}\033[0m")
        return None

def parse_resolution(cif_file):
    resolution = None
    method = ""
    symmetry_type_found = False
    
    file_path = os.path.join(output_dir, cif_file)
    with open(file_path, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        stripped_line = line.strip()
        if not stripped_line:
            continue
        tokens = stripped_line.split()
        if tokens[0] == '_reflns.d_resolution_high':
            if len(tokens) > 1:
                resolution = tokens[1]
                method = "Xtal"
                break
        elif tokens[0] == '_em_3d_reconstruction.resolution':
            if len(tokens) > 1:
                resolution = tokens[1]
                method = "EM"
                continue  # Still check symmetry_type
        elif stripped_line == "_exptl.method 'Solution NMR'":
            method = "NMR"
            break

        elif tokens[0] == '_em_3d_reconstruction.symmetry_type':
            symmetry_type_found = True
            continue
        elif symmetry_type_found:  # Look at the next line after symmetry_type, if resolution wasn't already found
            if len(tokens) >= 7 and resolution is None:
                resolution = tokens[6] #7th column
                method = "EM"
            break
    # 2 decimal places for resolution
    if resolution is not None:
        try:
            resolution = f"{float(resolution):.2f}"
        except ValueError:
            resolution = None

    return resolution, method

def parse_hetatms(cif_file):
    hetatms = set()
    file_path = os.path.join(output_dir, cif_file)
    with open(file_path, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('HETATM'):
            parts = line.split()
            hetatm_name = parts[5]  # HETATM name (e.g., GLC, ZN)
            hetatms.add(hetatm_name)  # Use a set to avoid duplicates
    return sorted(hetatms)  # Return a sorted list

def fetch_uniprot_id(uniprot_ac):
    url = f'https://rest.uniprot.org/uniprotkb/{uniprot_ac}.json'
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        return data.get('uniProtkbId', None)
    return None

def parse_uniprot_mapping(cif_file):
    uniprot_mapping = {}
    uniprot_ids = {}
    chain_mapping = {}
    ref_id_to_entity_id = {}

    fallback_accession = None
    fallback_ID = None
    fallback_chain = None

    ## for streamlit app ##
    file_path = os.path.join(output_dir, cif_file)
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Scan for fallbacks
    for i, line in enumerate(lines):
        line = line.strip()
        if line.startswith('_struct_ref.pdbx_db_accession'):
            fallback_accession = (line.split()[1] if len(line.split()) > 1
                                  else lines[i + 1].strip().split()[0])
        elif line.startswith('_struct_ref.db_code'):
            fallback_ID = (line.split()[1] if len(line.split()) > 1
                           else lines[i + 1].strip().split()[0])
        elif line.startswith('_struct_ref_seq.pdbx_strand_id'):
            parts = line.split()
            if len(parts) > 1:
                fallback_chain_candidate = parts[1]
                if len(fallback_chain_candidate) == 1 and fallback_chain_candidate.isalnum():
                    fallback_chain = fallback_chain_candidate
            else:
                k = i + 1
                while k < len(lines):
                    next_line = lines[k].strip()
                    if not next_line or next_line.startswith('_'):
                        break
                    fallback_chain_candidate = next_line.split()[0]
                    if len(fallback_chain_candidate) == 1 and fallback_chain_candidate.isalnum():
                        fallback_chain = fallback_chain_candidate
                        break
                    k += 1

    i = 0
    process_sifts = False
    sifts_mapping = {}
    sifts_chain_mapping = {}
    struct_ref_count = 0
    struct_ref_seq_count = 0
    sifts_count = 0

    while i < len(lines):
        line = lines[i].strip()

        if line.startswith('loop_'):
            struct_ref_headers = []
            header_line_idx = i + 1
            while header_line_idx < len(lines):
                header_line = lines[header_line_idx].strip()
                if header_line.startswith('_struct_ref.'):
                    struct_ref_headers.append(header_line)
                    header_line_idx += 1
                elif header_line.startswith('_') and not header_line.startswith('_struct_ref.'):
                    break
                else:
                    break

            if struct_ref_headers:
                db_name_idx = db_code_idx = entity_id_idx = accession_idx = ref_id_idx = None

                for idx, header in enumerate(struct_ref_headers):
                    if header == '_struct_ref.id':
                        ref_id_idx = idx
                    elif header == '_struct_ref.db_name':
                        db_name_idx = idx
                    elif header == '_struct_ref.db_code':
                        db_code_idx = idx
                    elif header == '_struct_ref.entity_id':
                        entity_id_idx = idx
                    elif header == '_struct_ref.pdbx_db_accession':
                        accession_idx = idx

                data_line_idx = header_line_idx
                while data_line_idx < len(lines):
                    data_line = lines[data_line_idx].strip()
                    if not data_line or data_line.startswith('#'):
                        data_line_idx += 1
                        continue
                    if data_line.startswith('loop_') or data_line.startswith('_'):
                        break

                    parts = data_line.split()
                    struct_ref_count += 1

                    if (None not in (ref_id_idx, db_name_idx, db_code_idx, entity_id_idx, accession_idx) and
                        len(parts) > max(ref_id_idx, db_name_idx, db_code_idx, entity_id_idx, accession_idx)):
                        
                        if parts[db_name_idx] == 'UNP':
                            ref_id = parts[ref_id_idx]
                            entity_id = parts[entity_id_idx]
                            uniprot_id = parts[db_code_idx]
                            uniprot_accession = parts[accession_idx]

                            ref_id_to_entity_id[ref_id] = entity_id

                            if len(uniprot_accession) >= 6 and '-' not in uniprot_accession and '_' not in uniprot_accession:
                                uniprot_mapping[entity_id] = uniprot_accession
                                uniprot_ids[entity_id] = uniprot_id  # Store as dict for easier lookup

                    data_line_idx += 1
                i = data_line_idx
                continue

        elif line.startswith('_struct_ref_seq.align_id'):
            i += 1
            while i < len(lines):
                next_line = lines[i].strip()
                if not next_line or next_line.startswith(('loop_', '#')):
                    break
                parts = next_line.split()
                if len(parts) > 8:
                    struct_ref_seq_count += 1
                    ref_id = parts[1]
                    chain_id = parts[3]
                    uniprot_accession = parts[8]

                    entity_id = ref_id_to_entity_id.get(ref_id)
                    if entity_id is None:
                        entity_id = ref_id

                    if entity_id:
                        chain_mapping.setdefault(entity_id, set()).add(chain_id)
                        if (entity_id not in uniprot_mapping and
                            len(uniprot_accession) >= 6 and 
                            '-' not in uniprot_accession and 
                            '_' not in uniprot_accession):
                            uniprot_mapping[entity_id] = uniprot_accession
                            if entity_id not in uniprot_ids:
                                uniprot_ids[entity_id] = None
                i += 1
            continue

        elif line.startswith('_pdbx_sifts_unp_segments.identity'):
            process_sifts = True

        elif process_sifts:
            if line.startswith(('loop_', '#')):
                process_sifts = False
            elif line:
                parts = line.split()
                if len(parts) > 3:
                    sifts_count += 1
                    entity_id, chain_id, sifts_accession = parts[0], parts[1], parts[2]
                    if len(sifts_accession) >= 6 and '-' not in sifts_accession:
                        sifts_mapping[entity_id] = sifts_accession
                        sifts_chain_mapping.setdefault(entity_id, set()).add(chain_id)

        i += 1

    # Use fallback if no mapping found
    if not uniprot_mapping and not sifts_mapping:
        fallback_entity = "1"

        if fallback_accession:
            uniprot_mapping[fallback_entity] = fallback_accession
            uniprot_ids[fallback_entity] = fallback_ID

            chains_for_entity = chain_mapping.get(fallback_entity)
            if chains_for_entity:
                chain_mapping[fallback_entity] = chains_for_entity
            else:
                chain_mapping[fallback_entity] = set([fallback_chain] if fallback_chain else [])

    # Final reconciliation - Use SIFTS chains only if struct_ref_seq chains are NOT available
    final_uniprot_mapping = {}
    final_uniprot_ids = []
    final_chain_mapping = {}

    all_entities = set(chain_mapping.keys()) | set(uniprot_mapping.keys()) | set(sifts_mapping.keys()) | set(sifts_chain_mapping.keys())
    
    for entity in all_entities:
        original_ac = uniprot_mapping.get(entity, None)
        uniprot_id = uniprot_ids.get(entity, None)
        sifts_ac = sifts_mapping.get(entity, None)

        final_ac = sifts_ac if sifts_ac else original_ac

        chains_from_struct_ref_seq = list(chain_mapping.get(entity, set()))
        chains_from_sifts = list(sifts_chain_mapping.get(entity, set()))
        
        if chains_from_struct_ref_seq:
            final_chains = sorted(chains_from_struct_ref_seq)
        elif chains_from_sifts:
            final_chains = sorted(chains_from_sifts)
        else:
            final_chains = []
        
        if final_ac:
            final_uniprot_mapping[entity] = final_ac

            uid = uniprot_ids.get(entity)
            if uid is None and fallback_ID:
                uid = fallback_ID
        
            final_uniprot_ids.append((entity, uid))
            final_chain_mapping[entity] = final_chains

    for entity, chains in final_chain_mapping.items():
        uid = None
        for ent, uniprot_id in final_uniprot_ids:
            if ent == entity:
                uid = uniprot_id
                break
        if fallback_chain and uid == fallback_ID:
            if chains and (fallback_chain not in chains):
                final_chain_mapping[entity] = [fallback_chain]

    return final_uniprot_mapping, final_uniprot_ids, final_chain_mapping

def sort_parsed_cif(csv_file):
    # Read the parsed CSV file and store the rows
    file_path = os.path.join(output_dir, csv_file)
    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = list(reader)
    
    # Sort the rows by 'PDB' and 'Entity'
    rows.sort(key=lambda x: (x['PDB'], x['Entity']))

    # Save the sorted data back to parsed_cif.csv
    with open(file_path, 'w', newline='') as csvfile:
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    return rows
    
def check_chimera(table_data):
    chain_uniprot_map = defaultdict(set)
    for row in table_data:
        pdb_code = row['PDB'].upper()
        chain_id = row['Chain']
        uniprot_accession = row['UniProt AC']
        chain_uniprot_map[chain_id].add(uniprot_accession)
    warnings = []
    for chain_id, uniprot_accessions in chain_uniprot_map.items():
        if len(uniprot_accessions) > 1:
            warnings.append(f"{pdb_code} chain {chain_id} is associated with mutiple UniProt ACs: {', '.join(uniprot_accessions)}.")
    return warnings
    
def process_cif_file(pdb_code):
    cif_file = download_cif_file(pdb_code)
    if not cif_file:
        print(f"\033[1mFailed to download CIF file for {pdb_code}. Skipping...\033[0m")
        return

    resolution, method = parse_resolution(cif_file)    
    hetatms = parse_hetatms(cif_file)
    uniprot_mapping, uniprot_ids, chain_mapping = parse_uniprot_mapping(cif_file)

    hetatm_names = ', '.join(hetatms)
    table_data = []
    
    printed_chains = set()
    for entity_id, uniprot_accession in uniprot_mapping.items():
        uniprot_id = next((id for ent_id, id in uniprot_ids if ent_id == entity_id), None)
        chain_ids = chain_mapping.get(entity_id, "A")
        for chain_id in chain_ids:
            # Create a more complex key to check duplicates (pdb_code, uniprot_accession, chain_id, entity_id)
            unique_key = (pdb_code.lower(), uniprot_accession, chain_id, entity_id)
            if unique_key in printed_chains:
                continue
            printed_chains.add(unique_key)
            table_data.append({
                'PDB': pdb_code.lower(),
                'Entity': entity_id,
                'Chain': chain_id,
                'UniProt AC': uniprot_accession,
                'UniProt ID': uniprot_id,
                'HETATM': hetatm_names,
                'Resolution': resolution,
                'Method': method
            })

    # Check for warnings and print them
    warnings = check_chimera(table_data)   # e.g. PDB:7F83
    if warnings:
        st.warning("⚠️ **Chimera warning!!**")
        for warning in warnings:
            st.write(f"🔴 {warning}")
        print("\n\033[1;31mChimera warning!!\033[0m")
        for warning in warnings:
            print(f"\033[1;31m{warning}\033[0m")
    
    return table_data

def pymol_command(csv_file):
    # Output path for the PyMOL commands
    output_file = os.path.join(os.path.dirname(csv_file), "pdb_pymol.txt")
    
    # List of predefined colors
    colors = [
        "carbon", "cyan", "lightmagenta", "yellow", "salmon", "slate", "orange",
        "deepteal", "violetpurple", "hydrogen", "marine", "olive", "smudge",
        "teal", "wheat", "lightpink", "skyblue"
    ]    
    uniprot_chains = defaultdict(lambda: defaultdict(list))

    with open(csv_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            pdb_id = row['PDB']
            chain_id = row['Chain']
            uniprot_id = row['UniProt ID']
            # Add the chain to the respective UniProt ID and PDB ID
            uniprot_chains[uniprot_id][pdb_id].append(chain_id)
    
    # Open the output file for writing
    with open(output_file, "w") as outfile:
        outfile.write("util.cbaw;\n")
        color_index = 0  # Start with the first color in the list
        
        reformatted_selections = {}
        for uniprot_id, pdb_dict in uniprot_chains.items():
            selections = []
            for pdb_id, chains in pdb_dict.items():
                chains_str = "+".join(chains)
                selections.append(f"(chain {chains_str} and {pdb_id})")
            
            reformatted_selections[uniprot_id] = " or ".join(selections)
        
        for uniprot_id, selection in reformatted_selections.items():
            color = colors[color_index % len(colors)]
            outfile.write(f"sele {uniprot_id}, {selection}; color {color}, {uniprot_id};\n")
            color_index += 1  # Move to the next color
        
        outfile.write("deselect; util.cnc") 

    return output_file  # Return the generated PyMOL command file path
    
def print_simplified_output(csv_file):
    # defaultdict with list to collect chain IDs
    data = defaultdict(list)
    
    # Read the parsed CSV file and group the chain IDs by UniProt accession and PDB ID
    file_path = os.path.join(output_dir, csv_file)
    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            pdb_id = row['PDB']
            entity_id = row['Entity']
            chain_id = row['Chain']
            uniprot_ac = row['UniProt AC']
            uniprot_id = row['UniProt ID']
            hetatm = row['HETATM']
            resolution = row['Resolution'] 
            method = row['Method']

            # Group chain IDs by PDB ID and UniProt accession
            data[(pdb_id, entity_id, uniprot_ac, uniprot_id, hetatm, resolution, method)].append(chain_id)
    
    # Prepare the simplified output
    simplified_data = []
    for (pdb_id, entity_id, uniprot_ac, uniprot_id, hetatm, resolution, method), chain_ids in data.items():
        # Combine the chain IDs into a single string
        chain_ids_str = '/'.join(sorted(chain_ids))
        simplified_data.append({
            'PDB': pdb_id,
            'Entity': entity_id,
            'Chain': chain_ids_str,
            'UniProt AC': uniprot_ac,
            'UniProt ID': uniprot_id,
            'HETATM': hetatm,
            'Resolution': resolution,
            'Method': method
        })
    
    # Sort the simplified data by PDB and Entity
    simplified_data.sort(key=lambda x: (x['PDB'], x['Entity']))

    # Save simplified data to a new CSV file (parsed_cif_simplified.csv)
    simplified_file_path = os.path.join(output_dir, "parsed_cif_simplified.csv")
    with open(simplified_file_path, 'w', newline='') as csvfile:
        fieldnames = ['PDB', 'Entity', 'Chain', 'UniProt AC', 'UniProt ID', 'HETATM', 'Resolution', 'Method']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(simplified_data)
    
    return simplified_data

def alternate_colors(df):
    color_1 = 'background-color: #e2fffb;'
    color_2 = ''
    row_styles = []
    last_pdb_code = None
    alternate = False

    for idx, row in df.iterrows():
        pdb_code = row['PDB']
        
        if pdb_code != last_pdb_code:
            alternate = not alternate  # Toggle color for a new PDB code
        if alternate:
            row_styles.append([color_1] * len(df.columns))
        else:
            row_styles.append([color_2] * len(df.columns))
        last_pdb_code = pdb_code

    def apply_row_styles(row):
        return row_styles[row.name]
    
    return df.style.apply(apply_row_styles, axis=1)

# Replaces stmol.showmol
def showmol(view, height=500, width=500):
    html = view._make_html()
    components.html(html, height=height, width=width)

# Function to display the structure visualization
def viewpdb(structure, sequence, pdb_code, uniprot_ac, show_all=True):
    view = py3Dmol.view(height=450, width=600)
    view.addModel(structure, 'mmcif')
    view.setBackgroundColor('#eeeeee')
    view.spin(False)

    try:
        chain_info, _ = pdb_chainID.get_chain_ids(pdb_code, uniprot_ac)
        chain_ids_list = [entry['Chain ID'] for entry in chain_info if 'Chain ID' in entry]
        chain_ids = ', '.join(chain_ids_list)      
    except Exception as e:
        chain_ids = []
        st.error(f"Failed to get chain IDs: {e}")

    if not chain_ids:
        chain_ids = ['A']

    if 'chain_ids' not in st.session_state:
        st.session_state.chain_ids = chain_ids
    else:
        st.session_state.chain_ids = chain_ids

    if show_all:
        view.setStyle({'model': 0}, {'cartoon': {'color': 'powderblue'}})
        for chain_id in chain_ids:
            view.setStyle({'chain': chain_id}, {'cartoon': {'color': '#830592'}})
    else:
        for chain_id in chain_ids:
            view.setStyle({'chain': chain_id}, {'cartoon': {'color': '#830592'}})

    view.zoomTo()
    showmol(view, height=450, width=600)

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

def fetch_uniprot_sequence(uniprot_ac):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        sequence = response.text.split('\n', 1)[1].replace('\n', '').replace('\r', '')
        return sequence
    else:
        st.error(f"Failed to fetch sequence for UniProt ID: {uniprot_ac}")
        return None
    
# Section 1: Download and parse cif files for display
# -----------------------------------------------------------------------------
progress_bar = st.progress(0)
status_text = st.empty()
all_data = []
pdb_codes = []

if fetch_button:
    pdb_codes = re.split(r'[,\s]+', pdb_input.strip().upper())
    pdb_codes = [code.strip() for code in pdb_codes if code]
    st.session_state.pdb_codes = pdb_codes
    total_pdbs = len(pdb_codes)
    all_data = []

    for i, pdb_code in enumerate(pdb_codes, start=1):
        status_text.text(f"Analysing {pdb_code} ({i}/{total_pdbs})...")
        table_data = process_cif_file(pdb_code)
        if table_data:
            all_data.extend(table_data)
        progress_bar.progress(i / total_pdbs)

    status_text.text("Parsing data for output...")

status_text.empty()

# Save consolidated table to CSV file
if all_data:
    all_data = [dict(t) for t in {tuple(d.items()) for d in all_data}]
    st.session_state.all_data = all_data
    file_path = os.path.join(output_dir, "parsed_cif.csv")
    with open(file_path, 'w', newline='') as csvfile:
        fieldnames = ['PDB', 'Entity', 'Chain', 'UniProt AC', 'UniProt ID', 'HETATM', 'Resolution', 'Method']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_data)
    st.session_state.csv_path = file_path
    sorted_data = sort_parsed_cif('parsed_cif.csv')

# Display simplified data if all_data exists
if "all_data" in st.session_state and st.session_state.all_data:
    file_path = st.session_state.get("csv_path", "parsed_cif.csv")  # fallback
    sorted_data = sort_parsed_cif(file_path)
    simplified_data = print_simplified_output(file_path)
    df = pd.DataFrame(simplified_data)
    styled_df = alternate_colors(df)
    st.markdown("<a name='simplified-data'></a>", unsafe_allow_html=True)
    st.markdown("#### Proteins")
    st.dataframe(styled_df, hide_index=True, use_container_width=True)

# Section 2: Display the structure
# -----------------------------------------------------------------------------
    
# Initialize variables
sequence = None
pdb_structure = None
all_data = st.session_state.get("all_data", [])
pdb_codes = st.session_state.get("pdb_codes", [])

st.markdown("<a name='structure'></a>", unsafe_allow_html=True)

if all_data and pdb_codes:
    st.markdown("#### Structures")  
    col1, col2 = st.columns([1, 3])

    with col1:
        # Only show UniProt ID from the selected PDB
        pdb_code = st.selectbox("Select a PDB code:", pdb_codes, index=0)
        filtered_rows = [row for row in all_data if row.get('PDB', '').upper() == pdb_code.upper()]
        
        uniprot_id_to_ac = {row['UniProt ID']: row['UniProt AC'] for row in filtered_rows}
        unique_uniprot_ids = sorted(uniprot_id_to_ac.keys())
        
        if unique_uniprot_ids:
            selected_uniprot_id = st.selectbox("UniProt ID: (purple)", unique_uniprot_ids, index=0)
            uniprot_ac = uniprot_id_to_ac[selected_uniprot_id]
            sequence = fetch_uniprot_sequence(uniprot_ac)
            pdb_structure = fetch_pdb_structure(pdb_code)
            show_all = st.checkbox("Show all chains", value=True)
            st.markdown("")
            st.caption("Only UniProt entries can be selected. Ligands and solvents are not shown.")
            
            # Save to session state
            st.session_state.pdb_code = pdb_code
            st.session_state.uniprot_ac = uniprot_ac
            st.session_state.sequence = sequence
            st.session_state.pdb_structure = pdb_structure
            
            with col2:
                viewpdb(pdb_structure, sequence, pdb_code, uniprot_ac, show_all)
                chain_ids = st.session_state.get("chain_ids", "N/A")
                chain_ids_str = None  # Initialize to None
                
                if isinstance(chain_ids, str):
                    filtered_chains = [ch.strip() for ch in chain_ids.split(',') if ch.strip()]
                    chain_ids_str = ', '.join(filtered_chains)
                
                if chain_ids_str:
                    st.caption(f"PDB:{pdb_code.upper()} ({chain_ids_str}) highlighting UniProt entry {uniprot_ac.upper()}.")
                else:
                    st.caption(f"PDB:{pdb_code.upper()} highlighting UniProt entry {uniprot_ac.upper()}. Mapping is incomplete in this entry.")
                
# Section 3: PyMOL Commands
# -----------------------------------------------------------------------------
st.markdown("<a name='pymol-commands'></a>", unsafe_allow_html=True)
if all_data:
    st.markdown("#### PyMOL Commands")
    st.code("fetch " + "; fetch ".join(str(pdb_code).lower() for pdb_code in pdb_codes))

    pymol_command(os.path.join(output_dir, 'parsed_cif.csv'))
    pdb_pymol_file_path = os.path.join(output_dir, 'pdb_pymol.txt')
    with open(pdb_pymol_file_path, 'r') as file:
        pymol_commands = file.read()
        st.code(pymol_commands)

    if st.button("Other PyMOL Commands"):
        st.switch_page("pages/Pymol_commands.py")

# End-of-run commands
# -----------------------------------------------------------------------------
progress_bar.empty()

## Delete cif file after analysis
for file in glob.glob(os.path.join(output_dir, "*.cif")):
        os.remove(file)


