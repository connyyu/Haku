import csv
import re
import os
import glob
import requests
import pandas as pd
import streamlit as st
from collections import defaultdict

# Sidebar, title, parameters
# -----------------------------------------------------------------------------
st.set_page_config(page_title="Haku - find UniProt", page_icon="💮")
st.sidebar.title("Navigation")
st.sidebar.markdown("📝 [New query](#top_title)")
st.sidebar.markdown("[Proteins and ligands](#simplified-data)")
st.sidebar.markdown("[UniProt AC list](#uniprot-ac)")
st.sidebar.markdown("[PyMOL Fetch Command](#pymol-fetch-command)")
st.sidebar.markdown("[PyMOL Colour Command](#pymol-colour-command)")

st.markdown("<a name='top_title'></a>", unsafe_allow_html=True)
st.markdown("#### List all the proteins (UniProt AC) in the structures.")

script_dir = os.path.dirname(os.path.abspath(__file__)) # location of pages directory
output_dir = os.path.join(script_dir, "chain_info")

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
        if line.startswith('_reflns.d_resolution_high'):
            parts = line.split()
            if len(parts) > 1:
                resolution = parts[1]
                method = "Xtal"
                break
        elif line.startswith('_em_3d_reconstruction.resolution'):
            parts = line.split()
            if len(parts) > 1:
                resolution = parts[1]
                method = "EM"
                break
        elif line.startswith("_exptl.method 'Solution NMR'"):
            method = "NMR"
            break
        elif line.startswith('_em_3d_reconstruction.symmetry_type'):
            symmetry_type_found = True
        elif symmetry_type_found:  # Look at the next line after symmetry_type
            parts = line.split()
            if len(parts) >= 7:
                resolution = parts[6]  # Take the 7th column as resolution
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
    uniprot_ids = []
    chain_mapping = {}
    seen_chain_ids = {}  # To track which chains are associated with each UniProt accession
    
    file_path = os.path.join(output_dir, cif_file)
    with open(file_path, 'r') as f:
        lines = f.readlines()

    i = 0
    process_sifts = False
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('_struct_ref.pdbx_db_accession'):
            found_valid_line = False
            j = i + 1
            while j < len(lines):
                if 'UNP' in lines[j]:  # Check if 'UNP' is in the line
                    parts = lines[j].strip().split()
                    if 'UNP' in parts and parts[1] == 'UNP':  # Ensure 'UNP' is in the expected position
                        for part in parts:
                            if len(part) == 6:
                                uniprot_accession = part  # Get UniProt accession
                                uniprot_id = parts[2]  # Get UniProt ID
                                entity_id = parts[0]  # Get entity_id
                                if entity_id not in uniprot_mapping:
                                    uniprot_mapping[entity_id] = uniprot_accession
                                    uniprot_ids.append((entity_id, uniprot_id))
                                    found_valid_line = True
                j += 1
            if not found_valid_line:
                for k in range(len(lines)):
                    if lines[k].startswith('_struct_ref.entity_id'):
                        entity_id = lines[k].strip().split()[-1]
                    if lines[k].startswith('_struct_ref.pdbx_db_accession'):
                        uniprot_accession = lines[k].strip().split()[-1]
                    if lines[k].startswith('_struct_ref.db_code'):
                        uniprot_id = lines[k].strip().split()[-1]
                if entity_id and uniprot_accession:
                    if entity_id not in uniprot_mapping:
                        uniprot_mapping[entity_id] = uniprot_accession
                        uniprot_ids.append((entity_id, uniprot_id))
        elif line.startswith('_struct_ref_seq.align_id'):
            i += 1  # Skip the header line
            while i < len(lines) and not lines[i].strip().startswith('loop_'):
                parts = lines[i].strip().split()
                if len(parts) > 3:
                    entity_id = parts[1]
                    chain_id = parts[3]
                    if entity_id in uniprot_mapping:
                        if entity_id not in chain_mapping:
                            chain_mapping[entity_id] = set()  # Initialize as a set to avoid duplicates
                        chain_mapping[entity_id].add(chain_id)  # Add the chain ID to the set of chains for this entity ID
                i += 1
            continue
        # Update the UniProt ID and AC with the new mappings in SIFTS
        elif line.startswith('_pdbx_sifts_unp_segments.identity'):
            process_sifts = True
        elif process_sifts and line.startswith('loop'):
            process_sifts = False
        elif process_sifts:
            parts = line.strip().split()
            if len(parts) > 3:
                entity_id = parts[0]
                uniprot_accession = parts[2]
                if len(uniprot_accession) == 6 and entity_id in uniprot_mapping:
                    uniprot_mapping[entity_id] = uniprot_accession
                    new_uniprot_id = fetch_uniprot_id(uniprot_accession)
                    if new_uniprot_id:
                        # Update existing uniprot_id for the entity_id
                        for idx, (eid, uid) in enumerate(uniprot_ids):
                            if eid == entity_id:
                                uniprot_ids[idx] = (entity_id, new_uniprot_id)
                                break  # Stop after finding the match

        i += 1
    
    chain_mapping = {entity_id: sorted(list(chains)) for entity_id, chains in chain_mapping.items()}
    
    return uniprot_mapping, uniprot_ids, chain_mapping

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
        chain_id = row['Chain ID']
        uniprot_accession = row['UniProt AC']
        chain_uniprot_map[chain_id].add(uniprot_accession)
    warnings = []
    for chain_id, uniprot_accessions in chain_uniprot_map.items():
        if len(uniprot_accessions) > 1:
            warnings.append(f"{pdb_code} chain {chain_id} is associated with mutiple UniProt AC: {', '.join(uniprot_accessions)}.")
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
                'Chain ID': chain_id,
                'UniProt AC': uniprot_accession,
                'UniProt ID': uniprot_id,
                'HETATM': hetatm_names,
                'Resolution': resolution,
                'Method': method
            })

    # Check for warnings and print them
    warnings = check_chimera(table_data)
    if warnings:
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
            chain_id = row['Chain ID']
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
            chain_id = row['Chain ID']
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
            'Chain ID': chain_ids_str,
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
        fieldnames = ['PDB', 'Entity', 'Chain ID', 'UniProt AC', 'UniProt ID', 'HETATM', 'Resolution', 'Method']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(simplified_data)
    
    return simplified_data

# Function to alternate between two colors
import pandas as pd

# Function to alternate between two colors
import pandas as pd

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

# Section 1: Download and parse cif files for display
# -----------------------------------------------------------------------------
progress_bar = st.progress(0)
status_text = st.empty()
all_data = []
pdb_codes = []

# Input form
with st.form(key='form'):
    col1, col2 = st.columns([3, 1])
    with col1:
        pdb_input = st.text_input("Enter PDB codes (comma or space-separated):",
                                  "9eii 9eij 9eih").strip()
    with col2:
        st.markdown("<h5 style='font-size: 10px;'></h5>", unsafe_allow_html=True)
        fetch_button = st.form_submit_button("Analyse structures")

if fetch_button:
    pdb_codes = re.split(r'[,\s]+', pdb_input.strip().upper())
    pdb_codes = [code.strip() for code in pdb_codes if code]
    total_pdbs = len(pdb_codes)
    all_data = []

    for i, pdb_code in enumerate(pdb_codes, start=1):
        status_text.text(f"Analysing {pdb_code} ({i}/{total_pdbs})...")
        table_data = process_cif_file(pdb_code)
        if table_data:
            all_data.extend(table_data)
        progress_bar.progress(i / total_pdbs)

    status_text.text("Parsing data for output...")

# Save PDB codes to file
pdb_code = st.session_state.get("pdb_codes", None)
if pdb_code:
    file_path = os.path.join(output_dir, "pdb_codes.txt")
    with open(file_path, 'w') as f:
        f.write('\n'.join(pdb_codes))
    for pdb_code in pdb_codes:
        table_data = process_cif_file(pdb_code)
        if table_data:
            all_data.extend(table_data)

status_text.empty()

# Save consolidated table to CSV file and display data
if all_data:
    all_data = [dict(t) for t in {tuple(d.items()) for d in all_data}]
    file_path = os.path.join(output_dir, "parsed_cif.csv")
    with open(file_path, 'w', newline='') as csvfile:
        fieldnames = ['PDB', 'Entity', 'Chain ID', 'UniProt AC', 'UniProt ID', 'HETATM', 'Resolution', 'Method']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_data)
    sorted_data = sort_parsed_cif('parsed_cif.csv')

    simplified_data = print_simplified_output('parsed_cif.csv')
    df = pd.DataFrame(simplified_data)
    styled_df = alternate_colors(df)
    st.markdown("<a name='simplified-data'></a>", unsafe_allow_html=True)
    st.markdown("#### Proteins and ligands:")
    st.dataframe(styled_df, hide_index=True, use_container_width=True)

# Section 2: List of UniProt AC
# -----------------------------------------------------------------------------
    st.markdown("<a name='uniprot-ac'></a>", unsafe_allow_html=True)
    st.markdown("#### UniProt AC in these entries:")
    unique_uniprot_acs = sorted(set(row['UniProt AC'] for row in all_data if row['UniProt AC']))
    st.code(" ".join(unique_uniprot_acs))

# Section 3: PyMOL Fetch Command
# -----------------------------------------------------------------------------
    st.markdown("<a name='pymol-fetch-command'></a>", unsafe_allow_html=True)
    st.markdown("#### PyMOL Fetch Command:")
    st.code("fetch " + "; fetch ".join(str(pdb_code).lower() for pdb_code in pdb_codes))

# Section 4: PyMOL Colouring Command
# -----------------------------------------------------------------------------
    st.markdown("<a name='pymol-colour-command'></a>", unsafe_allow_html=True)
    st.markdown("#### PyMOL Colouring Command:")
    pymol_command(os.path.join(output_dir, 'parsed_cif.csv'))
    pdb_pymol_file_path = os.path.join(output_dir, 'pdb_pymol.txt')
    with open(pdb_pymol_file_path, 'r') as file:
        pymol_commands = file.read()
        st.code(pymol_commands)

progress_bar.empty()

## Delete cif file at the end of run
for file in glob.glob(os.path.join(output_dir, "*.cif")):
        os.remove(file)


