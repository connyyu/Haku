import streamlit as st
import re
import pandas as pd
import subprocess
import time
import asyncio
import aiohttp
import concurrent.futures
import io
import os
import sys
from pathlib import Path
import requests

# Sidebar, title, parameters
# -----------------------------------------------------------------------------
st.set_page_config(page_title="Haku - Residual interactions", page_icon="💮", layout="wide")

st.markdown("<a name='top_title'></a>", unsafe_allow_html=True)
st.markdown("#### Show all the residual interactions of a specific ligand.")

default_pdb = '7UPI 7TYG 7TXH 7TVG 7TVF 7T7A 7SD1 7SD0'
script_dir = os.path.dirname(os.path.abspath(__file__)) # location of pages directory
pdb_dir = os.path.join(script_dir, "scripts")
script_path = os.path.join(pdb_dir, "pdb_unpID.py") # script to determine UniProt ID

# Initialize session state
if 'df_interactions' not in st.session_state:
    st.session_state.df_interactions = None
if 'processing_complete' not in st.session_state:
    st.session_state.processing_complete = False

def get_ligand_name(hetatm_code):
    """Fetch ligand name from PDBe API"""
    try:
        url = f"https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/{hetatm_code}"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            if hetatm_code in data:
                return data[hetatm_code][0].get('name', '')
    except Exception as e:
        st.error(f"Error fetching ligand name for {hetatm_code}: {e}")
    return ""

async def fetch_bound_molecules_async(pdb_code, session):
    url = f"https://www.ebi.ac.uk/pdbe/graph-api/pdb/bound_molecules/{pdb_code}"
    
    try:
        async with session.get(url) as response:
            if response.status != 200:
                return []
            
            data = await response.json()
            if pdb_code not in data:
                return []
            
            records = []
            for entry in data[pdb_code]:
                bm_id = entry.get("bm_id", "Unknown")
                ligands = entry.get("composition", {}).get("ligands", [])
                for ligand in ligands:
                    chem_comp_id = ligand.get("chem_comp_id", "Unknown")
                    records.append((pdb_code.upper(), bm_id, chem_comp_id))
            
            return records
    except Exception as e:
        st.error(f"Error fetching bound molecules for {pdb_code}: {e}")
        return []

async def fetch_interacting_residues_async(pdb_code, bm_id, session):
    url = f"https://www.ebi.ac.uk/pdbe/graph-api/pdb/bound_molecule_interactions/{pdb_code}/{bm_id}"
    
    try:
        async with session.get(url) as response:
            if response.status != 200:
                return []
            
            data = await response.json()
            if pdb_code not in data:
                return []
            
            records = []
            for entry in data[pdb_code]:
                if entry.get("bm_id") != bm_id:
                    continue
                interactions = entry.get("interactions", [])
                for interaction in interactions:
                    end_info = interaction.get("end", {})
                    chain_id = end_info.get("chain_id", "Unknown")
                    residue_number = end_info.get("author_residue_number", "Unknown")
                    chem_comp_id = end_info.get("chem_comp_id", "Unknown")
                    records.append((pdb_code.upper(), bm_id, chain_id, residue_number, chem_comp_id))
            
            return records
    except Exception as e:
        st.error(f"Error fetching interactions for {pdb_code} {bm_id}: {e}")
        return []

def get_uniprot_id(pdb_code, chain_id, cache):
    """Runs pdb_unpID.py to get the UniProt ID, using cache to avoid duplicate calls."""
    if (pdb_code, chain_id) in cache:
        return cache[(pdb_code, chain_id)]
    
    try:
        result = subprocess.run(["python", script_path, pdb_code, chain_id], capture_output=True, text=True)
        uniprot_id = result.stdout.strip()
        cache[(pdb_code, chain_id)] = uniprot_id  # Store result in cache
        return uniprot_id
    except Exception as e:
        st.error(f"Error fetching UniProt ID for {pdb_code} {chain_id}: {e}")
        return "Unknown"

async def get_uniprot_id_threaded(executor, pdb_code, chain_id, cache):
    if (pdb_code, chain_id) in cache:
        return cache[(pdb_code, chain_id)]
    
    uniprot_id = await asyncio.get_event_loop().run_in_executor(
        executor, 
        get_uniprot_id, 
        pdb_code, 
        chain_id,
        cache
    )
    return uniprot_id

async def process_pdb_code(pdb_code, session, uniprot_cache, executor):
    molecule_records = await fetch_bound_molecules_async(pdb_code, session)

    all_interaction_records = []
    interaction_tasks = []
    
    for pdb, bm_id, ligand in molecule_records:
        task = fetch_interacting_residues_async(pdb.lower(), bm_id, session)
        interaction_tasks.append((pdb, bm_id, ligand, task))
    
    uniprot_tasks = []
    chain_cache = set()
    
    for pdb, bm_id, ligand, task in interaction_tasks:
        interaction_records = await task
        
        for record in interaction_records:
            pdb_chain = (record[0], record[2])
            
            if pdb_chain not in uniprot_cache and pdb_chain not in chain_cache:
                chain_cache.add(pdb_chain)
                task = get_uniprot_id_threaded(executor, record[0], record[2], uniprot_cache)
                uniprot_tasks.append((pdb_chain, task))
            
            all_interaction_records.append((pdb, bm_id, ligand, record[2], None, record[4], record[3]))
    
    if uniprot_tasks:
        for i, ((pdb, chain), task) in enumerate(uniprot_tasks):
            uniprot_id = await task
            uniprot_cache[(pdb, chain)] = uniprot_id
    
    for i in range(len(all_interaction_records)):
        pdb = all_interaction_records[i][0]
        chain = all_interaction_records[i][3]
        uniprot_id = uniprot_cache.get((pdb, chain), "Unknown")
        all_interaction_records[i] = (
            all_interaction_records[i][0],
            all_interaction_records[i][1],
            all_interaction_records[i][2],
            all_interaction_records[i][3],
            uniprot_id,
            all_interaction_records[i][5],
            all_interaction_records[i][6]
        )
    return molecule_records, all_interaction_records

async def process_pdb_codes(pdb_codes, progress_bar, status_text):
    all_molecule_records = []
    all_interaction_records = []
    uniprot_cache = {}
    
    completed = 0
    total = len(pdb_codes)
    
    tcp_connector = aiohttp.TCPConnector(
        limit=50,
        ttl_dns_cache=300,
        keepalive_timeout=60
    )
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        timeout = aiohttp.ClientTimeout(total=60, connect=10)
        async with aiohttp.ClientSession(
            connector=tcp_connector, 
            timeout=timeout,
            headers={"Connection": "keep-alive"}
        ) as session:
            
            async def process_with_progress(pdb_code):
                nonlocal completed
                result = await process_pdb_code(pdb_code, session, uniprot_cache, executor)
                completed += 1
                progress = completed / total
                progress_bar.progress(progress)
                status_text.text(f"Completed processing {pdb_code} ({completed}/{total})...")
                return result
            
            tasks = [process_with_progress(pdb_code) for pdb_code in pdb_codes]
            results = await asyncio.gather(*tasks)
            
            for molecule_records, interaction_records in results:
                all_molecule_records.extend(molecule_records)
                all_interaction_records.extend(interaction_records)
    
    # Create dataframes
    df_molecules = pd.DataFrame(all_molecule_records, columns=["PDB", "bmID", "Ligand"])
    df_molecules.drop_duplicates(inplace=True)
    ligand_list = df_molecules["Ligand"].unique().tolist()
    
    df_interactions = pd.DataFrame(all_interaction_records, columns=["PDB", "bmID", "Ligand", "Chain", "UniProt ID", "Residue", "Position"])
    # Remove ligand-ligand interactions and interactions with water (HOH)
    df_interactions = df_interactions[~df_interactions["Residue"].isin(ligand_list) & (df_interactions["Residue"] != "HOH")]
    df_interactions.drop_duplicates(inplace=True)
    
    return df_interactions

def filter_and_display_data(df_interactions, selected_ligand):
    """Filter data by ligand and display results"""
    if selected_ligand and selected_ligand != "All Ligands":
        df_filtered = df_interactions[df_interactions["Ligand"].str.upper() == selected_ligand.upper()]
    else:
        df_filtered = df_interactions
    
    # Group and merge data
    df_merged = df_filtered.groupby(["Ligand", "UniProt ID", "Residue", "Position"]).agg({
        "bmID": lambda x: "/".join(sorted(set(x))),
        "Chain": lambda x: "/".join(sorted(set(x))),
        "PDB": lambda x: "/".join(sorted(set(x)))
    }).reset_index()
    
    df_merged = df_merged.sort_values(by="Position")
    
    return df_merged

# Sidebar for input
with st.sidebar:
    st.header("📝 New query")
    pdb_input = st.text_area("Enter PDB codes (comma or space-separated):", default_pdb)
    
    process_button = st.button("Process structures")

# Navigation bookmarks
st.sidebar.markdown("[Binding sites](#top_title)")
st.sidebar.markdown("[Ligands/Modifications](#ligands)")
st.sidebar.markdown("[PDB ligand bound pages](#pdb)")

# Main content area
if process_button and pdb_input:
    # Parse PDB codes
    pdb_codes = re.split(r'[,\s]+', pdb_input.strip().lower())
    pdb_codes = [code.strip() for code in pdb_codes if code.strip()]
    
    if pdb_codes:
        # st.info(f"Processing {len(pdb_codes)} PDB structure(s): {', '.join(pdb_codes)}")
        
        # Progress tracking
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        # Run async processing
        try:
            start_time = time.time()
            df_interactions = asyncio.run(process_pdb_codes(pdb_codes, progress_bar, status_text))
            end_time = time.time()
            
            st.session_state.df_interactions = df_interactions
            st.session_state.processing_complete = True
            
        except Exception as e:
            st.error(f"Error during processing: {e}")
    
# Display results if data is available
if st.session_state.processing_complete and st.session_state.df_interactions is not None:
    df_interactions = st.session_state.df_interactions
    
    # Get unique ligands
    ligand_list = [lig.upper() for lig in df_interactions["Ligand"].unique().tolist()]
           
    # Ligand selection
    selected_ligand = st.selectbox(
        "Select a ligand for analysis:",
        ["All Ligands"] + sorted(ligand_list),
        help="Choose a specific ligand or view all interactions"
    )
    
    # Filter and display data
    df_display = filter_and_display_data(df_interactions, selected_ligand)
    # Section 1: Binding sites
    # -----------------------------------------------------------------------------
    if not df_display.empty:
        st.markdown(f"#### Binding sites for {selected_ligand if selected_ligand != 'All Ligands' else 'All Ligands'}")
        
        # Display as table
        st.dataframe(
            df_display[["Ligand", "UniProt ID", "Residue", "Position", "Chain", "PDB"]],
            use_container_width=True,
            hide_index=True
        )
        
        # Download button
        csv_buffer = io.StringIO()
        df_display.to_csv(csv_buffer, index=False)
        st.download_button(
            label="Download Results as CSV",
            data=csv_buffer.getvalue(),
            file_name=f"binding_sites_{selected_ligand.lower() if selected_ligand != 'All Ligands' else 'all'}.csv",
            mime="text/csv"
        )
        
        # Show additional sections only when a specific ligand is selected
        if selected_ligand != "All Ligands":
            # Section 2: List of ligands (only the selected one)
            # -----------------------------------------------------------------------------
            st.markdown("<a name='ligands'></a>", unsafe_allow_html=True)
            st.markdown("#### Ligands / Modifications")
            
            with st.spinner("Fetching ligand information..."):
                ligand_name = get_ligand_name(selected_ligand)
            
            # Link to PDBe-KB Ligand with name
            if ligand_name:
                link = f"[{selected_ligand}](https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/{selected_ligand}): {ligand_name}"
            else:
                link = f"[{selected_ligand}](https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/{selected_ligand}): (Ligand name not available)"
            
            st.markdown(link, unsafe_allow_html=True)
            
            # Section 3: List of PDBs (only those containing the selected ligand)
            # -----------------------------------------------------------------------------
            st.markdown("<a name='pdb'></a>", unsafe_allow_html=True)
            st.markdown("#### PDB ligand bound pages")
            
            # Get unique PDBs from the filtered display data
            pdb_codes_filtered = sorted(set(
                pdb.strip() 
                for pdb_string in df_display["PDB"].tolist() 
                for pdb in pdb_string.split("/")
            ))
            
            pdb_str = " ".join(p.lower() for p in pdb_codes_filtered)
            available_links = ", ".join([f"[{pdb_code}](https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_code.lower()}/bound/{selected_ligand})" for pdb_code in pdb_codes_filtered])
            st.markdown(f"{available_links}")
     
# Instructions
if not st.session_state.processing_complete:
    st.info("""
    #### Instruction
    
    1. **Enter PDB code(s)** in the sidebar (comma or space-separated)
    2. **Process structures** to fetch binding site data from UniProt and PDBe
    3. **Select a ligand** to filter the interactions
    4. **Download results** as CSV for further analysis
    """)