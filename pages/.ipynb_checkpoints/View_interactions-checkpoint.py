import streamlit as st
from streamlit.components.v1 import html
import re
import pandas as pd
import subprocess
import time
import asyncio
import aiohttp
import concurrent.futures
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type
import io
import os
import glob
import sys
from pathlib import Path
import requests

# Sidebar, title, parameters
# -----------------------------------------------------------------------------
st.set_page_config(page_title="Haku - View the interactions", page_icon="💮")

# Initialize session state
if 'df_interactions' not in st.session_state:
    st.session_state.df_interactions = None
if 'processing_complete' not in st.session_state:
    st.session_state.processing_complete = False

st.session_state.guide = True
default_pdb = '7UPI 7TYG 7TXH 7TVG 7TVF 7T7A 7SD1 7SD0'

with st.sidebar:
    st.header("📝 New query")
    pdb_input = st.text_area("Enter PDB codes (comma or space-separated):", default_pdb)
    process_button = st.button("Process structures")
    if process_button:
        st.session_state.guide_bind = False
    st.sidebar.markdown("[Binding sites](#top_title)")
    st.sidebar.markdown("[Ligand viewer](#2D-viewer)")
    guide_bind = st.checkbox("Instructions", value=st.session_state.get("guide_bind", True))
    st.session_state.guide_bind = guide_bind

st.markdown("<a name='top_title'></a>", unsafe_allow_html=True)
st.markdown("#### View all the ligand-protein interactions in the structure(s).")

# Get the directory where this script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
# Add the scripts subdirectory to Python path
scripts_dir = os.path.join(script_dir, "scripts")
if scripts_dir not in sys.path:
    sys.path.insert(0, scripts_dir)

# Import the pdb_unpID function from the scripts folder
try:
    import pdb_unpID
except ImportError as e:
    st.error(f"Failed to import pdb_unpID module: {e}")

# Instructions
# -----------------------------------------------------------------------------
if st.session_state.get('guide_bind'):
    st.info("""  
    1. **Enter PDB code(s)** in the sidebar (comma or space-separated)  
    2. **Process structures** to fetch interaction data from UniProt and PDBe  
    3. **Select a ligand** to filter by ligand and display the ligand viewer. The results can be further filtered by structure or by protein. 
    4. **Download results** as CSV, including a full list of bmID and ligand positions.
    """)

# Define functions
# -----------------------------------------------------------------------------
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
                    ligand_position = interaction.get("begin", {}).get("author_residue_number", "Unknown")
                    records.append((pdb_code.upper(), bm_id, chain_id, residue_number, chem_comp_id, ligand_position))
            
            return records
    except Exception as e:
        st.error(f"Error fetching interactions for {pdb_code} {bm_id}: {e}")
        return []

@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=1, max=10),
    retry=retry_if_exception_type((subprocess.CalledProcessError, subprocess.TimeoutExpired, ValueError))
)

def get_uniprot_id_with_retry(pdb_code, chain_id):
    """Call the pdb_unpID function with retry logic"""
    try:
        # Call the function from the imported module
        uniprot_id = pdb_unpID.pdb_unpID(pdb_code, chain_id)

        # Convert list to string
        if isinstance(uniprot_id, list):
            if not uniprot_id:
                raise ValueError("Empty UniProt ID list")
            uniprot_id = uniprot_id[0]
        
        # Validate the result
        if uniprot_id in ['none', 'unknown', 'error', '']:
            raise ValueError(f"Invalid UniProt ID result: {uniprot_id}")
        
        return uniprot_id

    except Exception as e:
        print(f"Error fetching UniProt ID for {pdb_code} {chain_id}: {e}")
        raise

def get_uniprot_id(pdb_code, chain_id, cache):
    """Get UniProt ID using the pdb_unpID function, with caching to avoid duplicate calls."""
    if (pdb_code, chain_id) in cache:
        return cache[(pdb_code, chain_id)]
    
    try:
        uniprot_id = get_uniprot_id_with_retry(pdb_code, chain_id)
        cache[(pdb_code, chain_id)] = uniprot_id  # Store result in cache
        return uniprot_id
    except Exception as e:
        print(f"Failed to get UniProt ID for {pdb_code} {chain_id} after retries: {e}")
        cache[(pdb_code, chain_id)] = "Unknown"  # Cache the failure to avoid retrying
        return "Unknown"

async def get_uniprot_id_threaded(executor, pdb_code, chain_id, cache, semaphore, file_lock):
    """Async wrapper for getting UniProt IDs using thread pool"""
    if (pdb_code, chain_id) in cache:
        return cache[(pdb_code, chain_id)]
    
    # Use semaphore to limit concurrent UniProt requests
    async with semaphore:
        async with file_lock:
            # Run the blocking function in a thread pool
            uniprot_id = await asyncio.get_event_loop().run_in_executor(
                executor, 
                get_uniprot_id, 
                pdb_code, 
                chain_id,
                cache
            )
            return uniprot_id

async def process_pdb_code(pdb_code, session, uniprot_cache, executor, uniprot_semaphore, file_lock):
    molecule_records = await fetch_bound_molecules_async(pdb_code, session)

    # Process all interactions for each molecule in parallel
    all_interaction_records = []
    interaction_tasks = []
    
    for pdb, bm_id, ligand in molecule_records:
        task = fetch_interacting_residues_async(pdb.lower(), bm_id, session)
        interaction_tasks.append((pdb, bm_id, ligand, task))
    
    # Process the interaction results and prepare UniProt lookups
    uniprot_tasks = []
    chain_cache = set()  # Track chains already queued for lookup
    
    for pdb, bm_id, ligand, task in interaction_tasks:
        interaction_records = await task
        
        for record in interaction_records:
            pdb_chain = (record[0], record[2])
            
            # Only queue UniProt lookups for unseen chains
            if pdb_chain not in uniprot_cache and pdb_chain not in chain_cache:
                chain_cache.add(pdb_chain)
                task = get_uniprot_id_threaded(executor, record[0], record[2], uniprot_cache, uniprot_semaphore, file_lock)
                uniprot_tasks.append((pdb_chain, task))
            
            # Add the interactions to the results in the next step once we have all UniProt IDs
            all_interaction_records.append((pdb, bm_id, ligand, record[2], None, record[4], record[3], record[5]))
    
    if uniprot_tasks:
        for i, ((pdb, chain), task) in enumerate(uniprot_tasks):
            uniprot_id = await task
            uniprot_cache[(pdb, chain)] = uniprot_id
    
    for i in range(len(all_interaction_records)):
        pdb = all_interaction_records[i][0]
        chain = all_interaction_records[i][3]
        uniprot_id = uniprot_cache.get((pdb, chain), "Unknown")
        all_interaction_records[i] = (
            all_interaction_records[i][0],  # PDB
            all_interaction_records[i][1],  # bmID
            all_interaction_records[i][2],  # Ligand
            all_interaction_records[i][3],  # Chain
            uniprot_id,
            all_interaction_records[i][5],  # Residue
            all_interaction_records[i][6],  # Position
            all_interaction_records[i][7]   # Ligand position
        )

    return molecule_records, all_interaction_records

async def process_pdb_codes(pdb_codes, progress_bar, status_text):
    start_time = time.time()
    all_molecule_records = []
    all_interaction_records = []
    uniprot_cache = {}
    
    completed = 0
    total = len(pdb_codes)
    
    tcp_connector = aiohttp.TCPConnector(
        limit=50,
        ttl_dns_cache=300,  # Cache DNS results for 5 minutes
        keepalive_timeout=60,  # Keep connections alive for 60 seconds
        limit_per_host=10
    )
    
    # Create semaphore to limit concurrent UniProt requests
    uniprot_semaphore = asyncio.Semaphore(5)
    file_lock = asyncio.Lock()  # File access lock to prevent concurrent file operations
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        # Create a shared session for all HTTP requests with optimized settings
        timeout = aiohttp.ClientTimeout(total=10, connect=3, sock_connect=3, sock_read=5)
        async with aiohttp.ClientSession(
            connector=tcp_connector, 
            timeout=timeout,
            headers={"Connection": "keep-alive"}
        ) as session:
            
            async def process_with_progress(pdb_code):
                nonlocal completed
                result = await process_pdb_code(pdb_code, session, uniprot_cache, executor, uniprot_semaphore, file_lock)
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
 
            # Clear progress bar
            progress_bar.empty()
            status_text.empty()   
            
    # Create dataframes
    df_molecules = pd.DataFrame(all_molecule_records, columns=["PDB", "bmID", "Ligand"])
    df_molecules.drop_duplicates(inplace=True)
    ligand_list = df_molecules["Ligand"].unique().tolist()
    
    df_interactions = pd.DataFrame(all_interaction_records, columns=["PDB", "bmID", "Ligand", "Chain", "UniProt ID", "Residue", "Position", "Ligand Position"])
    if df_interactions.empty:
        st.warning("No ligand found for the provided structure(s).")

    # Remove ligand-ligand interactions and interactions with water (HOH)
    df_interactions = df_interactions[~df_interactions["Residue"].isin(ligand_list) & (df_interactions["Residue"] != "HOH")]
    df_interactions.drop_duplicates(inplace=True)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Script executed in {elapsed_time:.2f} seconds.")   
    
    return df_interactions

def calculate_ligand_numbers(df_merged):
    """Calculate smart Indexs based on position-PDB overlap"""
    # Build ligand groups based on position-PDB overlap
    ligand_groups = []
    ligand_number = 1
    
    for idx, row in df_merged.iterrows():
        ligand_pos = row["Ligand Position"]
        pdb_list = row["PDB"]
        
        # Create position-PDB pairs for current row
        current_positions = str(ligand_pos).split("/")
        current_pdbs = pdb_list.split("/")
        current_pairs = set(zip(current_positions, current_pdbs))
        
        # Check if this ligand overlaps with any existing group
        assigned_number = None
        for group_num, group_pairs in ligand_groups:
            if current_pairs.intersection(group_pairs):
                assigned_number = group_num
                # Update the group with new pairs
                ligand_groups[group_num - 1] = (group_num, group_pairs.union(current_pairs))
                break
        
        # If no overlap found, create new group
        if assigned_number is None:
            assigned_number = ligand_number
            ligand_groups.append((ligand_number, current_pairs))
            ligand_number += 1
        
        # Store the assignment
        df_merged.loc[idx, "Index"] = assigned_number
    
    return df_merged

def filter_and_display_data(df_interactions, selected_ligand):
    """Filter data by ligand and display results with Indexs"""
    if selected_ligand and selected_ligand != "All Ligands":
        df_filtered = df_interactions[df_interactions["Ligand"].str.upper() == selected_ligand.upper()].copy()
    else:
        df_filtered = df_interactions.copy()

    required_columns = ["Ligand", "UniProt ID", "Residue", "Position", "bmID", "Chain", "PDB", "Ligand Position"]
    missing = [col for col in required_columns if col not in df_filtered.columns]
    
    if missing:
        st.error(f"Missing required column(s) in filtered data: {', '.join(missing)}")
        return pd.DataFrame()  # Return empty DataFrame if columns are missing

    # Group by residue-related fields and align Ligand Position with PDB and Chain
    group_cols = ["Ligand", "UniProt ID", "Residue", "Position"]
    # Convert Ligand Position column to string dtype to avoid future warnings
    df_filtered = df_filtered.astype({"Ligand Position": "string"})
    
    def aggregate_aligned(group):
        aligned = list(zip(group["Ligand Position"], group["Chain"], group["PDB"]))
        ligand_positions, chains, pdbs = zip(*aligned)
        return pd.Series({
            "Ligand Position": "/".join(ligand_positions),
            "Chain": "/".join(chains),
            "PDB": "/".join(pdbs),
            "bmID": "/".join(sorted(set(group["bmID"])))
        })
    
    df_merged = df_filtered.groupby(group_cols, sort=False, group_keys=False)[["Ligand Position", "Chain", "PDB", "bmID"]].apply(aggregate_aligned).reset_index()
    df_merged = df_merged.sort_values(by="Position")
    
    # Calculate Indexs
    df_merged = calculate_ligand_numbers(df_merged)
    
    return df_merged

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
           
    # Section 1: Select a ligand
    # -----------------------------------------------------------------------------
    selected_ligand = st.selectbox(
        "Select a ligand for analysis:",
        ["All Ligands"] + sorted(ligand_list),
        help="Select a ligand to display the ligand viewer"
    )
    
    # Section 2: Filter and display data
    # -----------------------------------------------------------------------------
    df_display = filter_and_display_data(df_interactions, selected_ligand)
    
    if not df_display.empty:
        st.markdown(f"#### Binding sites for {selected_ligand if selected_ligand != 'All Ligands' else 'All Ligands'}")
    
        with st.expander("🔍 Filter results"):
            # Extract all individual PDB codes from splitting
            individual_pdbs = set()
            for entry in df_display["PDB"]:
                individual_pdbs.update(entry.split('/'))
            individual_pdbs = sorted(individual_pdbs)
    
            # Get all original combined entries
            combined_pdbs = sorted(df_display["PDB"].unique())
    
            # Combine for filter options: 'All' + combined + individual
            pdb_options = ['All'] + combined_pdbs + individual_pdbs
    
            uniprot_options = ['All'] + sorted(df_display["UniProt ID"].unique())
    
            selected_pdbs = st.multiselect("Filter by structure (PDB)", options=pdb_options, default=['All'])
            selected_uniprots = st.multiselect("Filter by protein (UniProt ID)", options=uniprot_options, default=['All'])
    
            df_filtered = df_display.copy()
    
            if 'All' not in selected_pdbs:
                def pdb_match(pdb_entry):
                    # If exact combined selected, match exact
                    if any(sel == pdb_entry for sel in selected_pdbs):
                        return True
                    # If any individual code selected matches a part in pdb_entry
                    parts = pdb_entry.split('/')
                    return any(sel in parts for sel in selected_pdbs)
    
                df_filtered = df_filtered[df_filtered["PDB"].apply(pdb_match)]
    
            if 'All' not in selected_uniprots:
                df_filtered = df_filtered[df_filtered["UniProt ID"].isin(selected_uniprots)]
    
        # Columns to display (excluding 'Ligand Position')
        display_columns = ["Ligand", "Index", "UniProt ID", "Residue", "Position", "Chain", "PDB"]
        
        # Display dataframe
        if "Index" in df_filtered.columns:
            df_filtered["Index"] = df_filtered["Index"].astype(int)
            df_display = df_filtered[display_columns]
        else:
            df_display = df_filtered[["Ligand", "UniProt ID", "Residue", "Position", "Chain", "PDB"]]
        
        st.dataframe(
            df_display,
            use_container_width=True,
            hide_index=True
        )
        
        # Prepare CSV export
        csv_buffer = io.StringIO()
        df_csv = df_filtered.copy()
        
        # Rename 'Index' and 'Position' in CSV
        if "Index" in df_csv.columns:
            df_csv = df_csv.rename(columns={"Index": "Ligand Index"})
        if "Position" in df_csv.columns:
            df_csv = df_csv.rename(columns={"Position": "Protein Position"})
        
        df_csv.to_csv(csv_buffer, index=False)
        
        # Download button
        st.download_button(
            label="Download Results as CSV",
            data=csv_buffer.getvalue(),
            file_name=f"binding_sites_{selected_ligand.lower() if selected_ligand != 'All Ligands' else 'all'}.csv",
            mime="text/csv"
        )
        
        # Showing data in 2 columns
        col1, col2 = st.columns([2, 1])
        with col1:
            # Show additional sections only when a specific ligand is selected
            if selected_ligand != "All Ligands":
                # Section 3: 2D Ligand Environment Viewer
                # -----------------------------------------------------------------------------
                st.markdown("<a name='2D-viewer'></a>", unsafe_allow_html=True)
                st.markdown("#### Ligand viewer")
                
                # Get the first PDB and chain from the filtered results for the viewer
                first_match = df_interactions[df_interactions["Ligand"].str.upper() == selected_ligand.upper()].iloc[0]
                first_pdb = first_match["PDB"].lower()
                first_chain = first_match["Chain"]
                first_position = first_match["Ligand Position"]
                
                # Create the 2D viewer HTML with corrected parameters
                viewer_html = f"""
                <!DOCTYPE html>
                <html>
                <head>
                    <meta charset="utf-8">
                    <meta name="viewport" content="width=device-width, initial-scale=1.0">
                    <script src="https://d3js.org/d3.v5.min.js"></script>
                    <script src="https://code.jquery.com/jquery-3.3.1.min.js"></script>
                    <link href="https://ebi.emblstatic.net/web_guidelines/EBI-Icon-fonts/v1.3/fonts.css" type="text/css" rel="stylesheet" media="screen, projection" />
                    <link href="https://www.ebi.ac.uk/pdbe/pdb-component-library/css/pdb-ligand-env-svg.css" type="text/css" rel="stylesheet" media="screen, projection" />
                    <link href="https://cdn.jsdelivr.net/npm/pdbe-molstar@3.3.2/build/pdbe-molstar-light.css" type="text/css" rel="stylesheet" media="screen, projection" />
                    <script src="https://www.ebi.ac.uk/pdbe/pdb-component-library/js/pdbe-molstar-plugin-1.1.0.js"></script>
                    <script type="module" src="https://www.ebi.ac.uk/pdbe/pdb-component-library/js/pdb-ligand-env-component-1.0.0-min.js"></script>
                    
                    <style>
                        pdb-ligand-env {{
                            display: block;
                            width: 400px;
                            height: 400px;
                            border: 1px solid #ccc9c0;
                            border-radius: 5px;
                        }}
                        .loading {{
                            display: flex;
                            justify-content: center;
                            align-items: center;
                            height: 400px;
                            font-family: Arial, sans-serif;
                            color: #666;
                            background-color: #f8f9fa;
                        }}
                    </style>
                </head>
                <body>
                    <div style="width: 400px; height: 400px; position: top;">
                        <div class="loading" id="loading">Loading 2D viewer...</div>
                        <pdb-ligand-env
                            pdb-id="{first_pdb}" 
                            pdb-chain-id="{first_chain}"
                            pdb-res-id="{first_position}"
                            environment="production"
                            zoom-on>
                        </pdb-ligand-env>
                    </div>
                    
                    <script>
                        // Wait for component to load
                        setTimeout(function() {{
                            const loading = document.getElementById('loading');
                            const component = document.querySelector('pdb-ligand-env');
                            
                            if (component) {{
                                loading.style.display = 'none';
                                console.log('PDB Ligand Environment component initialized');
                            }}
                        }}, 3000);
                        
                        // Error handling
                        window.addEventListener('error', function(e) {{
                            console.error('Error loading 2D viewer:', e);
                            document.getElementById('loading').innerHTML = 'Error loading 2D viewer. Please check the PDB ID and parameters.';
                        }});
                    </script>
                    
                    <p style="margin-top: 12px; font-size: 14px; color: #666;">
                        Ligand {selected_ligand} in {first_pdb.upper()} (Chain: {first_chain}, Ligand position: {first_position})
                    </p>
                </body>
                </html>
                """
                
                # Display the 2D viewer
                html(viewer_html, height=450)
                
                with col2:                
                    # Section 4: List of ligands (only the selected one)
                    # -----------------------------------------------------------------------------
                    st.markdown("")
                    st.markdown("")
                    st.markdown("##### Ligand name")
                    
                    with st.spinner("Fetching ligand information..."):
                        ligand_name = get_ligand_name(selected_ligand)
                    
                    # Link to PDBe-KB Ligand with name
                    if ligand_name:
                        link = f"[{selected_ligand}](https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/{selected_ligand}): {ligand_name}"
                    else:
                        link = f"[{selected_ligand}](https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/{selected_ligand}): (Ligand name not available)"
                    
                    st.markdown(link, unsafe_allow_html=True)

                    copies = df_filtered["Index"].max()
                    if copies > 1:
                        st.markdown("##### Number of copies")
                        st.markdown(f"**{copies}** (indexed by position)")
    
                    # Section 5: PDB ligand bound pages (only those containing the selected ligand)
                    # -----------------------------------------------------------------------------
                    st.markdown("##### View other structures")
                    
                    # Get unique PDBs from the filtered display data
                    pdb_codes_filtered = sorted(set(
                        pdb.strip() 
                        for pdb_string in df_display["PDB"].tolist() 
                        for pdb in pdb_string.split("/")
                    ))
                    available_links = "<br>".join([
                        f"[{pdb_code}](https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_code.lower()}/bound/{selected_ligand})"
                        for pdb_code in pdb_codes_filtered
                    ])
                    st.markdown(f"{available_links}", unsafe_allow_html=True)

    # Section 6: Acknowledgement
    # -----------------------------------------------------------------------------
    
    st.markdown(
        """
        ---
        The ligand viewer is generated using the **PDBe** [pdb-ligand-env](https://gitlab.ebi.ac.uk/pdbe/web-components/ligand-env).  
        The interactions are calculated from the structures using [Arpeggio](https://www.sciencedirect.com/science/article/pii/S0022283616305332) from Blundell et al (2017) and downloaded from **PDBe** API.
        """
    )
