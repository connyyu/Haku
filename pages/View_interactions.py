import streamlit as st
from streamlit.components.v1 import html
import re
import pandas as pd
import time
import asyncio
import aiohttp
import concurrent.futures
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type
import io
import os
import sys
import requests

from utils import get_ligand_name

st.set_page_config(page_title="Haku - View the interactions", page_icon="💮")

if 'df_interactions' not in st.session_state:
    st.session_state.df_interactions = None
if 'processing_complete' not in st.session_state:
    st.session_state.processing_complete = False

# Fix #4: removed the line that always reset guide to True
default_pdb = '7UPI 7TYG 7TXH 7TVG 7TVF 7T7A 7SD1 7SD0'

script_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.join(script_dir, "scripts")
if scripts_dir not in sys.path:
    sys.path.insert(0, scripts_dir)

try:
    import pdb_unpID
except ImportError as e:
    st.error(f"Failed to import pdb_unpID module: {e}")

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

if st.session_state.get('guide_bind'):
    st.info("""
    1. **Enter PDB code(s)** in the sidebar (comma or space-separated)
    2. **Process structures** to fetch interaction data from UniProt and PDBe
    3. **Select a ligand** to filter by ligand and display the ligand viewer.
    4. **Download results** as CSV, including a full list of bmID and ligand positions.
    """)

# Functions
# -----------------------------------------------------------------------------

async def fetch_bound_molecules_async(pdb_code, session):
    url = f"https://www.ebi.ac.uk/pdbe/graph-api/pdb/bound_molecules/{pdb_code}"
    try:
        async with session.get(url) as response:
            if response.status != 200:
                return []
            data = await response.json()
            records = []
            for entry in data.get(pdb_code, []):
                bm_id = entry.get("bm_id", "Unknown")
                for ligand in entry.get("composition", {}).get("ligands", []):
                    records.append((pdb_code.upper(), bm_id, ligand.get("chem_comp_id", "Unknown")))
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
            records = []
            for entry in data.get(pdb_code, []):
                if entry.get("bm_id") != bm_id:
                    continue
                for interaction in entry.get("interactions", []):
                    end = interaction.get("end", {})
                    records.append((
                        pdb_code.upper(), bm_id,
                        end.get("chain_id", "Unknown"),
                        end.get("author_residue_number", "Unknown"),
                        end.get("chem_comp_id", "Unknown"),
                        interaction.get("begin", {}).get("author_residue_number", "Unknown"),
                    ))
            return records
    except Exception as e:
        st.error(f"Error fetching interactions for {pdb_code} {bm_id}: {e}")
        return []


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=1, max=10),
    retry=retry_if_exception_type((ValueError,)),
)
def get_uniprot_id_with_retry(pdb_code, chain_id):
    uniprot_id = pdb_unpID.pdb_unpID(pdb_code, chain_id)
    if isinstance(uniprot_id, list):
        if not uniprot_id:
            raise ValueError("Empty UniProt ID list")
        uniprot_id = uniprot_id[0]
    if uniprot_id in ['none', 'unknown', 'error', '']:
        raise ValueError(f"Invalid UniProt ID: {uniprot_id}")
    return uniprot_id


def get_uniprot_id(pdb_code, chain_id, cache):
    if (pdb_code, chain_id) in cache:
        return cache[(pdb_code, chain_id)]
    try:
        uid = get_uniprot_id_with_retry(pdb_code, chain_id)
    except Exception:
        uid = "Unknown"
    cache[(pdb_code, chain_id)] = uid
    return uid


async def get_uniprot_id_threaded(executor, pdb_code, chain_id, cache, semaphore, file_lock):
    if (pdb_code, chain_id) in cache:
        return cache[(pdb_code, chain_id)]
    async with semaphore:
        async with file_lock:
            return await asyncio.get_event_loop().run_in_executor(
                executor, get_uniprot_id, pdb_code, chain_id, cache
            )


async def process_pdb_code(pdb_code, session, uniprot_cache, executor, uniprot_semaphore, file_lock):
    molecule_records = await fetch_bound_molecules_async(pdb_code, session)
    all_interaction_records = []
    interaction_tasks = [(pdb, bm_id, ligand, fetch_interacting_residues_async(pdb.lower(), bm_id, session))
                         for pdb, bm_id, ligand in molecule_records]
    uniprot_tasks = []
    chain_cache = set()
    for pdb, bm_id, ligand, task in interaction_tasks:
        records = await task
        for rec in records:
            pdb_chain = (rec[0], rec[2])
            if pdb_chain not in uniprot_cache and pdb_chain not in chain_cache:
                chain_cache.add(pdb_chain)
                uniprot_tasks.append((pdb_chain,
                    get_uniprot_id_threaded(executor, rec[0], rec[2], uniprot_cache, uniprot_semaphore, file_lock)))
            all_interaction_records.append((pdb, bm_id, ligand, rec[2], None, rec[4], rec[3], rec[5]))
    for (pdb, chain), task in uniprot_tasks:
        uniprot_cache[(pdb, chain)] = await task
    for i, rec in enumerate(all_interaction_records):
        uid = uniprot_cache.get((rec[0], rec[3]), "Unknown")
        all_interaction_records[i] = (rec[0], rec[1], rec[2], rec[3], uid, rec[5], rec[6], rec[7])
    return molecule_records, all_interaction_records


async def process_pdb_codes(pdb_codes, progress_bar, status_text):
    all_molecule_records, all_interaction_records, uniprot_cache = [], [], {}
    completed = 0
    connector = aiohttp.TCPConnector(limit=50, ttl_dns_cache=300, keepalive_timeout=60, limit_per_host=10)
    semaphore = asyncio.Semaphore(5)
    file_lock = asyncio.Lock()
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        timeout = aiohttp.ClientTimeout(total=10, connect=3, sock_connect=3, sock_read=5)
        async with aiohttp.ClientSession(connector=connector, timeout=timeout,
                                         headers={"Connection": "keep-alive"}) as session:
            async def run(pdb_code):
                nonlocal completed
                result = await process_pdb_code(pdb_code, session, uniprot_cache, executor, semaphore, file_lock)
                completed += 1
                progress_bar.progress(completed / len(pdb_codes))
                status_text.text(f"Completed {pdb_code} ({completed}/{len(pdb_codes)})...")
                return result
            for mol, inter in await asyncio.gather(*[run(p) for p in pdb_codes]):
                all_molecule_records.extend(mol)
                all_interaction_records.extend(inter)
    progress_bar.empty()
    status_text.empty()
    df_molecules = pd.DataFrame(all_molecule_records, columns=["PDB", "bmID", "Ligand"]).drop_duplicates()
    ligand_list = df_molecules["Ligand"].unique().tolist()
    df_interactions = pd.DataFrame(all_interaction_records,
                                   columns=["PDB", "bmID", "Ligand", "Chain", "UniProt ID", "Residue", "Position", "Ligand Position"])
    if df_interactions.empty:
        st.warning("No ligand found for the provided structure(s).")
        return df_interactions
    df_interactions = df_interactions[
        ~df_interactions["Residue"].isin(ligand_list) & (df_interactions["Residue"] != "HOH")
    ].drop_duplicates()
    return df_interactions


def calculate_ligand_numbers(df_merged):
    ligand_groups = []
    ligand_number = 1
    for idx, row in df_merged.iterrows():
        current_pairs = set(zip(str(row["Ligand Position"]).split("/"), str(row["PDB"]).split("/")))
        assigned = None
        for gn, gp in ligand_groups:
            if current_pairs & gp:
                assigned = gn
                ligand_groups[gn - 1] = (gn, gp | current_pairs)
                break
        if assigned is None:
            assigned = ligand_number
            ligand_groups.append((ligand_number, current_pairs))
            ligand_number += 1
        df_merged.loc[idx, "Index"] = assigned
    return df_merged


def filter_and_display_data(df_interactions, selected_ligand):
    df_filtered = (df_interactions[df_interactions["Ligand"].str.upper() == selected_ligand.upper()].copy()
                   if selected_ligand != "All Ligands" else df_interactions.copy())
    required = ["Ligand", "UniProt ID", "Residue", "Position", "bmID", "Chain", "PDB", "Ligand Position"]
    missing = [c for c in required if c not in df_filtered.columns]
    if missing:
        st.error(f"Missing columns: {', '.join(missing)}")
        return pd.DataFrame()
    df_filtered = df_filtered.astype({"Ligand Position": "string"})
    def agg(g):
        aligned = list(zip(g["Ligand Position"], g["Chain"], g["PDB"]))
        lp, ch, pb = zip(*aligned)
        return pd.Series({"Ligand Position": "/".join(lp), "Chain": "/".join(ch),
                          "PDB": "/".join(pb), "bmID": "/".join(sorted(set(g["bmID"])))})
    df_merged = (df_filtered.groupby(["Ligand", "UniProt ID", "Residue", "Position"],
                                     sort=False, group_keys=False)
                 [["Ligand Position", "Chain", "PDB", "bmID"]].apply(agg).reset_index()
                 .sort_values("Position"))
    return calculate_ligand_numbers(df_merged)


# Main
# -----------------------------------------------------------------------------
if process_button and pdb_input:
    pdb_codes = [c for c in re.split(r'[,\s]+', pdb_input.strip().lower()) if c]
    if pdb_codes:
        progress_bar = st.progress(0)
        status_text = st.empty()
        try:
            df_interactions = asyncio.run(process_pdb_codes(pdb_codes, progress_bar, status_text))
            st.session_state.df_interactions = df_interactions
            st.session_state.processing_complete = True
        except Exception as e:
            st.error(f"Error during processing: {e}")

if st.session_state.processing_complete and st.session_state.df_interactions is not None:
    df_interactions = st.session_state.df_interactions
    ligand_list = [l.upper() for l in df_interactions["Ligand"].unique()]

    selected_ligand = st.selectbox("Select a ligand for analysis:",
                                   ["All Ligands"] + sorted(ligand_list),
                                   help="Select a ligand to display the ligand viewer")

    df_display = filter_and_display_data(df_interactions, selected_ligand)

    if not df_display.empty:
        st.markdown(f"#### Binding sites for {selected_ligand if selected_ligand != 'All Ligands' else 'All Ligands'}")

        with st.expander("🔍 Filter results"):
            individual_pdbs = sorted({p for entry in df_display["PDB"] for p in entry.split('/')})
            combined_pdbs = sorted(df_display["PDB"].unique())
            pdb_options = ['All'] + combined_pdbs + individual_pdbs
            uniprot_options = ['All'] + sorted(df_display["UniProt ID"].unique())
            selected_pdbs = st.multiselect("Filter by structure (PDB)", pdb_options, default=['All'])
            selected_uniprots = st.multiselect("Filter by protein (UniProt ID)", uniprot_options, default=['All'])
            df_filtered = df_display.copy()
            if 'All' not in selected_pdbs:
                def pdb_match(e): return any(s == e or s in e.split('/') for s in selected_pdbs)
                df_filtered = df_filtered[df_filtered["PDB"].apply(pdb_match)]
            if 'All' not in selected_uniprots:
                df_filtered = df_filtered[df_filtered["UniProt ID"].isin(selected_uniprots)]

        display_columns = ["Ligand", "Index", "UniProt ID", "Residue", "Position", "Chain", "PDB"]
        if "Index" in df_filtered.columns:
            df_filtered["Index"] = df_filtered["Index"].astype(int)
        st.dataframe(df_filtered[display_columns], use_container_width=True, hide_index=True)

        csv_buf = io.StringIO()
        df_csv = df_filtered.copy().rename(columns={"Index": "Ligand Index", "Position": "Protein Position"})
        df_csv.to_csv(csv_buf, index=False)
        fname = f"binding_sites_{selected_ligand.lower() if selected_ligand != 'All Ligands' else 'all'}.csv"
        st.download_button("Download Results as CSV", data=csv_buf.getvalue(), file_name=fname, mime="text/csv")

        col1, col2 = st.columns([2, 1])
        with col1:
            if selected_ligand != "All Ligands":
                st.markdown("<a name='2D-viewer'></a>", unsafe_allow_html=True)
                st.markdown("#### Ligand viewer")
                first = df_interactions[df_interactions["Ligand"].str.upper() == selected_ligand.upper()].iloc[0]
                viewer_html = f"""<!DOCTYPE html><html><head>
                    <meta charset="utf-8">
                    <script src="https://d3js.org/d3.v5.min.js"></script>
                    <script src="https://code.jquery.com/jquery-3.3.1.min.js"></script>
                    <link href="https://ebi.emblstatic.net/web_guidelines/EBI-Icon-fonts/v1.3/fonts.css" rel="stylesheet"/>
                    <link href="https://www.ebi.ac.uk/pdbe/pdb-component-library/css/pdb-ligand-env-svg.css" rel="stylesheet"/>
                    <script type="module" src="https://www.ebi.ac.uk/pdbe/pdb-component-library/js/pdb-ligand-env-component-1.0.0-min.js"></script>
                    <style>pdb-ligand-env{{display:block;width:400px;height:400px;border:1px solid #ccc9c0;border-radius:5px;}}</style>
                    </head><body>
                    <pdb-ligand-env pdb-id="{first['PDB'].lower()}" pdb-chain-id="{first['Chain']}" pdb-res-id="{first['Ligand Position']}" environment="production" zoom-on></pdb-ligand-env>
                    <p style="margin-top:12px;font-size:14px;color:#666;">
                    Ligand {selected_ligand} in {first['PDB'].upper()} (Chain: {first['Chain']}, Position: {first['Ligand Position']})</p>
                    </body></html>"""
                html(viewer_html, height=450)

                with col2:
                    st.markdown("")
                    st.markdown("##### Ligand name")
                    with st.spinner("Fetching ligand information..."):
                        ligand_name = get_ligand_name(selected_ligand)
                    desc = ligand_name if ligand_name else "(Ligand name not available)"
                    st.markdown(f"[{selected_ligand}](https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/{selected_ligand}): {desc}", unsafe_allow_html=True)

                    copies = df_filtered["Index"].max() if "Index" in df_filtered.columns else 1
                    if copies > 1:
                        st.markdown("##### Number of copies")
                        st.markdown(f"**{int(copies)}** (indexed by position)")

                    st.markdown("##### View other structures")
                    pdb_codes_filtered = sorted({p for e in df_display["PDB"] for p in e.split("/")})
                    st.markdown("<br>".join(
                        f"[{p}](https://www.ebi.ac.uk/pdbe/entry/pdb/{p.lower()}/bound/{selected_ligand})"
                        for p in pdb_codes_filtered
                    ), unsafe_allow_html=True)

    st.markdown("""
        ---
        The ligand viewer is generated using **PDBe** [pdb-ligand-env](https://gitlab.ebi.ac.uk/pdbe/web-components/ligand-env).  
        The interactions are calculated using [Arpeggio](https://www.sciencedirect.com/science/article/pii/S0022283616305332) from Blundell et al (2017) and downloaded from **PDBe** API.
    """)
