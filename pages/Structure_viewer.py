"""
Structure_viewer.py — Navigate protein structures in a 3D viewer.

Changes from v1:
  - Removed ~300 lines of duplicated CIF parsing (now in utils.py)
  - Fixed race condition: uses tempfile.mkdtemp() instead of shared chain_info/
  - Session state keys prefixed sv_ to avoid cross-page clashes
  - Uses @st.cache_data via utils for API calls
"""

import re
import os
import sys
import shutil
import tempfile
import requests
import pandas as pd
import streamlit as st
import py3Dmol
import streamlit.components.v1 as components
from collections import defaultdict

from utils import (
    parse_resolution,
    parse_hetatms,
    parse_uniprot_mapping,
    generate_pymol_commands,
    alternate_colors,
)

# ---------------------------------------------------------------------------
# Page config & session state
# ---------------------------------------------------------------------------
st.set_page_config(page_title="Haku - Structure viewer", page_icon="💮")

if "sv_all_data" not in st.session_state:
    st.session_state.sv_all_data = []
if "sv_pdb_codes" not in st.session_state:
    st.session_state.sv_pdb_codes = []

# Add scripts dir to path so pdb_chainID can be imported
script_dir = os.path.dirname(os.path.abspath(__file__))
scripts_dir = os.path.join(script_dir, "scripts")
if scripts_dir not in sys.path:
    sys.path.insert(0, scripts_dir)

try:
    import pdb_chainID
except ImportError as e:
    st.error(f"Failed to import pdb_chainID module: {e}")

# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
default_pdb = "9eii 9eij 9eih"

with st.sidebar:
    st.header("📝 New query")
    pdb_input = st.text_area("Enter PDB codes (comma or space-separated):", default_pdb).strip()
    fetch_button = st.button("Analyse structure")
    if fetch_button:
        st.session_state.sv_guide = False
    st.sidebar.markdown("[Proteins](#simplified-data)")
    st.sidebar.markdown("[Structures](#structure)")
    st.sidebar.markdown("[PyMOL commands](#pymol-commands)")
    sv_guide = st.checkbox("Instructions", value=st.session_state.get("sv_guide", True))
    st.session_state.sv_guide = sv_guide

st.markdown("<a name='top_title'></a>", unsafe_allow_html=True)
st.markdown("#### Navigate protein structures in a 3D viewer.")

if st.session_state.get("sv_guide"):
    st.info("""
    Input:&emsp;**one or multiple PDB codes**  
      
    Output:
    * all protein components mapped in UniProt
    * resolution and method of each structure
    * **selected protein shown in a 3D viewer**
    * PyMOL command to colour and label each protein in the structures
    """)

# ---------------------------------------------------------------------------
# CIF processing — uses utils parsers, temp dir per run (no race condition)
# ---------------------------------------------------------------------------

def download_cif(pdb_code: str, output_dir: str) -> str | None:
    url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_code.lower()}_updated.cif"
    try:
        r = requests.get(url, timeout=15)
        if r.status_code == 200:
            path = os.path.join(output_dir, f"{pdb_code.lower()}.cif")
            with open(path, "wb") as f:
                f.write(r.content)
            return path
    except requests.exceptions.RequestException:
        pass
    return None


def check_chimera(table_data: list) -> list:
    chain_map = defaultdict(set)
    pdb_code = ""
    for row in table_data:
        pdb_code = row["PDB"].upper()
        chain_map[row["Chain"]].add(row["UniProt AC"])
    return [
        f"{pdb_code} chain {ch} maps to multiple UniProt ACs: {', '.join(acs)}."
        for ch, acs in chain_map.items()
        if len(acs) > 1
    ]


def process_cif(pdb_code: str, tmp_dir: str) -> list:
    cif_path = download_cif(pdb_code, tmp_dir)
    if not cif_path:
        st.warning(f"Could not download CIF for {pdb_code.upper()}.")
        return []

    resolution, method = parse_resolution(cif_path)
    hetatms = parse_hetatms(cif_path)
    uniprot_mapping, uniprot_ids, chain_mapping = parse_uniprot_mapping(cif_path)

    hetatm_str = ", ".join(hetatms)
    table_data = []
    seen = set()

    for entity_id, ac in uniprot_mapping.items():
        uid = next((u for e, u in uniprot_ids if e == entity_id), None)
        for chain_id in chain_mapping.get(entity_id, ["A"]):
            key = (pdb_code.lower(), ac, chain_id, entity_id)
            if key in seen:
                continue
            seen.add(key)
            table_data.append({
                "PDB": pdb_code.lower(),
                "Entity": entity_id,
                "Chain": chain_id,
                "UniProt AC": ac,
                "UniProt ID": uid,
                "HETATM": hetatm_str,
                "Resolution": resolution,
                "Method": method,
            })

    for w in check_chimera(table_data):
        st.warning(f"⚠️ Chimera: {w}")

    return table_data


def simplify(all_data: list) -> list:
    grouped = defaultdict(list)
    for row in all_data:
        key = (row["PDB"], row["Entity"], row["UniProt AC"],
               row["UniProt ID"], row["HETATM"], row["Resolution"], row["Method"])
        grouped[key].append(row["Chain"])
    result = []
    for (pdb, ent, ac, uid, hetatm, res, method), chains in grouped.items():
        result.append({
            "PDB": pdb, "Entity": ent,
            "Chain": "/".join(sorted(chains)),
            "UniProt AC": ac, "UniProt ID": uid,
            "HETATM": hetatm, "Resolution": res, "Method": method,
        })
    return sorted(result, key=lambda x: (x["PDB"], x["Entity"]))


# ---------------------------------------------------------------------------
# 3D viewer helpers
# ---------------------------------------------------------------------------

def showmol(view, height=450, width=600):
    components.html(view._make_html(), height=height, width=width)


@st.cache_data(show_spinner=False)
def fetch_pdb_structure(pdb_code: str) -> str | None:
    url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_code.lower()}_updated.cif"
    try:
        r = requests.get(url, timeout=15)
        if r.status_code == 200:
            return r.text
    except requests.exceptions.RequestException:
        pass
    st.error(f"Failed to fetch structure for {pdb_code.upper()}.")
    return None


@st.cache_data(show_spinner=False)
def fetch_uniprot_sequence(uniprot_ac: str) -> str | None:
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}.fasta"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            return r.text.split("\n", 1)[1].replace("\n", "").replace("\r", "")
    except requests.exceptions.RequestException:
        pass
    st.error(f"Failed to fetch sequence for {uniprot_ac}.")
    return None


def viewpdb(structure: str, pdb_code: str, uniprot_ac: str, show_all: bool = True):
    view = py3Dmol.view(height=450, width=600)
    view.addModel(structure, "mmcif")
    view.setBackgroundColor("#eeeeee")
    view.spin(False)

    chain_ids_str = ""
    try:
        chain_info, _ = pdb_chainID.get_chain_ids(pdb_code, uniprot_ac)
        chain_ids_list = [e["Chain ID"] for e in chain_info if "Chain ID" in e]
        chain_ids_str = ", ".join(chain_ids_list)
        st.session_state.sv_chain_ids = chain_ids_str
    except Exception as e:
        st.error(f"Failed to get chain IDs: {e}")
        chain_ids_list = ["A"]

    if show_all:
        view.setStyle({"model": 0}, {"cartoon": {"color": "powderblue"}})
    for chain_id in chain_ids_list:
        view.setStyle({"chain": chain_id}, {"cartoon": {"color": "#830592"}})

    view.zoomTo()
    showmol(view)

    if chain_ids_str:
        st.caption(f"PDB:{pdb_code.upper()} ({chain_ids_str}) — {uniprot_ac.upper()} highlighted in purple.")
    else:
        st.caption(f"PDB:{pdb_code.upper()} — {uniprot_ac.upper()} highlighted. Chain mapping incomplete.")


# ---------------------------------------------------------------------------
# Section 1: Fetch and parse CIF files
# ---------------------------------------------------------------------------
progress_bar = st.progress(0)
status_text = st.empty()

if fetch_button:
    pdb_codes = [c for c in re.split(r"[,\s]+", pdb_input.strip().upper()) if c]
    st.session_state.sv_pdb_codes = pdb_codes
    all_data = []

    tmp_dir = tempfile.mkdtemp()
    try:
        for i, pdb_code in enumerate(pdb_codes, 1):
            status_text.text(f"Analysing {pdb_code} ({i}/{len(pdb_codes)})...")
            rows = process_cif(pdb_code, tmp_dir)
            all_data.extend(rows)
            progress_bar.progress(i / len(pdb_codes))
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    status_text.empty()
    st.session_state.sv_all_data = [dict(t) for t in {tuple(d.items()) for d in all_data}]

progress_bar.empty()
status_text.empty()

all_data = st.session_state.sv_all_data
pdb_codes = st.session_state.sv_pdb_codes

# ---------------------------------------------------------------------------
# Section 2: Proteins table
# ---------------------------------------------------------------------------
if all_data:
    simplified = simplify(all_data)
    df = pd.DataFrame(simplified)

    st.markdown("<a name='simplified-data'></a>", unsafe_allow_html=True)
    st.markdown("#### Proteins")
    st.dataframe(alternate_colors(df), hide_index=True, use_container_width=True)

# ---------------------------------------------------------------------------
# Section 3: 3D Structure viewer
# ---------------------------------------------------------------------------
st.markdown("<a name='structure'></a>", unsafe_allow_html=True)

if all_data and pdb_codes:
    st.markdown("#### Structures")
    col1, col2 = st.columns([1, 3])

    with col1:
        pdb_code = st.selectbox("Select a PDB code:", pdb_codes, index=0)
        filtered_rows = [r for r in all_data if r["PDB"].upper() == pdb_code.upper()]
        uid_to_ac = {r["UniProt ID"]: r["UniProt AC"] for r in filtered_rows}
        unique_uids = sorted(uid_to_ac.keys())

        if unique_uids:
            selected_uid = st.selectbox("UniProt ID: (purple)", unique_uids, index=0)
            uniprot_ac = uid_to_ac[selected_uid]
            show_all = st.checkbox("Show all chains", value=True)
            st.markdown("")
            st.caption("Only UniProt entries can be selected. Ligands and solvents are not shown.")

            pdb_structure = fetch_pdb_structure(pdb_code)

            with col2:
                if pdb_structure:
                    viewpdb(pdb_structure, pdb_code, uniprot_ac, show_all)

# ---------------------------------------------------------------------------
# Section 4: PyMOL commands
# ---------------------------------------------------------------------------
st.markdown("<a name='pymol-commands'></a>", unsafe_allow_html=True)

if all_data:
    st.markdown("#### PyMOL Commands")
    st.code("fetch " + "; fetch ".join(p.lower() for p in pdb_codes))
    st.code(generate_pymol_commands(all_data))

    if st.button("Other PyMOL Commands"):
        st.switch_page("pages/Pymol_commands.py")
