import re
import tempfile
import requests
import pandas as pd
import streamlit as st
from collections import defaultdict

from utils import (
    download_cif_file, parse_resolution, parse_hetatms, parse_uniprot_mapping,
    get_gene_name_from_uniprot, get_protein_name_from_uniprot,
    get_ligand_name, get_structure_title,
    generate_pymol_commands, alternate_colors,
)

st.set_page_config(page_title="Haku - List the proteins", page_icon="💮")

if 'list_all_data' not in st.session_state:
    st.session_state.list_all_data = []
if 'list_pdb_codes' not in st.session_state:
    st.session_state.list_pdb_codes = []

default_pdb = '9eii 9eij 9eih'

with st.sidebar:
    st.header("📝 New query")
    pdb_input = st.text_area("Enter PDB codes (comma or space-separated):", default_pdb).strip()
    fetch_button = st.button("Analyse structure")
    if fetch_button:
        st.session_state.guide_protein = False
    st.sidebar.markdown("[Proteins and ligands](#simplified-data)")
    st.sidebar.markdown("[UniProt list](#uniprot-ac)")
    st.sidebar.markdown("[PDB list](#pdb)")
    st.sidebar.markdown("[Ligands/Modifications](#ligands)")
    st.sidebar.markdown("[PyMOL commands](#pymol-commands)")
    guide_protein = st.checkbox("Instructions", value=st.session_state.get("guide_protein", True))
    st.session_state.guide_protein = guide_protein

st.markdown("<a name='top_title'></a>", unsafe_allow_html=True)
st.markdown("#### List all the UniProt protein entries in the structure(s).")

if st.session_state.get('guide_protein'):
    st.info("""
    Input:&emsp;**one or multiple PDB codes**\n\n
    Output:\n
    * all protein components mapped in UniProt\n
    * resolution and method of each structure\n
    * name of each protein and links to UniProt\n
    * all the bound ligands / modifications\n
    * PyMOL command to colour and label each protein
    """)

def check_chimera(table_data):
    chain_uniprot_map = defaultdict(set)
    for row in table_data:
        chain_uniprot_map[row['Chain']].add(row['UniProt AC'])
    pdb = table_data[0]['PDB'].upper() if table_data else '?'
    return [
        f"{pdb} chain {ch} is associated with multiple UniProt ACs: {', '.join(acs)}."
        for ch, acs in chain_uniprot_map.items() if len(acs) > 1
    ]

def process_cif_file(pdb_code, tmp_dir):
    # Fix #3: tmp_dir is per-run so concurrent users never collide
    cif_path = download_cif_file(pdb_code, tmp_dir)
    if not cif_path:
        st.warning(f"Could not download CIF for {pdb_code.upper()}.")
        return []
    resolution, method = parse_resolution(cif_path)
    hetatms = parse_hetatms(cif_path)
    uniprot_mapping, uniprot_ids, chain_mapping = parse_uniprot_mapping(cif_path)
    hetatm_names = ', '.join(hetatms)
    table_data = []
    seen = set()
    for entity_id, acc in uniprot_mapping.items():
        uid = next((u for e, u in uniprot_ids if e == entity_id), None)
        for chain_id in chain_mapping.get(entity_id, ["A"]):
            key = (pdb_code.lower(), acc, chain_id, entity_id)
            if key in seen:
                continue
            seen.add(key)
            table_data.append({
                'PDB': pdb_code.lower(), 'Entity': entity_id, 'Chain': chain_id,
                'UniProt AC': acc, 'UniProt ID': uid,
                'HETATM': hetatm_names, 'Resolution': resolution, 'Method': method,
            })
    for w in check_chimera(table_data):
        st.warning("⚠️ **Chimera warning!!**")
        st.write(f"🔴 {w}")
    return table_data

def build_simplified(all_data):
    grouped = defaultdict(list)
    for row in all_data:
        key = (row['PDB'], row['Entity'], row['UniProt AC'], row['UniProt ID'],
               row['HETATM'], row['Resolution'], row['Method'])
        grouped[key].append(row['Chain'])
    result = [{'PDB': k[0], 'Entity': k[1], 'Chain': '/'.join(sorted(v)),
               'UniProt AC': k[2], 'UniProt ID': k[3], 'HETATM': k[4],
               'Resolution': k[5], 'Method': k[6]} for k, v in grouped.items()]
    return sorted(result, key=lambda x: (x['PDB'], x['Entity']))

progress_bar = st.progress(0)
status_text = st.empty()

if fetch_button:
    pdb_codes = [c for c in re.split(r'[,\s]+', pdb_input.strip().upper()) if c]
    st.session_state.list_pdb_codes = pdb_codes
    all_data = []
    with tempfile.TemporaryDirectory() as tmp_dir:
        for i, pdb_code in enumerate(pdb_codes, 1):
            status_text.text(f"Analysing {pdb_code} ({i}/{len(pdb_codes)})...")
            all_data.extend(process_cif_file(pdb_code, tmp_dir))
            progress_bar.progress(i / len(pdb_codes))
    all_data = [dict(t) for t in {tuple(d.items()) for d in all_data}]
    st.session_state.list_all_data = all_data

status_text.empty()
progress_bar.empty()

all_data = st.session_state.list_all_data
pdb_codes = st.session_state.list_pdb_codes

if all_data:
    df = pd.DataFrame(build_simplified(all_data))

    st.markdown("<a name='simplified-data'></a>", unsafe_allow_html=True)
    st.markdown("#### Proteins and ligands")
    st.dataframe(alternate_colors(df), hide_index=True, use_container_width=True)

    # Fix #8: get_gene_name and get_protein_name are now cached in utils
    st.markdown("<a name='uniprot-ac'></a>", unsafe_allow_html=True)
    st.markdown("#### UniProt list")
    unique_acs = sorted({row['UniProt AC'] for row in all_data if row['UniProt AC']})
    st.code(" ".join(unique_acs))
    for ac in unique_acs:
        gene = get_gene_name_from_uniprot(ac)
        name = get_protein_name_from_uniprot(ac)
        st.write(f"[{ac}](https://www.uniprot.org/uniprotkb/{ac}) {gene}: {name}")

    st.markdown("<a name='ligands'></a>", unsafe_allow_html=True)
    hetatms = sorted({item.strip()
                      for sublist in df['HETATM'].str.split(', ')
                      for item in sublist if item.strip()})
    if hetatms:
        st.markdown("#### Ligands / Modifications")
        with st.spinner("Fetching ligand information..."):
            ligand_names = {h: get_ligand_name(h) for h in hetatms}
        html = "".join(
            f"[{h}](https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/{h}): "
            f"{ligand_names.get(h, '') or '(Ligand name not available)'}<br>"
            for h in hetatms
        )
        st.markdown(html, unsafe_allow_html=True)

    st.markdown("<a name='pdb'></a>", unsafe_allow_html=True)
    st.markdown("#### PDB list")
    st.code(" ".join(p.lower() for p in pdb_codes))
    for pdb_code in pdb_codes:
        title = get_structure_title(pdb_code)
        st.write(f"[{pdb_code}](https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_code.lower()}): {title}")

    st.markdown("<a name='pymol-commands'></a>", unsafe_allow_html=True)
    st.markdown("#### PyMOL Commands")
    st.code("fetch " + "; fetch ".join(p.lower() for p in pdb_codes))
    st.code(generate_pymol_commands(all_data))

    if st.button("Other PyMOL Commands"):
        st.switch_page("pages/Pymol_commands.py")
