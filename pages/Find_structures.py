"""
Find_structures.py — Find all structures for a UniProt AC.

Changes from v1:
  - All API calls moved to utils.py (cached, retried)
  - Fixed double PMC call per PMID (#2) — get_publication_info now cached
  - Session state keys prefixed fs_ to avoid cross-page clashes (#1)
"""

import streamlit as st
import pandas as pd
import requests

from utils import (
    get_protein_name_from_uniprot,
    is_pmid_curated_in_uniprot,
)

# ---------------------------------------------------------------------------
# Local helpers that aren't shared across pages (kept here, not in utils)
# ---------------------------------------------------------------------------

@st.cache_data(ttl=3600, show_spinner=False)
def get_pdb_codes_from_pdbe(uniprot_ac: str) -> list:
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{uniprot_ac.upper()}"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code != 200:
            return []
        data = r.json()
        return sorted({e["pdb_id"].lower() for e in data.get(uniprot_ac.upper(), []) if "pdb_id" in e})
    except Exception:
        return []


@st.cache_data(ttl=3600, show_spinner=False)
def get_pdb_codes_from_uniprot(uniprot_ac: str) -> list:
    url = f"https://rest.uniprot.org/uniprot/{uniprot_ac}.xml"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code != 200:
            return []
        return sorted({
            line.split('id="')[1].split('"')[0].lower()
            for line in r.text.splitlines()
            if 'dbReference type="PDB"' in line
        })
    except Exception:
        return []


@st.cache_data(ttl=3600, show_spinner=False)
def _doi_to_pmid(doi: str) -> str | None:
    if not doi:
        return None
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:{doi}&resultType=core&format=json"
    try:
        r = requests.get(url, timeout=10)
        results = r.json().get("resultList", {}).get("result", [])
        return results[0].get("pmid") if results else None
    except Exception:
        return None


@st.cache_data(ttl=3600, show_spinner=False)
def get_pdb_pmids(pdb_code: str) -> list:
    url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/publications/{pdb_code.lower()}"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code != 200:
            return []
        pmids = set()
        for pub in r.json().get(pdb_code.lower(), []):
            pmid = pub.get("pubmed_id") or _doi_to_pmid(pub.get("doi", ""))
            if pmid:
                pmids.add(str(pmid))
        return sorted(pmids)
    except Exception:
        return []


@st.cache_data(ttl=3600, show_spinner=False)
def get_publication_info(pmid: str) -> dict | None:
    """Cached — each PMID fetched at most once per session. Fixes issue #2."""
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={pmid}&format=json"
    try:
        r = requests.get(url, timeout=10)
        results = r.json().get("resultList", {}).get("result", [])
        for result in results:
            if result.get("pmid") == str(pmid):
                return {
                    "PMID": result.get("pmid"),
                    "Title": result.get("title"),
                    "Year": result.get("pubYear"),
                    "Journal": result.get("journalTitle"),
                    "Authors": result.get("authorString"),
                }
    except Exception:
        pass
    return None


# ---------------------------------------------------------------------------
# Page config & sidebar
# ---------------------------------------------------------------------------
st.set_page_config(page_title="Haku - Find the structures", page_icon="💮")

default_unp = "Q9NR09"

with st.sidebar:
    st.header("📝 New query")
    uniprot_ac = st.text_input("Enter UniProt AC:", default_unp).strip()
    fetch_button = st.button("Fetch data")
    if fetch_button:
        st.session_state.fs_guide = False
    st.sidebar.markdown("[PDB and PMID](#pdb-and-pmid)")
    st.sidebar.markdown("[Structures in PDBe](#structures)")
    st.sidebar.markdown("[Associated literature](#associated-literature-pmids)")
    st.sidebar.markdown("[PyMOL commands](#pymol-fetch-command)")
    st.sidebar.markdown("[3D viewer ↗](https://haku3dviewer.netlify.app/)")
    fs_guide = st.checkbox("Instructions", value=st.session_state.get("fs_guide", True))
    st.session_state.fs_guide = fs_guide

st.markdown("#### Find all the structures of a protein using its UniProt AC.")

if st.session_state.get("fs_guide"):
    st.info("""
    Input:&emsp;**[UniProt AC](https://www.uniprot.org/)**  
      
    Output:
    * all associated structures (PDBs)
    * associated literature (PMIDs)
    * PyMOL command to download the structures
    
    Open **[Haku 3D viewer](https://haku3dviewer.netlify.app/)** for visualisation in Mol*.
    """)

# ---------------------------------------------------------------------------
# Data fetch
# ---------------------------------------------------------------------------
if fetch_button and uniprot_ac:
    progress_bar = st.progress(0)

    protein_name = get_protein_name_from_uniprot(uniprot_ac)
    st.session_state.fs_protein_name = protein_name
    progress_bar.progress(20)

    pdb_codes_pdbe = set(get_pdb_codes_from_pdbe(uniprot_ac))
    progress_bar.progress(40)

    pdb_codes_uniprot = set(get_pdb_codes_from_uniprot(uniprot_ac))
    progress_bar.progress(60)

    if not pdb_codes_pdbe:
        st.warning("No PDB structure available for this UniProt entry.")
        progress_bar.empty()
        st.session_state.fs_data = None
        st.session_state.fs_all_pmids = None
    else:
        data = []
        all_pmids = set()

        for idx, pdb_code in enumerate(sorted(pdb_codes_pdbe)):
            pmids = get_pdb_pmids(pdb_code)
            pmid_str = ", ".join(pmids) if pmids else "Not available"
            available = "Yes" if pdb_code in pdb_codes_uniprot else "No"
            data.append({
                "PDB": pdb_code.upper(),
                "PMID": pmid_str,
                "Available in UniProt": available,
            })
            all_pmids.update(pmids)
            progress_bar.progress(60 + int(30 * (idx + 1) / len(pdb_codes_pdbe)))

        # Curation column — reuses cached UniProt text via is_pmid_curated_in_uniprot
        for row in data:
            pmids_list = [p for p in row["PMID"].split(", ") if p != "Not available"]
            curated = ["Yes" if is_pmid_curated_in_uniprot(uniprot_ac, p) else "No"
                       for p in pmids_list]
            row["Curated in UniProt"] = ", ".join(curated) if curated else "Not available"

        progress_bar.progress(100)
        progress_bar.empty()

        st.session_state.fs_data = data
        st.session_state.fs_all_pmids = all_pmids

# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------
if st.session_state.get("fs_data"):
    data = st.session_state.fs_data
    all_pmids = st.session_state.fs_all_pmids

    if st.session_state.get("fs_protein_name"):
        st.markdown(f"**Protein Name**: {st.session_state.fs_protein_name}")

    def highlight_not_available(row):
        return ["background-color: #fffdcd" if row["Available in UniProt"] == "No" else "" for _ in row]

    def highlight_curated(row):
        return ["background-color: #fffede" if row["Curated in UniProt"] == "Yes" else "" for _ in row]

    # Section 1: PDB + PMID table
    st.markdown("<a name='pdb-and-pmid'></a>", unsafe_allow_html=True)
    st.markdown("#### Structures (PDB) and literature (PMID)")
    df = pd.DataFrame(data).sort_values("PDB", ascending=False)
    st.dataframe(df.style.apply(highlight_not_available, axis=1), hide_index=True, use_container_width=True)

    # Section 2: PDB links
    st.markdown("<a name='structures'></a>", unsafe_allow_html=True)
    st.markdown("#### Structures in PDBe")
    available_pdbs = [r["PDB"] for r in data if r["Available in UniProt"] == "Yes"]
    unavailable_pdbs = [r["PDB"] for r in data if r["Available in UniProt"] == "No"]
    col1, col2 = st.columns([2, 1])
    with col1:
        if available_pdbs:
            links = ", ".join(f"[{p}](https://www.ebi.ac.uk/pdbe/entry/pdb/{p.lower()})" for p in available_pdbs)
            st.markdown(f"Available in UniProt:<br>{links}", unsafe_allow_html=True)
        if unavailable_pdbs:
            links = ", ".join(f"[{p}](https://www.ebi.ac.uk/pdbe/entry/pdb/{p.lower()})" for p in unavailable_pdbs)
            st.markdown(f"Not yet available in UniProt:<br>{links}", unsafe_allow_html=True)
    with col2:
        if available_pdbs:
            st.code(" ".join(p.lower() for p in available_pdbs))
        if unavailable_pdbs:
            st.code(" ".join(p.lower() for p in unavailable_pdbs))

    # Section 3: Literature
    st.markdown("<a name='associated-literature-pmids'></a>", unsafe_allow_html=True)
    st.markdown("#### Associated literature (PMIDs)")
    curated_pmids = {r["PMID"] for r in data if r.get("Curated in UniProt") == "Yes"}
    if all_pmids:
        pub_rows = [get_publication_info(p) for p in all_pmids]
        pub_rows = [r for r in pub_rows if r]
        if pub_rows:
            pmid_df = pd.DataFrame(pub_rows)
            pmid_df["Curated in UniProt"] = pmid_df["PMID"].apply(
                lambda p: "Yes" if p in curated_pmids else "No"
            )
            pmid_df = pmid_df.sort_values("PMID", ascending=False)
            st.dataframe(pmid_df.style.apply(highlight_curated, axis=1),
                         hide_index=True, use_container_width=True)
        pmid_links = ", ".join(
            f"[{p}](https://europepmc.org/abstract/MED/{p})"
            for p in sorted(all_pmids, reverse=True)
        )
        st.markdown(f"Europe PMC: {pmid_links}")

    # Section 4: PyMOL
    st.markdown("<a name='pymol-fetch-command'></a>", unsafe_allow_html=True)
    if available_pdbs:
        st.markdown("#### PyMOL Fetch Command")
        st.code("; ".join(f"fetch {p.lower()}" for p in available_pdbs))

    st.markdown("""
        ---
All the structures (PDB) mapped to a UniProt entry are listed, using the UniProt and PDBe API.  
This includes structures recently released by PDB but are not yet available in UniProt.  
Literature associated to each structure (PMID) is listed using the Europe PMC API.  
Highlighted entries are curated in UniProt.  
    """)
