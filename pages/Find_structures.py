import streamlit as st
import pandas as pd
import requests
import re

# Sidebar, title, parameters
# -----------------------------------------------------------------------------
# layout="wide" is not working without scrollbar
st.set_page_config(page_title="Haku - Find the structures", page_icon="💮")

default_unp = 'Q9NR09'

# Initialize session state
with st.sidebar:
    st.header("📝 New query")
    uniprot_ac = st.text_input("Enter UniProt AC:", default_unp).strip()
    fetch_button = st.button("Fetch data")
    if fetch_button:
        st.session_state.guide_structure = False
    st.sidebar.markdown("[PDB and PMID](#pdb-and-pmid)")
    st.sidebar.markdown("[Structures in PDBe](#structures)")
    st.sidebar.markdown("[Associated literature](#associated-literature-pmids)")
    st.sidebar.markdown("[PyMOL commands](#pymol-fetch-command)")
    guide_structure = st.checkbox("Instructions", value=st.session_state.get("guide_structure", True))
    st.session_state.guide_structure = guide_structure

st.markdown("#### Find all the structures of a protein using its UniProt AC.")

# Instructions
# -----------------------------------------------------------------------------
if st.session_state.get('guide_structure'):
    st.info("""
    Input:&emsp;**[UniProt AC](https://www.uniprot.org/)**  
      
    Output:
    * all associated structures (PDBs)
    * associated literature (PMIDs)
    * PyMOL command to download the structures
    
    """)

# Define functions
# -----------------------------------------------------------------------------

# Function to check if a PMID is curated in UniProt
def is_pmid_curated_in_uniprot(uniprot_ac, pmid):
    if not pmid or pmid == "Not available":
        return False
    uniprot_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}?format=txt"
    response = requests.get(uniprot_url)
    if response.status_code != 200:
        print(f"Error fetching data from UniProt for AC {uniprot_ac}")
        return False
    uniprot_text = response.text
    first_line = uniprot_text.split('\n')[0]
    if "Unreviewed;" in first_line:
        return False  # Not a SwissProt entry, ignore PMID check
    return pmid in uniprot_text

# Function to fetch protein name from UniProt
def get_protein_name_from_uniprot(uniprot_ac):
    uniprot_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}?format=txt"
    response = requests.get(uniprot_url)
    if response.status_code != 200:
        st.error(f"Error fetching data from UniProt for AC {uniprot_ac}")
        return None
    uniprot_text = response.text
    # Extract the protein name from the response
    for line in uniprot_text.splitlines():
        if line.startswith("DE   RecName: Full="):
            # Extract the part before any curly braces
            protein_name = line.split("=")[1].strip()
            protein_name_cleaned = re.sub(r"\{.*?\}", "", protein_name).strip()
            protein_name_cleaned = re.sub(r"[;.,]+$", "", protein_name_cleaned).strip()
            return protein_name_cleaned
    return "Protein name not available."

# Function to fetch PDB codes from PDBe
def get_pdb_codes_from_pdbe(uniprot_ac):
    uniprot_ac = uniprot_ac.upper()
    pdbe_search_url = f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{uniprot_ac}"
    response = requests.get(pdbe_search_url)
    if response.status_code != 200:
        print(f"\033[1mError fetching data from PDBe for AC {uniprot_ac}\033[0m")
        return []
    data = response.json()
    pdb_codes = set()
    if uniprot_ac in data:
        for entry in data[uniprot_ac]:
            if "pdb_id" in entry:
                pdb_codes.add(entry["pdb_id"].lower())
    if not pdb_codes:
        print(f"No PDB codes found for UniProt AC {uniprot_ac} in PDBe.")
    return list(pdb_codes)

# Function to fetch PDB codes from UniProt
def get_pdb_codes_from_uniprot(uniprot_ac):
    uniprot_url = f"https://rest.uniprot.org/uniprot/{uniprot_ac}.xml"
    response = requests.get(uniprot_url)
    if response.status_code != 200:
        st.error(f"Error fetching data from UniProt for AC {uniprot_ac}")
        return []
    xml_data = response.text
    pdb_codes = {line.split('id="')[1].split('"')[0].lower() for line in xml_data.splitlines() if 'dbReference type="PDB"' in line}
    return list(pdb_codes)

# Function to convert DOI to PMID using Europe PMC API
def convert_doi_to_pmid(doi):
    if not doi:
        return None
    europe_pmc_url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:{doi}&resultType=core&format=json"
    response = requests.get(europe_pmc_url)
    if response.status_code != 200:
        return None
    data = response.json()
    results = data.get('resultList', {}).get('result', [])
    return results[0].get('pmid') if results else None

# Function to fetch PMIDs from PDBe and convert DOI to PMID if necessary
def get_pdb_pmids_from_pdbe(pdb_code):
    pdbe_url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/publications/{pdb_code}"
    response = requests.get(pdbe_url)
    if response.status_code != 200:
        return []
    pdb_data = response.json()
    pmids = set()
    if pdb_code in pdb_data:
        for publication in pdb_data[pdb_code]:
            pmid = publication.get('pubmed_id') or convert_doi_to_pmid(publication.get('doi'))
            if pmid:
                pmids.add(pmid)
    return list(pmids)

# Function to fetch publication info from EuropePMC
def get_info_from_PMC(pmid):
    if not pmid:
        return None
    europe_pmc_url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={pmid}&format=json"
    response = requests.get(europe_pmc_url)
    if response.status_code != 200:
        return None
    data = response.json()
    results = data.get('resultList', {}).get('result', [])
    if not results:
        return None
    # Find the result with the exact PMID match
    for result in results:
        if result.get("pmid") == str(pmid):
            return {
                "PMID": result.get("pmid"),
                "Title": result.get("title"),
                "Year": result.get("pubYear"),
                "Journal": result.get("journalTitle"),
                "Authors": result.get("authorString")
            }
    return None

# Function to generate dataframe on info from EuropePMC
def generate_pmid_dataframe(all_pmids):
    data = [get_info_from_PMC(pmid) for pmid in all_pmids if get_info_from_PMC(pmid)]
    return pd.DataFrame(data)

## The script begins here - Data Fetching
# -----------------------------------------------------------------------------

if fetch_button and uniprot_ac:
    progress_bar = st.progress(0)
    progress_status = 0
    total_steps = 5  # Define total steps for smooth progress updates

    protein_name = get_protein_name_from_uniprot(uniprot_ac)
    if protein_name:
        st.session_state.protein_name = protein_name
    progress_status += 1
    progress_bar.progress(int(progress_status / total_steps * 100))
    
    pdb_codes_pdbe = set(get_pdb_codes_from_pdbe(uniprot_ac))
    progress_status += 1
    progress_bar.progress(int(progress_status / total_steps * 100))
    
    pdb_codes_uniprot = set(get_pdb_codes_from_uniprot(uniprot_ac))
    progress_status += 1
    progress_bar.progress(int(progress_status / total_steps * 100))
    
    if not pdb_codes_pdbe:
        st.warning("No PDB structure available for this UniProt entry.")
        progress_bar.empty()
        st.session_state.data = None
        st.session_state.all_pmids = None
    else:
        pdb_not_in_uniprot = pdb_codes_pdbe - pdb_codes_uniprot
        pdb_in_uniprot = pdb_codes_pdbe & pdb_codes_uniprot

        data = []
        all_pmids = set()

        for idx, pdb_code in enumerate(sorted(pdb_codes_pdbe)):
            pmids = get_pdb_pmids_from_pdbe(pdb_code)
            pmid_str = ", ".join(pmids) if pmids else "Not available"
            available = "Yes" if pdb_code in pdb_codes_uniprot else "No"
            data.append({"PDB": pdb_code.upper(), "PMID": pmid_str, "Available in UniProt": available})
            all_pmids.update(pmids)
            progress_bar.progress(int((progress_status + (idx + 1) / len(pdb_codes_pdbe)) / total_steps * 100))

        progress_status += 1  # Completed PDB processing

        # Add the Curated in UniProt column
        curated_column = []
        for row in data:
            pmids = row["PMID"].split(", ") if row["PMID"] != "Not available" else []
            curated_results = [("Yes" if is_pmid_curated_in_uniprot(uniprot_ac, pmid) else "No") for pmid in pmids]
            curated_column.append(", ".join(curated_results) if curated_results else "Not available")

        for idx, row in enumerate(data):
            row["Curated in UniProt"] = curated_column[idx]
            progress_bar.progress(int((progress_status + (idx + 1) / len(data)) / total_steps * 100))
        
        progress_status += 1  # Completed curation check

        # Ensure progress reaches 100%
        progress_bar.progress(100)
        progress_bar.empty()
        
        # Store data in session state for persistence
        st.session_state.data = data
        st.session_state.all_pmids = all_pmids

## Display Logic - Shows persisted data
# -----------------------------------------------------------------------------

if hasattr(st.session_state, 'data') and st.session_state.data is not None:
    data = st.session_state.data
    all_pmids = st.session_state.all_pmids
    
    # Display protein name if available
    if hasattr(st.session_state, 'protein_name'):
        st.markdown(f"**Protein Name**: {st.session_state.protein_name}")
    
    # Functions to highlight rows in df tables
    def highlight_not_available(row):
        return ['background-color: #fffdcd' if row['Available in UniProt'] == 'No' else '' for _ in row]
    def highlight_curated(row):
        return ['background-color: #fffede' if row['Curated in UniProt'] == 'Yes' else '' for _ in row]

    # Section 1: PDB and PMID Mapping Table
    # -----------------------------------------------------------------------------
    st.markdown("<a name='pdb-and-pmid'></a>", unsafe_allow_html=True)
    st.markdown("#### Structures (PDB) and literature (PMID)")
    df = pd.DataFrame(data)
    df = df.sort_values(by="PDB", ascending=False).style.hide(axis='index').apply(highlight_not_available, axis=1)
    st.dataframe(df, hide_index=True, use_container_width=True)

    # Section 2: PDBs Available in UniProt
    # -----------------------------------------------------------------------------
    st.markdown("<a name='structures'></a>", unsafe_allow_html=True)
    st.markdown("#### Structures in PDBe")
    available_pdbs = [entry["PDB"] for entry in data if entry["Available in UniProt"] == "Yes"]
    unavailable_pdbs = [entry["PDB"] for entry in data if entry["Available in UniProt"] == "No"]
    col1, col2 = st.columns([2, 1])
    with col1:
        if available_pdbs:
            available_links = ", ".join([f"[{pdb_code}](https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_code.lower()})" for pdb_code in available_pdbs])
            st.markdown(f"Available in UniProt:<br>{available_links}", unsafe_allow_html=True)
        if unavailable_pdbs:
            unavailable_links = ", ".join([f"[{pdb_code}](https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_code.lower()})" for pdb_code in unavailable_pdbs])
            st.markdown(f"Not yet available in UniProt:<br>{unavailable_links}", unsafe_allow_html=True)
    with col2:
        if available_pdbs:
            pdb_str = " ".join(p.lower() for p in available_pdbs)
            st.code(pdb_str)
        if unavailable_pdbs:
            pdb_str = " ".join(p.lower() for p in unavailable_pdbs)
            st.code(pdb_str)

    # Section 3: Associated Literature (PMIDs)
    # -----------------------------------------------------------------------------
    st.markdown("<a name='associated-literature-pmids'></a>", unsafe_allow_html=True)
    st.markdown("#### Associated literature (PMIDs)")
    curated_pmids = [entry["PMID"] for entry in data if entry["Curated in UniProt"] == "Yes"]
    if all_pmids:
        pmid_df = pd.DataFrame(generate_pmid_dataframe(all_pmids))
        pmid_df["Curated in UniProt"] = pmid_df["PMID"].apply(lambda pmid: "Yes" if pmid in curated_pmids else "No")
        pmid_df = pmid_df.sort_values(by="PMID", ascending=False).style.hide(axis='index').apply(highlight_curated, axis=1)
        html_table = pmid_df.to_html(escape=False, index=False)
        st.dataframe(pmid_df, hide_index=True, use_container_width=True)
        
        pmid_links = ", ".join([f"[{pmid}](https://europepmc.org/abstract/MED/{pmid})" for pmid in sorted(all_pmids, reverse=True)])
        st.markdown(f"Europe PMC: {pmid_links}")
        
    # Section 4: PyMOL Fetch Command
    # -----------------------------------------------------------------------------
    st.markdown("<a name='pymol-fetch-command'></a>", unsafe_allow_html=True)
    if available_pdbs:
        st.markdown("#### PyMOL Fetch Command")
        st.code("; ".join(f"fetch {p.lower()}" for p in available_pdbs))

    # Section 5: Coder's Note
    # -----------------------------------------------------------------------------
    st.markdown(
        """
        ---
All the structures (PDB) mapped to a UniProt entry are listed, using the UniProt and PDBe API.  
This includes structures recently released by PDB but are not yet available in UniProt.  
Literature associated to each structure (PMID) is listed using the Europe PMC API.  
Higlighted entries are curated in UniProt.  
        """
    )