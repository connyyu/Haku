### some clash with Find_structures.py using the same parameters, to fix.

import streamlit as st
import pandas as pd
import requests
import re
import xml.etree.ElementTree as ET

# Sidebar, title, parameters
# -----------------------------------------------------------------------------
# layout="wide" is not working without scrollbar
st.set_page_config(page_title="Haku - Similar structures", page_icon="💮")

default_unp = 'Q13002'

# Initialize session state
with st.sidebar:
    st.header("📝 New query")
    uniprot_ac_query = st.text_input("Enter UniProt AC:", default_unp).strip()
    fetch_button = st.button("Fetch data")
    if fetch_button:
        st.session_state.guide_homolog = False
    st.sidebar.markdown("[Homologs](#homologs)")
    st.sidebar.markdown("[UniRef90](#uniref90)")
    st.sidebar.markdown("[UniRef50](#uniref50)")
    guide_homolog = st.checkbox("Instructions", value=st.session_state.get("guide_homolog", True))
    st.session_state.guide_homolog = guide_homolog

if "fetched" not in st.session_state:
    st.session_state.fetched = False
    
st.markdown("#### Find similar proteins with experimental structures.")

# Instructions
# -----------------------------------------------------------------------------
if st.session_state.get('guide_homolog'):
    st.info("""
    Input:&emsp;**[UniProt AC](https://www.uniprot.org/)**  
      
    Output:
    * Homologs (entries with the same UniProt ID prefix)
    * UniProt entries with experimentally determined structures from UniRef90 and UniRef50 clusters
    """)

# Define functions
# -----------------------------------------------------------------------------
# Fetch UniRef clusters from UniProt
def get_clusters_from_uniprot(uniprot_ac):
    search_url = f"https://rest.uniprot.org/uniref/search?query={uniprot_ac}&format=json&size=500"
    response = requests.get(search_url)
    if response.status_code != 200:
        st.error(f"Error fetching data from UniProt for AC {uniprot_ac}")
        return [], [], []
    data = response.json()
    results = data.get("results", [])
    uniref100 = [entry["id"] for entry in results if entry["id"].startswith("UniRef100_")]
    uniref90 = [entry["id"] for entry in results if entry["id"].startswith("UniRef90_")]
    uniref50 = [entry["id"] for entry in results if entry["id"].startswith("UniRef50_")]
    return uniref100, uniref90, uniref50

# Function to fetch UniRef clusters entries with structures
def get_clusters_entries_with_structures(cluster_ids, uniprot_ac):
    if not cluster_ids:
        return []

    cluster_with_structure = []
    for cluster_id in cluster_ids:
        if cluster_id == "N/A":
            continue
        # Determine cluster type dynamically
        if cluster_id.startswith("UniRef90_"):
            cluster_type = "90"
        elif cluster_id.startswith("UniRef50_"):
            cluster_type = "50"
        elif cluster_id.startswith("UniRef100_"):
            cluster_type = "100"
        else:
            st.warning(f"Unknown cluster type for {cluster_id}, skipping")
            continue
        kb_query = f"(uniref_cluster_{cluster_type}:{cluster_id}) NOT (accession:{uniprot_ac}) AND structure_3d=true"
        kb_url = f"https://rest.uniprot.org/uniprotkb/search?query={kb_query}&format=json&size=500"
        kb_resp = requests.get(kb_url)
        if kb_resp.status_code != 200:
            st.error(f"Error fetching data from UniRef{cluster_type} cluster: {cluster_id}")
            continue
        kb_data = kb_resp.json()
        for entry in kb_data.get("results", []):
            entry['cluster_type'] = cluster_type
            entry['cluster_id'] = cluster_id
            cluster_with_structure.append(entry)
    return cluster_with_structure

# Function to fetch gene name from UniProt
def get_gene_name_from_uniprot(uniprot_ac):
    uniprot_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}?format=txt"
    response = requests.get(uniprot_url)
    if response.status_code != 200:
        st.error(f"Error fetching data from UniProt for AC {uniprot_ac}")
        return None
    uniprot_text = response.text
    for line in uniprot_text.splitlines():
        if line.startswith("GN   Name="):
            gene_name = line.split("=")[1].strip()
            gene_name_cleaned = re.sub(r"\{.*?\}", "", gene_name).strip()
            gene_name_cleaned = re.split(r"[;,.\s]", gene_name, 1)[0].strip()
            return gene_name_cleaned
    return "Gene name not available."

# Function to fetch protein name from UniProt
def get_protein_name_from_uniprot(uniprot_ac):
    uniprot_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}?format=txt"
    response = requests.get(uniprot_url)
    if response.status_code != 200:
        st.error(f"Error fetching data from UniProt for AC {uniprot_ac}")
        return None
    uniprot_text = response.text
    for line in uniprot_text.splitlines():
        if line.startswith("DE   RecName: Full="):
            protein_name = line.split("=")[1].strip()
            protein_name_cleaned = re.sub(r"\{.*?\}", "", protein_name).strip()
            protein_name_cleaned = re.sub(r"[;.,]+$", "", protein_name_cleaned).strip()
            return protein_name_cleaned
    return "Protein name not available."

# Functions to fetch taxonomy from UniProt (to display mammalian entries first)
def get_uniprot_xml(uniprot_ac):
    url = f"https://rest.uniprot.org/uniprot/{uniprot_ac}.xml"
    resp = requests.get(url)
    if resp.status_code != 200:
        return None
    return resp.text

def is_mammal(uniprot_ac):
    """Check if a UniProt entry belongs to mammals by searching the lineage"""
    xml_data = get_uniprot_xml(uniprot_ac)
    if not xml_data:
        return False
    try:
        root = ET.fromstring(xml_data)
        namespaces = {'uniprot': 'http://uniprot.org/uniprot'}
        taxons = root.findall('.//uniprot:lineage/uniprot:taxon', namespaces)
        for taxon in taxons:
            if taxon.text == "Mammalia":
                return True
        return False
    except ET.ParseError:
        return False

# Function to fetch PDB codes from cluster with structures
def get_pdb_codes_from_cluster_entries(cluster_with_structure):
    results_list = []
    seen = set()
    for entry in cluster_with_structure:
        ac = entry.get("primaryAccession")
        if not ac or (ac, entry['cluster_type']) in seen:
            continue
        seen.add((ac, entry['cluster_type']))
        pdb_ids = [xref["id"] for xref in entry.get("uniProtKBCrossReferences", []) if xref.get("database") == "PDB"]
        if not pdb_ids:
            continue
        cluster_type = entry.get('cluster_type', 'N/A')
        cluster_id = entry.get('cluster_id', 'N/A')
        results_list.append({
            "Cluster Type": f"UniRef{cluster_type}",
            "Cluster ID": cluster_id,
            "UniProt AC": ac,
            "Protein Name": entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "N/A"),
            "Organism": entry.get("organism", {}).get("scientificName", "N/A"),
            "PDB IDs": ", ".join(pdb_ids)
        })
    return results_list

# Function to fetch PDB codes from UniProt
def get_pdb_codes_from_uniprot(uniprot_ac):
    uniprot_url = f"https://rest.uniprot.org/uniprot/{uniprot_ac}.xml"
    response = requests.get(uniprot_url)
    if response.status_code != 200:
        st.error(f"Error fetching data from UniProt for AC {uniprot_ac}")
        return []
    xml_data = response.text
    pdb_codes = {line.split('id="')[1].split('"')[0].upper() for line in xml_data.splitlines() if 'dbReference type="PDB"' in line}
    return list(pdb_codes)

# Function to display each cluster type
def display_cluster(cluster_prefix, title):
    cluster_df = df[df["Cluster Type"] == cluster_prefix]
    if not cluster_df.empty:
        st.markdown(f"<a name='{cluster_prefix.lower()}'></a>", unsafe_allow_html=True)
        st.markdown(f"#### {title}")
        cluster_df["num_pdb"] = cluster_df["PDB IDs"].apply(
            lambda x: len(x.split(",")) if isinstance(x, str) and x.strip() else 0
        )
        cluster_df = cluster_df.sort_values(by="num_pdb", ascending=False).drop(columns=["Cluster Type", "num_pdb"])
        st.dataframe(cluster_df, hide_index=True, use_container_width=True)

        st.markdown(f"##### UniProt entries")
        unique_acs = sorted(set(cluster_df['UniProt AC'].dropna()))
        for ac in unique_acs:
            gene_name = get_gene_name_from_uniprot(ac)
            protein_name = get_protein_name_from_uniprot(ac)
            st.write(f"[{ac}](https://www.uniprot.org/uniprotkb/{ac}) {gene_name}: {protein_name}")

        st.markdown(f"##### PDB entries")
        pdb_codes = sorted(set(",".join(cluster_df['PDB IDs'].dropna()).split(",")))
        pdb_links = ", ".join([
            f"[{pdb_code.strip()}](https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_code.strip().lower()})"
            for pdb_code in pdb_codes
        ])
        st.markdown(pdb_links, unsafe_allow_html=True)

def get_uniprot_id_from_ac(uniprot_ac):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}.json"
    resp = requests.get(url)
    if resp.status_code != 200:
        st.warning(f"Could not fetch UniProt ID for {uniprot_ac}")
        return None
    data = resp.json()
    return data.get("uniProtkbId")

def get_id_prefix_from_id(uniprot_id):
    if "_" in uniprot_id:
        return uniprot_id.rsplit("_", 1)[0] + "_"
    return uniprot_id + "_"

def search_homologs_by_id_prefix(uniprot_ac):
    uniprot_id = get_uniprot_id_from_ac(uniprot_ac)
    if not uniprot_id:
        return []
    prefix = get_id_prefix_from_id(uniprot_id)
    kb_query = f"id:{prefix}* AND NOT accession:{uniprot_ac}"
    kb_url = f"https://rest.uniprot.org/uniprotkb/search?query={kb_query}&format=json&size=500"
    resp = requests.get(kb_url)
    if resp.status_code != 200:
        st.warning(f"Error fetching homologs for prefix {prefix}")
        return []

    data = resp.json()
    results = data.get("results", [])
    homologs = []
    for entry in results:
        ac = entry.get("primaryAccession")
        uid = entry.get("uniProtkbId")
        if ac and uid and ac != uniprot_ac:
            # Filter by same ID prefix only
            if uid.startswith(prefix):
                homologs.append({"UniProt_AC": ac, "UniProt_ID": uid})
    return homologs

def get_homologs_structures(homologs):
    results = []
    for h in homologs:
        ac = h["UniProt_AC"]
        uid = h["UniProt_ID"]
        protein_name = get_protein_name_from_uniprot(ac)
        organism = "N/A"
        pdb_ids = []
        try:
            url = f"https://rest.uniprot.org/uniprotkb/{ac}.json"
            resp = requests.get(url)
            if resp.status_code == 200:
                data = resp.json()
                organism = data.get("organism", {}).get("scientificName", "N/A")
                for xref in data.get("uniProtKBCrossReferences", []):
                    if xref.get("database") == "PDB":
                        pdb_ids.append(xref.get("id").upper())
        except:
            pass
        results.append({
            "UniProt AC": ac,
            "UniProt ID": uid,
            "Protein Name": protein_name,
            "Organism": organism,
            "PDB IDs": ", ".join(pdb_ids)
        })
    return results

## The script begins here - Data Fetching
# -----------------------------------------------------------------------------

if fetch_button and uniprot_ac_query:
    st.session_state.fetched = True
    progress_bar = st.progress(0)
    progress_status = 0
    total_steps = 4

    # Step 1: Fetch Protein Name
    protein_name = get_protein_name_from_uniprot(uniprot_ac_query)
    st.session_state.protein_name_homolog = protein_name if protein_name else "N/A"
    progress_status += 1
    progress_bar.progress(int(progress_status / total_steps * 100))

    # Step 2: Fetch PDB codes for UniProt entry
    pdb_codes_direct = get_pdb_codes_from_uniprot(uniprot_ac_query)
    st.session_state.pdb_codes_direct = pdb_codes_direct
    progress_status += 1
    progress_bar.progress(int(progress_status / total_steps * 100))

    # Step 3: Fetch homologs with the same UniProt ID prefix
    homologs = search_homologs_by_id_prefix(uniprot_ac_query)
    st.session_state.homologs_table = get_homologs_structures(homologs) if homologs else []

    # Step 4: Fetch UniRef clusters
    _, uniref90, uniref50 = get_clusters_from_uniprot(uniprot_ac_query)

    # Flatten and convert all items to strings
    def string_only(cluster_list):
        flat_list = []
        for item in cluster_list:
            if isinstance(item, list) or isinstance(item, set):
                flat_list.extend(str(i) for i in item)
            else:
                flat_list.append(str(item))
        return flat_list
    uniref90 = string_only(uniref90)
    uniref50 = string_only(uniref50)
        
    if not uniref90 and not uniref50:
        st.warning("No UniRef clusters found for this UniProt entry.")
        st.session_state.data_homolog = []
        progress_bar.empty()
    else:

    # Step 5: Fetch cluster members with structures
        clusters_combined = uniref90 + uniref50
        cluster_with_structure = get_clusters_entries_with_structures(clusters_combined, uniprot_ac_query)
        progress_status += 1
        progress_bar.progress(int(progress_status / total_steps * 100))

        # Step 6: Fetch PDB codes from cluster entries
        results_list = get_pdb_codes_from_cluster_entries(cluster_with_structure)
        progress_status += 1
        progress_bar.progress(int(progress_status / total_steps * 100))
        st.session_state.data_homolog = results_list

        # Ensure progress reaches 100%
        progress_bar.progress(100)
        progress_bar.empty()

## Data display
# Section 1: About the query entry
# -----------------------------------------------------------------------------
if st.session_state.fetched:
    st.markdown(f"**Protein Name:** {st.session_state.protein_name_homolog}")

# Show direct experimental structures (PDBs for the UniProt entry itself)
if hasattr(st.session_state, 'pdb_codes_direct') and st.session_state.pdb_codes_direct:
    pdb_codes_sorted = sorted(st.session_state.pdb_codes_direct)
    pdb_links = ", ".join([
        f"[{pdb_code}](https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_code.lower()})" for pdb_code in pdb_codes_sorted
    ])
    st.markdown(f"**Available Structures:** {pdb_links}", unsafe_allow_html=True)

# Section 2: Look for entries with the same prefix
# -----------------------------------------------------------------------------
if st.session_state.fetched and hasattr(st.session_state, 'homologs_table') and st.session_state.homologs_table:
    st.markdown("<a name='homologs'></a>", unsafe_allow_html=True)
    st.markdown("#### Homologs")

    df_homologs = pd.DataFrame(st.session_state.homologs_table)
    df_homologs["mammal_priority"] = df_homologs["UniProt AC"].apply(
        lambda ac: 0 if is_mammal(ac) else 1
    )

    # Sort mammals first, then by UniProt AC
    df_homologs = df_homologs.sort_values(
        by=["mammal_priority", "UniProt AC"],
        ascending=[True, True]
    ).drop(columns=["mammal_priority"])
    st.dataframe(df_homologs, hide_index=True, use_container_width=True)

    st.markdown("##### UniProt entries")
    ordered_acs = df_homologs["UniProt AC"].tolist()
    ac_links = ", ".join(
        f"[{ac}](https://www.uniprot.org/uniprotkb/{ac})"
        for ac in ordered_acs
    )
    st.write(ac_links)

# Section 3: UniRef clusters
if st.session_state.get("fetched", False):
    data = st.session_state.get("data_homolog", [])

    if data:
        df = pd.DataFrame(data)
        display_cluster("UniRef90", "UniRef90 Cluster")
        display_cluster("UniRef50", "UniRef50 Cluster")
    else:
        st.info("No homolog structures found in UniRef clusters.")
        
# Section 4: Coder's Note
# -----------------------------------------------------------------------------
st.markdown(
    """
    ---
This script uses **UniProt Reference Clusters (UniRef)** to identify similar proteins. Details of these clusters can be found on the [UniProt help page](https://www.uniprot.org/help/uniref).

*UniRef100* - identical sequences and sub-fragments with > 11 residues  
*UniRef90* - sequences with > 90 % sequence identity and > 80 % overlap  
*UniRef50* - sequences with > 50 % sequence identity and > 80 % overlap  

    """
)
