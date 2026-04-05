import streamlit as st
import pandas as pd
import requests
import re
import xml.etree.ElementTree as ET

st.set_page_config(page_title="Haku - Similar structures", page_icon="💮")

default_unp = 'Q13002'

# Fix #1: use unambiguous session state keys (homolog_*) so they
# never clash with Find_structures.py which uses uniprot_ac / protein_name
with st.sidebar:
    st.header("📝 New query")
    homolog_uniprot_ac = st.text_input("Enter UniProt AC:", default_unp).strip()
    fetch_button = st.button("Fetch data")
    if fetch_button:
        st.session_state.guide_homolog = False
    st.sidebar.markdown("[Homologs](#homologs)")
    st.sidebar.markdown("[UniRef90](#uniref90)")
    st.sidebar.markdown("[UniRef50](#uniref50)")
    guide_homolog = st.checkbox("Instructions", value=st.session_state.get("guide_homolog", True))
    st.session_state.guide_homolog = guide_homolog

if "homolog_fetched" not in st.session_state:
    st.session_state.homolog_fetched = False

st.markdown("#### Find similar proteins with experimental structures.")

if st.session_state.get('guide_homolog'):
    st.info("""
    Input:&emsp;**[UniProt AC](https://www.uniprot.org/)**\n\n
    Output:
    * Homologs (entries with the same UniProt ID prefix)
    * UniProt entries with experimentally determined structures from UniRef90 and UniRef50 clusters
    """)

# Functions
# -----------------------------------------------------------------------------

@st.cache_data(ttl=3600, show_spinner=False)
def get_clusters_from_uniprot(uniprot_ac):
    url = f"https://rest.uniprot.org/uniref/search?query={uniprot_ac}&format=json&size=500"
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
    except requests.RequestException:
        st.error(f"Error fetching clusters for {uniprot_ac}")
        return [], [], []
    results = r.json().get("results", [])
    uniref100 = [e["id"] for e in results if e["id"].startswith("UniRef100_")]
    uniref90  = [e["id"] for e in results if e["id"].startswith("UniRef90_")]
    uniref50  = [e["id"] for e in results if e["id"].startswith("UniRef50_")]
    return uniref100, uniref90, uniref50


@st.cache_data(ttl=3600, show_spinner=False)
def get_clusters_entries_with_structures(cluster_ids, uniprot_ac):
    if not cluster_ids:
        return []
    cluster_with_structure = []
    for cluster_id in cluster_ids:
        if cluster_id == "N/A":
            continue
        if cluster_id.startswith("UniRef90_"):
            cluster_type = "90"
        elif cluster_id.startswith("UniRef50_"):
            cluster_type = "50"
        elif cluster_id.startswith("UniRef100_"):
            cluster_type = "100"
        else:
            continue
        kb_query = (f"(uniref_cluster_{cluster_type}:{cluster_id}) "
                    f"NOT (accession:{uniprot_ac}) AND structure_3d=true")
        url = f"https://rest.uniprot.org/uniprotkb/search?query={kb_query}&format=json&size=500"
        try:
            r = requests.get(url, timeout=10)
            r.raise_for_status()
        except requests.RequestException:
            continue
        for entry in r.json().get("results", []):
            ac = entry.get("primaryAccession", "")
            uid = entry.get("uniProtkbId", "")
            name = ""
            for desc in entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", ""):
                name = entry["proteinDescription"]["recommendedName"]["fullName"]["value"]
                break
            if not name:
                try:
                    name = entry["proteinDescription"]["recommendedName"]["fullName"]["value"]
                except (KeyError, TypeError):
                    name = ""
            cluster_with_structure.append({
                "UniProt AC": ac,
                "UniProt ID": uid,
                "Protein Name": name,
                "Cluster": cluster_id,
            })
    return cluster_with_structure


@st.cache_data(ttl=3600, show_spinner=False)
def get_homologs(uniprot_ac):
    # Homologs share the same UniProt ID prefix (before the underscore)
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}.json"
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
    except requests.RequestException:
        return []
    data = r.json()
    uid = data.get("uniProtkbId", "")
    if not uid or "_" not in uid:
        return []
    gene_prefix = uid.split("_")[0]
    search_url = f"https://rest.uniprot.org/uniprotkb/search?query=id:{gene_prefix}_*+AND+structure_3d:true&format=json&size=100"
    try:
        r2 = requests.get(search_url, timeout=10)
        r2.raise_for_status()
    except requests.RequestException:
        return []
    homologs = []
    for entry in r2.json().get("results", []):
        ac = entry.get("primaryAccession", "")
        if ac == uniprot_ac:
            continue
        homologs.append({
            "UniProt AC": ac,
            "UniProt ID": entry.get("uniProtkbId", ""),
        })
    return homologs


# Main
# -----------------------------------------------------------------------------

if fetch_button:
    st.session_state.homolog_fetched = True
    st.session_state.homolog_query_ac = homolog_uniprot_ac

if st.session_state.homolog_fetched:
    query_ac = st.session_state.get("homolog_query_ac", homolog_uniprot_ac)
    progress = st.progress(0)

    uniref100, uniref90, uniref50 = get_clusters_from_uniprot(query_ac)
    progress.progress(25)

    # Homologs
    st.markdown("<a name='homologs'></a>", unsafe_allow_html=True)
    st.markdown("#### Homologs")
    homologs = get_homologs(query_ac)
    if homologs:
        df_h = pd.DataFrame(homologs)
        df_h["UniProt AC"] = df_h["UniProt AC"].apply(
            lambda ac: f"[{ac}](https://www.uniprot.org/uniprotkb/{ac})"
        )
        st.dataframe(df_h, hide_index=True, use_container_width=True)
    else:
        st.info("No homologs with structures found.")
    progress.progress(50)

    # UniRef90
    st.markdown("<a name='uniref90'></a>", unsafe_allow_html=True)
    st.markdown("#### UniRef90 cluster members with structures")
    entries90 = get_clusters_entries_with_structures(uniref90, query_ac)
    progress.progress(75)
    if entries90:
        df90 = pd.DataFrame(entries90)
        st.dataframe(df90, hide_index=True, use_container_width=True)
    else:
        st.info("No UniRef90 cluster members with structures found.")

    # UniRef50
    st.markdown("<a name='uniref50'></a>", unsafe_allow_html=True)
    st.markdown("#### UniRef50 cluster members with structures")
    entries50 = get_clusters_entries_with_structures(uniref50, query_ac)
    progress.progress(100)
    progress.empty()
    if entries50:
        df50 = pd.DataFrame(entries50)
        st.dataframe(df50, hide_index=True, use_container_width=True)
    else:
        st.info("No UniRef50 cluster members with structures found.")
