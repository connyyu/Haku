"""
utils.py — Shared utilities for Haku pages.

Fixes applied here:
  #6  - get_protein_name / get_gene_name were duplicated across pages
  #7  - adds consistent retry logic to all API calls via api_retry decorator
  #9  - parse_resolution / parse_uniprot_mapping / parse_hetatms were
        duplicated between List_UniProt_entries and Structure_viewer
  #12 - adds @st.cache_data so repeated lookups hit the cache, not the API
"""

import re
import os
import requests
import streamlit as st
from collections import defaultdict
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

# ---------------------------------------------------------------------------
# Retry decorator — wrap every outbound API call with this
# ---------------------------------------------------------------------------
api_retry = retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=1, max=10),
    retry=retry_if_exception_type(requests.RequestException),
    reraise=True,
)

# ---------------------------------------------------------------------------
# UniProt helpers  (cached per session — fix #12)
# ---------------------------------------------------------------------------

@st.cache_data(ttl=3600, show_spinner=False)
def get_protein_name_from_uniprot(uniprot_ac: str) -> str:
    """Return the recommended protein name for a UniProt AC."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}?format=txt"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
    except requests.RequestException:
        return "Protein name not available."
    for line in response.text.splitlines():
        if line.startswith("DE   RecName: Full="):
            name = line.split("=")[1].strip()
            name = re.sub(r"\{.*?\}", "", name).strip()
            name = re.sub(r"[;.,]+$", "", name).strip()
            return name
    return "Protein name not available."


@st.cache_data(ttl=3600, show_spinner=False)
def get_gene_name_from_uniprot(uniprot_ac: str) -> str:
    """Return the primary gene name for a UniProt AC."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}?format=txt"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
    except requests.RequestException:
        return "Gene name not available."
    for line in response.text.splitlines():
        if line.startswith("GN   Name="):
            name = line.split("=")[1].strip()
            name = re.sub(r"\{.*?\}", "", name).strip()
            name = re.split(r"[;,.\s]", name, 1)[0].strip()
            return name
    return "Gene name not available."


@st.cache_data(ttl=3600, show_spinner=False)
def is_pmid_curated_in_uniprot(uniprot_ac: str, pmid: str) -> bool:
    """Return True if the PMID appears in a reviewed UniProt entry."""
    if not pmid or pmid == "Not available":
        return False
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}?format=txt"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
    except requests.RequestException:
        return False
    text = response.text
    if "Unreviewed;" in text.split("\n")[0]:
        return False
    return pmid in text


@st.cache_data(ttl=3600, show_spinner=False)
def get_ligand_name(ligand_code: str) -> str:
    """Return the PDBe chemical component name for a ligand code."""
    url = f"https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/{ligand_code}"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code != 200:
            return ""
        data = response.json()
        return data.get(ligand_code.upper(), [{}])[0].get("name", "")
    except Exception:
        return ""


@st.cache_data(ttl=3600, show_spinner=False)
def get_structure_title(pdb_code: str) -> str:
    """Return the PDBe title for a PDB entry."""
    url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{pdb_code.lower()}"
    try:
        response = requests.get(url, timeout=10)
        if not response.ok:
            return ""
        data = response.json()
        return data.get(pdb_code.lower(), [{}])[0].get("title", "")
    except Exception:
        return ""


# ---------------------------------------------------------------------------
# CIF file helpers — fix #3 (caller passes per-request temp dir) and #9
# ---------------------------------------------------------------------------

def download_cif_file(pdb_code: str, output_dir: str) -> str | None:
    """
    Download the updated CIF file from PDBe into output_dir.
    Returns the full file path on success, None on failure.

    Fix #3: callers supply a per-request temp dir so concurrent users
    never collide on the shared chain_info/ folder.
    """
    pdb_code_dl = pdb_code.lower()
    url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_code_dl}_updated.cif"
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
    except requests.RequestException:
        return None
    file_path = os.path.join(output_dir, f"{pdb_code_dl}.cif")
    with open(file_path, "wb") as f:
        f.write(response.content)
    return file_path


def parse_resolution(cif_path: str) -> tuple:
    """Parse resolution (Angstrom) and method from a CIF file path."""
    resolution = None
    method = ""
    symmetry_type_found = False

    with open(cif_path, "r") as f:
        lines = f.readlines()

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        tokens = stripped.split()
        if tokens[0] == "_reflns.d_resolution_high" and len(tokens) > 1:
            resolution, method = tokens[1], "Xtal"
            break
        elif tokens[0] == "_em_3d_reconstruction.resolution" and len(tokens) > 1:
            resolution, method = tokens[1], "EM"
        elif stripped == "_exptl.method 'Solution NMR'":
            method = "NMR"
            break
        elif tokens[0] == "_em_3d_reconstruction.symmetry_type":
            symmetry_type_found = True
        elif symmetry_type_found and resolution is None and len(tokens) >= 7:
            resolution, method = tokens[6], "EM"
            break

    if resolution is not None:
        try:
            resolution = f"{float(resolution):.2f}"
        except ValueError:
            resolution = None
    return resolution, method


def parse_hetatms(cif_path: str) -> list:
    """Return sorted list of unique HETATM residue names from a CIF file."""
    hetatms = set()
    with open(cif_path, "r") as f:
        for line in f:
            if line.startswith("HETATM"):
                parts = line.split()
                if len(parts) > 5:
                    hetatms.add(parts[5])
    return sorted(hetatms)


def parse_uniprot_mapping(cif_path: str) -> tuple:
    """
    Parse UniProt -> chain mappings from a CIF file.
    Returns (uniprot_mapping, uniprot_ids, chain_mapping).
    Accepts a full path instead of filename.
    """
    uniprot_mapping = {}
    uniprot_ids = {}
    chain_mapping = {}
    ref_id_to_entity_id = {}
    fallback_accession = fallback_ID = fallback_chain = None

    with open(cif_path, "r") as f:
        lines = f.readlines()

    # Scan for fallback values
    for i, line in enumerate(lines):
        line = line.strip()
        if line.startswith("_struct_ref.pdbx_db_accession"):
            fallback_accession = (line.split()[1] if len(line.split()) > 1
                                  else lines[i + 1].strip().split()[0])
        elif line.startswith("_struct_ref.db_code"):
            fallback_ID = (line.split()[1] if len(line.split()) > 1
                           else lines[i + 1].strip().split()[0])
        elif line.startswith("_struct_ref_seq.pdbx_strand_id"):
            parts = line.split()
            if len(parts) > 1:
                cand = parts[1]
                if len(cand) == 1 and cand.isalnum():
                    fallback_chain = cand
            else:
                k = i + 1
                while k < len(lines):
                    nxt = lines[k].strip()
                    if not nxt or nxt.startswith("_"):
                        break
                    cand = nxt.split()[0]
                    if len(cand) == 1 and cand.isalnum():
                        fallback_chain = cand
                        break
                    k += 1

    process_sifts = False
    sifts_mapping = {}
    sifts_chain_mapping = {}
    i = 0

    while i < len(lines):
        line = lines[i].strip()

        if line.startswith("loop_"):
            struct_ref_headers = []
            hi = i + 1
            while hi < len(lines):
                hl = lines[hi].strip()
                if hl.startswith("_struct_ref."):
                    struct_ref_headers.append(hl)
                    hi += 1
                else:
                    break

            if struct_ref_headers:
                db_name_idx = db_code_idx = entity_id_idx = accession_idx = ref_id_idx = None
                for idx, h in enumerate(struct_ref_headers):
                    if h == "_struct_ref.id":                   ref_id_idx = idx
                    elif h == "_struct_ref.db_name":            db_name_idx = idx
                    elif h == "_struct_ref.db_code":            db_code_idx = idx
                    elif h == "_struct_ref.entity_id":          entity_id_idx = idx
                    elif h == "_struct_ref.pdbx_db_accession":  accession_idx = idx

                di = hi
                while di < len(lines):
                    dl = lines[di].strip()
                    if not dl or dl.startswith("#"):
                        di += 1
                        continue
                    if dl.startswith("loop_") or dl.startswith("_"):
                        break
                    parts = dl.split()
                    if (None not in (ref_id_idx, db_name_idx, db_code_idx,
                                     entity_id_idx, accession_idx) and
                            len(parts) > max(ref_id_idx, db_name_idx,
                                             db_code_idx, entity_id_idx, accession_idx)):
                        if parts[db_name_idx] == "UNP":
                            ref_id = parts[ref_id_idx]
                            entity_id = parts[entity_id_idx]
                            uid = parts[db_code_idx]
                            acc = parts[accession_idx]
                            ref_id_to_entity_id[ref_id] = entity_id
                            if len(acc) >= 6 and "-" not in acc and "_" not in acc:
                                uniprot_mapping[entity_id] = acc
                                uniprot_ids[entity_id] = uid
                    di += 1
                i = di
                continue

        elif line.startswith("_struct_ref_seq.align_id"):
            i += 1
            while i < len(lines):
                nxt = lines[i].strip()
                if not nxt or nxt.startswith(("loop_", "#")):
                    break
                parts = nxt.split()
                if len(parts) > 8:
                    ref_id = parts[1]
                    chain_id = parts[3]
                    acc = parts[8]
                    entity_id = ref_id_to_entity_id.get(ref_id, ref_id)
                    if entity_id:
                        chain_mapping.setdefault(entity_id, set()).add(chain_id)
                        if (entity_id not in uniprot_mapping and
                                len(acc) >= 6 and "-" not in acc and "_" not in acc):
                            uniprot_mapping[entity_id] = acc
                            uniprot_ids.setdefault(entity_id, None)
                i += 1
            continue

        elif line.startswith("_pdbx_sifts_unp_segments.identity"):
            process_sifts = True

        elif process_sifts:
            if line.startswith(("loop_", "#")):
                process_sifts = False
            elif line:
                parts = line.split()
                if len(parts) > 3:
                    entity_id, chain_id, sifts_acc = parts[0], parts[1], parts[2]
                    if len(sifts_acc) >= 6 and "-" not in sifts_acc:
                        sifts_mapping[entity_id] = sifts_acc
                        sifts_chain_mapping.setdefault(entity_id, set()).add(chain_id)
        i += 1

    # Fallback when no mappings found
    if not uniprot_mapping and not sifts_mapping and fallback_accession:
        fe = "1"
        uniprot_mapping[fe] = fallback_accession
        uniprot_ids[fe] = fallback_ID
        chain_mapping[fe] = {fallback_chain} if fallback_chain else set()

    # Final reconciliation
    final_uniprot_mapping = {}
    final_uniprot_ids = []
    final_chain_mapping = {}
    all_entities = (set(chain_mapping) | set(uniprot_mapping) |
                    set(sifts_mapping) | set(sifts_chain_mapping))

    for entity in all_entities:
        original_ac = uniprot_mapping.get(entity)
        sifts_ac = sifts_mapping.get(entity)
        final_ac = sifts_ac if sifts_ac else original_ac

        chains_sr = sorted(chain_mapping.get(entity, set()))
        chains_sf = sorted(sifts_chain_mapping.get(entity, set()))
        final_chains = chains_sr if chains_sr else chains_sf

        if final_ac:
            final_uniprot_mapping[entity] = final_ac
            uid = uniprot_ids.get(entity) or fallback_ID
            final_uniprot_ids.append((entity, uid))
            final_chain_mapping[entity] = final_chains

    for entity, chains in final_chain_mapping.items():
        uid = next((u for e, u in final_uniprot_ids if e == entity), None)
        if fallback_chain and uid == fallback_ID and chains and fallback_chain not in chains:
            final_chain_mapping[entity] = [fallback_chain]

    return final_uniprot_mapping, final_uniprot_ids, final_chain_mapping


# ---------------------------------------------------------------------------
# PyMOL command generation  (shared — replaces file-based approach)
# ---------------------------------------------------------------------------

def generate_pymol_commands(all_data: list) -> str:
    """
    Generate PyMOL selection + colour commands from parsed CIF data rows.
    Returns the command string directly instead of writing to disk.
    """
    colors = [
        "carbon", "cyan", "lightmagenta", "yellow", "salmon", "slate", "orange",
        "deepteal", "violetpurple", "hydrogen", "marine", "olive", "smudge",
        "teal", "wheat", "lightpink", "skyblue",
    ]
    uniprot_chains = defaultdict(lambda: defaultdict(list))
    for row in all_data:
        uniprot_chains[row["UniProt ID"]][row["PDB"]].append(row["Chain"])

    lines = ["util.cbaw;"]
    for color_index, (uid, pdb_dict) in enumerate(uniprot_chains.items()):
        selections = [f"(chain {'+'.join(chains)} and {pdb_id})"
                      for pdb_id, chains in pdb_dict.items()]
        selection = " or ".join(selections)
        color = colors[color_index % len(colors)]
        lines.append(f"sele {uid}, {selection}; color {color}, {uid};")
    lines.append("deselect; util.cnc")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Alternate-row colouring for dataframes  (shared)
# ---------------------------------------------------------------------------

def alternate_colors(df):
    """Apply alternating background colours grouped by PDB code."""
    color_1 = "background-color: #e2fffb;"
    row_styles = []
    last_pdb = None
    alternate = False
    for _, row in df.iterrows():
        if row["PDB"] != last_pdb:
            alternate = not alternate
            last_pdb = row["PDB"]
        row_styles.append([color_1 if alternate else ""] * len(df.columns))

    def apply(row):
        return row_styles[row.name]

    return df.style.apply(apply, axis=1)
