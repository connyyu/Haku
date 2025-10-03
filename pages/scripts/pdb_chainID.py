import csv
import re
import os
import sys
import glob
import requests
from collections import defaultdict

script_dir = os.path.dirname(os.path.abspath(__file__))  # location of pages directory

def download_cif_file(pdb_code):
    pdb_code_dl = pdb_code.lower()
    url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_code_dl}_updated.cif"
    response = requests.get(url)
    if response.status_code == 200:
        with open(f"{pdb_code_dl}.cif", 'wb') as f:
            f.write(response.content)
        return f"{pdb_code_dl}.cif"
    else:
        print(f"\033[1mFailed to download CIF file for {pdb_code_dl}\033[0m")
        return None

def parse_resolution(cif_file):
    resolution = None
    method = ""
    symmetry_type_found = False

    with open(cif_file, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if line.startswith('_reflns.d_resolution_high'):
            parts = line.split()
            if len(parts) > 1:
                resolution = parts[1]
                method = "Xtal"
                break
        elif line.startswith('_em_3d_reconstruction.resolution'):
            parts = line.split()
            if len(parts) > 1:
                resolution = parts[1]
                method = "EM"
                break
        elif line.startswith("_exptl.method 'Solution NMR'"):
            method = "NMR"
            break
        elif line.startswith('_em_3d_reconstruction.symmetry_type'):
            symmetry_type_found = True
        elif symmetry_type_found:
            parts = line.split()
            if len(parts) >= 7:
                resolution = parts[6]
                method = "EM"
            break
    if resolution is not None:
        try:
            resolution = f"{float(resolution):.2f}"
        except ValueError:
            resolution = None

    return resolution, method

def parse_hetatms(cif_file):
    hetatms = set()
    with open(cif_file, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('HETATM'):
            parts = line.split()
            hetatm_name = parts[5]
            hetatms.add(hetatm_name)
    return sorted(hetatms)

def fetch_uniprot_id(uniprot_ac):
    url = f'https://rest.uniprot.org/uniprotkb/{uniprot_ac}.json'
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        return data.get('uniProtkbId', None)
    return None

def parse_uniprot_mapping(cif_file):
    uniprot_mapping = {}
    uniprot_ids = {}
    chain_mapping = {}
    ref_id_to_entity_id = {}

    fallback_accession = None
    fallback_ID = None
    fallback_chain = None

    with open(cif_file, 'r') as f:
        lines = f.readlines()

    # Scan for fallbacks
    for i, line in enumerate(lines):
        line = line.strip()
        if line.startswith('_struct_ref.pdbx_db_accession'):
            fallback_accession = (line.split()[1] if len(line.split()) > 1
                                  else lines[i + 1].strip().split()[0])
        elif line.startswith('_struct_ref.db_code'):
            fallback_ID = (line.split()[1] if len(line.split()) > 1
                           else lines[i + 1].strip().split()[0])
        elif line.startswith('_struct_ref_seq.pdbx_strand_id'):
            parts = line.split()
            if len(parts) > 1:
                fallback_chain_candidate = parts[1]
                if len(fallback_chain_candidate) == 1 and fallback_chain_candidate.isalnum():
                    fallback_chain = fallback_chain_candidate
            else:
                k = i + 1
                while k < len(lines):
                    next_line = lines[k].strip()
                    if not next_line or next_line.startswith('_'):
                        break
                    fallback_chain_candidate = next_line.split()[0]
                    if len(fallback_chain_candidate) == 1 and fallback_chain_candidate.isalnum():
                        fallback_chain = fallback_chain_candidate
                        break
                    k += 1

    i = 0
    process_sifts = False
    sifts_mapping = {}
    sifts_chain_mapping = {}
    struct_ref_count = 0
    struct_ref_seq_count = 0
    sifts_count = 0

    while i < len(lines):
        line = lines[i].strip()

        if line.startswith('loop_'):
            struct_ref_headers = []
            header_line_idx = i + 1
            while header_line_idx < len(lines):
                header_line = lines[header_line_idx].strip()
                if header_line.startswith('_struct_ref.'):
                    struct_ref_headers.append(header_line)
                    header_line_idx += 1
                elif header_line.startswith('_') and not header_line.startswith('_struct_ref.'):
                    break
                else:
                    break

            if struct_ref_headers:
                db_name_idx = db_code_idx = entity_id_idx = accession_idx = ref_id_idx = None

                for idx, header in enumerate(struct_ref_headers):
                    if header == '_struct_ref.id':
                        ref_id_idx = idx
                    elif header == '_struct_ref.db_name':
                        db_name_idx = idx
                    elif header == '_struct_ref.db_code':
                        db_code_idx = idx
                    elif header == '_struct_ref.entity_id':
                        entity_id_idx = idx
                    elif header == '_struct_ref.pdbx_db_accession':
                        accession_idx = idx

                data_line_idx = header_line_idx
                while data_line_idx < len(lines):
                    data_line = lines[data_line_idx].strip()
                    if not data_line or data_line.startswith('#'):
                        data_line_idx += 1
                        continue
                    if data_line.startswith('loop_') or data_line.startswith('_'):
                        break

                    parts = data_line.split()
                    struct_ref_count += 1

                    if (None not in (ref_id_idx, db_name_idx, db_code_idx, entity_id_idx, accession_idx) and
                        len(parts) > max(ref_id_idx, db_name_idx, db_code_idx, entity_id_idx, accession_idx)):
                        
                        if parts[db_name_idx] == 'UNP':
                            ref_id = parts[ref_id_idx]
                            entity_id = parts[entity_id_idx]
                            uniprot_id = parts[db_code_idx]
                            uniprot_accession = parts[accession_idx]

                            ref_id_to_entity_id[ref_id] = entity_id

                            if len(uniprot_accession) >= 6 and '-' not in uniprot_accession and '_' not in uniprot_accession:
                                uniprot_mapping[entity_id] = uniprot_accession
                                uniprot_ids[entity_id] = uniprot_id  # Store as dict for easier lookup

                    data_line_idx += 1
                i = data_line_idx
                continue

        elif line.startswith('_struct_ref_seq.align_id'):
            i += 1
            while i < len(lines):
                next_line = lines[i].strip()
                if not next_line or next_line.startswith(('loop_', '#')):
                    break
                parts = next_line.split()
                if len(parts) > 8:
                    struct_ref_seq_count += 1
                    ref_id = parts[1]
                    chain_id = parts[3]
                    uniprot_accession = parts[8]

                    entity_id = ref_id_to_entity_id.get(ref_id)
                    if entity_id is None:
                        entity_id = ref_id

                    if entity_id:
                        chain_mapping.setdefault(entity_id, set()).add(chain_id)
                        if (entity_id not in uniprot_mapping and
                            len(uniprot_accession) >= 6 and 
                            '-' not in uniprot_accession and 
                            '_' not in uniprot_accession):
                            uniprot_mapping[entity_id] = uniprot_accession
                            if entity_id not in uniprot_ids:
                                uniprot_ids[entity_id] = None
                i += 1
            continue

        elif line.startswith('_pdbx_sifts_unp_segments.identity'):
            process_sifts = True

        elif process_sifts:
            if line.startswith(('loop_', '#')):
                process_sifts = False
            elif line:
                parts = line.split()
                if len(parts) > 3:
                    sifts_count += 1
                    entity_id, chain_id, sifts_accession = parts[0], parts[1], parts[2]
                    if len(sifts_accession) >= 6 and '-' not in sifts_accession:
                        sifts_mapping[entity_id] = sifts_accession
                        sifts_chain_mapping.setdefault(entity_id, set()).add(chain_id)

        i += 1

    # Use fallback if no mapping found
    if not uniprot_mapping and not sifts_mapping:
        fallback_entity = "1"

        if fallback_accession:
            uniprot_mapping[fallback_entity] = fallback_accession
            uniprot_ids[fallback_entity] = fallback_ID

            chains_for_entity = chain_mapping.get(fallback_entity)
            if chains_for_entity:
                chain_mapping[fallback_entity] = chains_for_entity
            else:
                chain_mapping[fallback_entity] = set([fallback_chain] if fallback_chain else [])

    # Final reconciliation - Use SIFTS chains only if struct_ref_seq chains are NOT available
    final_uniprot_mapping = {}
    final_uniprot_ids = []
    final_chain_mapping = {}

    all_entities = set(chain_mapping.keys()) | set(uniprot_mapping.keys()) | set(sifts_mapping.keys()) | set(sifts_chain_mapping.keys())
    
    for entity in all_entities:
        original_ac = uniprot_mapping.get(entity, None)
        uniprot_id = uniprot_ids.get(entity, None)
        sifts_ac = sifts_mapping.get(entity, None)

        final_ac = sifts_ac if sifts_ac else original_ac

        chains_from_struct_ref_seq = list(chain_mapping.get(entity, set()))
        chains_from_sifts = list(sifts_chain_mapping.get(entity, set()))
        
        if chains_from_struct_ref_seq:
            final_chains = sorted(chains_from_struct_ref_seq)
        elif chains_from_sifts:
            final_chains = sorted(chains_from_sifts)
        else:
            final_chains = []
        
        if final_ac:
            final_uniprot_mapping[entity] = final_ac

            uid = uniprot_ids.get(entity)
            if uid is None and fallback_ID:
                uid = fallback_ID
        
            final_uniprot_ids.append((entity, uid))
            final_chain_mapping[entity] = final_chains

    for entity, chains in final_chain_mapping.items():
        uid = None
        for ent, uniprot_id in final_uniprot_ids:
            if ent == entity:
                uid = uniprot_id
                break
        if fallback_chain and uid == fallback_ID:
            if chains and (fallback_chain not in chains):
                final_chain_mapping[entity] = [fallback_chain]

    return final_uniprot_mapping, final_uniprot_ids, final_chain_mapping

def check_chimera(table_data):
    chain_uniprot_map = defaultdict(set)
    for row in table_data:
        pdb_code = row['PDB'].upper()
        chain_id = row['Chain ID']
        uniprot_accession = row['UniProt AC']
        chain_uniprot_map[chain_id].add(uniprot_accession)
    warnings = []
    for chain_id, uniprot_accessions in chain_uniprot_map.items():
        if len(uniprot_accessions) > 1:
            warnings.append(f"{pdb_code} chain {chain_id} is associated with mutiple UniProt AC: {', '.join(uniprot_accessions)}.")
    return warnings

def process_cif_file(pdb_code):
    cif_file = download_cif_file(pdb_code)
    if not cif_file:
        print(f"\033[1mFailed to download CIF file for {pdb_code}. Skipping...\033[0m")
        return None

    resolution, method = parse_resolution(cif_file)
    hetatms = parse_hetatms(cif_file)
    uniprot_mapping, uniprot_ids, chain_mapping = parse_uniprot_mapping(cif_file)

    hetatm_names = ', '.join(hetatms)
    table_data = []

    printed_chains = set()
    for entity_id, uniprot_accession in uniprot_mapping.items():
        uniprot_id = next((id for ent_id, id in uniprot_ids if ent_id == entity_id), None)
        chain_ids = chain_mapping.get(entity_id, "A")
        for chain_id in chain_ids:
            if (uniprot_accession, chain_id) in printed_chains:
                continue
            printed_chains.add((uniprot_accession, chain_id))
            table_data.append({
                'PDB': pdb_code.lower(),
                'Entity': entity_id,
                'Chain ID': chain_id,
                'UniProt AC': uniprot_accession,
                'UniProt ID': uniprot_id,
                'HETATM': hetatm_names,
                'Resolution': resolution,
                'Method': method
            })

    warnings = check_chimera(table_data)

    return table_data, warnings

def get_chain_ids(pdb_code, uniprot_ac=None):
    """
    Main function to retrieve chain IDs and UniProt info for a PDB code.
    Optionally filter results for a specific UniProt accession.
    Returns list of dictionaries of chain data and any warnings.
    """
    table_data, warnings = process_cif_file(pdb_code)
    if not table_data:
        return [], []

    if uniprot_ac:
        filtered = [row for row in table_data if row['UniProt AC'] == uniprot_ac]
    else:
        filtered = table_data

    return filtered, warnings

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <PDB_CODE> [UNIPROT_AC]")
        sys.exit(1)
    pdb_code = sys.argv[1]
    uniprot_ac = sys.argv[2] if len(sys.argv) > 2 else None
    results, warnings = get_chain_ids(pdb_code, uniprot_ac)

    for row in results:
        print(f"PDB: {row['PDB']}\tEntity: {row['Entity']}\tChain ID: {row['Chain ID']}\t"
              f"UniProt AC: {row['UniProt AC']}\tUniProt ID: {row['UniProt ID']}\t"
              f"HETATM: {row['HETATM']}\tResolution: {row['Resolution']}\tMethod: {row['Method']}")

    if warnings:
        print("\nWarnings:")
        for w in warnings:
            print(w)

if __name__ == "__main__":
    main()
