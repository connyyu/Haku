import re
import os
import sys
import glob
import csv
import requests
from collections import defaultdict
import tempfile
import shutil

def download_cif_file(pdb_code, output_dir=None):
    """Download CIF file for given PDB code."""
    if output_dir is None:
        output_dir = tempfile.gettempdir()
    
    pdb_code_dl = pdb_code.lower()
    url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_code_dl}_updated.cif"
    response = requests.get(url)
    if response.status_code == 200:
        cif_path = os.path.join(output_dir, f"{pdb_code_dl}.cif")
        with open(cif_path, 'wb') as f:
            f.write(response.content)        
        return cif_path
    else:
        return None

def parse_resolution(cif_file):
    """Parse resolution and method from CIF file."""
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
        elif symmetry_type_found:  # Look at the next line after symmetry_type
            parts = line.split()
            if len(parts) >= 7:
                resolution = parts[6]  # Take the 7th column as resolution
                method = "EM"
            break
    # 2 decimal places for resolution
    if resolution is not None:
        try:
            resolution = f"{float(resolution):.2f}"
        except ValueError:
            resolution = None

    return resolution, method

def parse_hetatms(cif_file):
    """Parse HETATM records from CIF file."""
    hetatms = set()
    with open(cif_file, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('HETATM'):
            parts = line.split()
            hetatm_name = parts[5]  # HETATM name (e.g., GLC, ZN)
            hetatms.add(hetatm_name)  # Use a set to avoid duplicates
    return sorted(hetatms)  # Return a sorted list

def fetch_uniprot_id(uniprot_ac):
    """Fetch UniProt ID from accession."""
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
    """Check for chimeric chains (chains with multiple UniProt accessions)."""
    chain_uniprot_map = defaultdict(set)
    pdb_code = ""
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
    
def process_cif_file(pdb_code, output_dir=None):
    """Process CIF file and extract data."""
    cif_file = download_cif_file(pdb_code, output_dir)
    if not cif_file:
        return None, []
    
    try:
        resolution, method = parse_resolution(cif_file)    
        hetatms = parse_hetatms(cif_file)
        uniprot_mapping, uniprot_ids, chain_mapping = parse_uniprot_mapping(cif_file)

        hetatm_names = ', '.join(hetatms)
        table_data = []
        
        printed_chains = set()
        for entity_id, uniprot_accession in uniprot_mapping.items():
            uniprot_id = next((id for ent_id, id in uniprot_ids if ent_id == entity_id), None)
            chain_ids = chain_mapping.get(entity_id, ["A"])
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

        # Check for warnings
        warnings = check_chimera(table_data)
        
        return table_data, warnings
    
    finally:
        # Clean up temporary CIF file
        if cif_file and os.path.exists(cif_file):
            try:
                os.remove(cif_file)
            except:
                pass  # Ignore cleanup errors

def pdb_unpID(pdb_code, chain_id, verbose=False, return_all_data=False):
    """
    Main function to get UniProt IDs for a given PDB code and chain ID.
    
    Args:
        pdb_code (str): PDB code (e.g., '1ABC')
        chain_id (str): Chain identifier (e.g., 'A')
        verbose (bool): If True, print warnings and additional information
        return_all_data (bool): If True, return all parsed data along with UniProt IDs
    
    Returns:
        If return_all_data is False:
            list: Sorted list of UniProt IDs for the specified chain
        If return_all_data is True:
            tuple: (uniprot_ids_list, all_table_data, warnings)
    """
    # Create temporary directory for processing
    temp_dir = tempfile.mkdtemp()
    
    try:
        # Process the CIF file
        table_data, warnings = process_cif_file(pdb_code, temp_dir)
        
        if not table_data:
            if verbose:
                print(f"No data found for PDB code {pdb_code}")
            if return_all_data:
                return [], [], []
            return []

        # Print warnings if verbose
        if warnings and verbose:
            print("WARNINGS:")
            for warning in warnings:
                print(f"  {warning}")

        # Find UniProt IDs for the specified chain
        uniprot_ids = set()
        found_entries = []
        
        for row in table_data:
            if row['PDB'].lower() == pdb_code.lower() and row['Chain ID'] == chain_id:
                found_entries.append(row)
                if row['UniProt ID']:
                    uniprot_ids.add(row['UniProt ID'])
        
        uniprot_ids_list = sorted(list(uniprot_ids))
        
        if not uniprot_ids_list and verbose:
            print(f"No UniProt IDs found for PDB code {pdb_code} and chain {chain_id}.")
        
        if return_all_data:
            return uniprot_ids_list, table_data, warnings
        else:
            return uniprot_ids_list
            
    finally:
        # Clean up temporary directory
        try:
            shutil.rmtree(temp_dir)
        except:
            pass  # Ignore cleanup errors

def print_simplified_output(table_data):
    """Print simplified tabular output of the data."""
    # defaultdict with list to collect chain IDs
    data = defaultdict(list)
    
    # Group the chain IDs by UniProt accession and PDB ID
    for row in table_data:
        pdb_id = row['PDB']
        entity_id = row['Entity']
        chain_id = row['Chain ID']
        uniprot_ac = row['UniProt AC']
        uniprot_id = row['UniProt ID']
        hetatm = row['HETATM']
        resolution = row['Resolution'] 
        method = row['Method']

        # Group chain IDs by PDB ID and UniProt accession
        data[(pdb_id, entity_id, uniprot_ac, uniprot_id, hetatm, resolution, method)].append(chain_id)
    
    # Print the simplified output
    print("")
    print(f"{'PDB':<5} {'Entity':<6} {'Chain ID':<20} {'UniProt AC':<11} {'UniProt ID':<15} {'HETATM':<30} {'Resolution':<12} {'Method':<8}")
    print("-" * 120)
    for (pdb_id, entity_id, uniprot_ac, uniprot_id, hetatm, resolution, method), chain_ids in data.items():
        # Combine the chain IDs into a single string
        chain_ids_str = '/'.join(sorted(chain_ids))
        print(f"{pdb_id:<5} {entity_id:<6} {chain_ids_str:<20} {uniprot_ac:<11} {uniprot_id:<15} {hetatm:<30} {resolution:<12} {method:<8}")

# Legacy function for backward compatibility
def main():
    """Legacy main function for command-line usage."""
    if len(sys.argv) != 3:
        print("Usage: python script.py <PDB_CODE> <CHAIN_ID>")
        sys.exit(1)
    
    pdb_code = sys.argv[1].strip()
    chain_id = sys.argv[2].strip()
    
    uniprot_ids = pdb_unpID(pdb_code, chain_id, verbose=True)
    
    if uniprot_ids:
        print(",".join(uniprot_ids))

if __name__ == "__main__":
    main()