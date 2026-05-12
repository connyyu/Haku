import streamlit as st
import requests
from Bio.Align import PairwiseAligner, substitution_matrices
import streamlit.components.v1 as components

st.set_page_config(page_title="Haku - Sequence alignment", page_icon="💮")
default_unp = "G0SGT6"
default_pdb = "9i1w"
default_chain = "La"
st.session_state.guide_align = True

with st.sidebar:
    st.header("📝 New query")
    pdb_id_in = st.text_input("PDB code", value=default_pdb).strip()
    chain_in = st.text_input("PDB Chain ID", value=default_chain).strip()
    uniprot_ac = st.text_input("UniProt AC", value=default_unp).strip()
    fetch_button = st.button("Align sequences")
    if fetch_button:
        st.session_state.guide_align = False
    guide_align = st.checkbox("Instructions", value=st.session_state.get("guide_align", True))
    st.session_state.guide_align = guide_align

st.markdown("<a name='top_title'></a>", unsafe_allow_html=True)
st.markdown("#### Align sequences from PDB structure to UniProt protein sequence.")

if st.session_state.get('guide_align'):
    st.info("""
    Input:
    * **PDB code**
    * **PDB chain ID**
    * **UniProt AC**\n\n
    Output:\n
    The protein sequence from the specified chain in the PDB structure is aligned to the canonical sequence in the UniProt entry.
            The alignment is performed using local pairwise algorithm BLOSUM62, with a gap open penalty of -12 and a gap extension penalty of -4.
    """)

# ── Haku-Style White Theme CSS ────────────────────────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;600&family=Syne:wght@400;600;800&display=swap');

html, body, [class*="css"], .stApp { 
    background-color: #ffffff !important; 
    color: #1f2937 !important;
}

h1, h2, h3 { font-family: 'Syne', sans-serif; font-weight: 800; color: #111827 !important; }

.stTextInput > div > div > input {
    background: #f9fafb !important;
    border: 1px solid #e5e7eb !important;
    border-radius: 4px !important;
    color: #111827 !important;
}

div[data-testid="stHorizontalBlock"] .stButton:first-child > button {
    background-color: #3b82f6 !important;
    border: 1px solid #2563eb !important;
    color: #ffffff !important;
    border-radius: 4px !important;
    font-family: 'Syne', sans-serif !important;
    font-weight: 600 !important;
}

.protein-card {
    background: #ffffff;
    border: 1px solid #e5e7eb;
    border-radius: 4px;
    padding: 18px;
    margin: 8px 0;
}
.protein-card .ptype {
    font-size: 0.65rem;
    letter-spacing: 0.1em;
    text-transform: uppercase;
    color: #9ca3af;
    margin-bottom: 4px;
}
.protein-card .pname { font-size: 1.05rem; font-weight: 700; color: #111827; }
.protein-card .psub { font-size: 0.85rem; color: #6b7280; font-family: 'JetBrains Mono', monospace; }

.stat-chip {
    background: #f3f4f6;
    border-radius: 4px;
    padding: 4px 12px;
    font-size: 0.8rem;
    font-family: 'JetBrains Mono', monospace;
    color: #4b5563;
    display: inline-block;
    margin-right: 8px;
    margin-bottom: 8px;
}
.stat-chip span { color: #2563eb; font-weight: 700; }

.aln-block {
    background: #f3f4f6 !important;
    border: 1px solid #e5e7eb !important;
    border-radius: 4px;
    padding: 24px;
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.85rem;
    line-height: 1.7;
    overflow-x: auto;
    color: #374151 !important;
    margin-top: 15px;
}
.aln-line { white-space: pre; display: block; }
.aln-name { color: #1d4ed8; font-weight: 600; } 
.aln-pos  { color: #9ca3af; } 
.aln-seq  { color: #111827; } 
.aln-match { color: #059669; font-weight: 700; white-space: pre; }

.section-divider { border-top: 1px solid #f3f4f6; margin: 24px 0; }

a { color: #3b82f6; text-decoration: none; }
a:hover { text-decoration: underline; }
</style>
""", unsafe_allow_html=True)

# ── Helper functions ──────────────────────────────────────────────────────────

def fetch_pdbe_data(pdb_id, chain_id):
    """Fetches sequence and entity_id for a specific chain."""
    pdb_id = pdb_id.lower()
    # Get sequence
    fasta_url = f"https://www.ebi.ac.uk/pdbe/api/v2/pdb/entry/{pdb_id}/fasta"
    fasta_res = requests.get(fasta_url, timeout=10)
    fasta_res.raise_for_status()
    
    current_chain, sequences = None, {}
    for line in fasta_res.text.strip().splitlines():
        if line.startswith(">"):
            parts = line.split("|")
            if len(parts) >= 3:
                current_chain = parts[2].strip()
                sequences[current_chain] = ""
        elif current_chain:
            sequences[current_chain] += line.strip()
    
    seq = sequences.get(chain_id)
    if not seq:
        raise ValueError(f"Chain {chain_id} not found in {pdb_id}")

    # Get entity ID mapping
    entity_url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdb_id}"
    entity_res = requests.get(entity_url, timeout=10).json()
    entity_id = None
    for molecule in entity_res.get(pdb_id, []):
        if chain_id in molecule.get('in_chains', []):
            entity_id = molecule.get('entity_id')
            break
            
    return seq, entity_id

def fetch_uniprot_sequence(accession):
    response = requests.get(f"https://rest.uniprot.org/uniprotkb/{accession}.fasta", timeout=10)
    response.raise_for_status()
    lines = response.text.strip().split("\n")
    return "".join(line.strip() for line in lines if not line.startswith(">"))

def fetch_pdbe_title(pdb_id):
    try:
        url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{pdb_id.lower()}"
        data = requests.get(url, timeout=8).json()
        return data[pdb_id.lower()][0].get("title", "").capitalize()
    except: return None

def fetch_uniprot_metadata(accession):
    """Fetches name and Entry Name (ID)."""
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
        data = requests.get(url, timeout=8).json()
        name = data.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
        entry_name = data.get("uniProtkbId", "") 
        return name, entry_name
    except: return None, None

def run_alignment(seq1, seq2):
    aligner = PairwiseAligner()
    aligner.mode, aligner.open_gap_score, aligner.extend_gap_score = "local", -12, -4
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    alignment = aligner.align(seq1, seq2)[0]
    blocks1, blocks2 = alignment.aligned[0], alignment.aligned[1]
    aln_seq1, aln_seq2 = "", ""
    for i, ((s1, e1), (s2, e2)) in enumerate(zip(blocks1, blocks2)):
        if i > 0:
            aln_seq1 += "-" * (s1 - blocks1[i-1][1])
            aln_seq2 += "-" * (s2 - blocks2[i-1][1])
        aln_seq1 += seq1[s1:e1]
        aln_seq2 += seq2[s2:e2]
    
    matches = sum(1 for a, b in zip(aln_seq1, aln_seq2) if a == b and a != "-")
    gaps = sum(1 for a, b in zip(aln_seq1, aln_seq2) if a == "-" or b == "-")
    return {
        "aln_seq1": aln_seq1, "aln_seq2": aln_seq2,
        "s1_start": int(blocks1[0][0]) + 1, "s1_end": int(blocks1[-1][1]),
        "s2_start": int(blocks2[0][0]) + 1, "s2_end": int(blocks2[-1][1]),
        "length": len(aln_seq1), "score": alignment.score,
        "matches": matches, "gap_freq": (gaps / len(aln_seq1) * 100) if len(aln_seq1) > 0 else 0,
        "identity": (matches / len(aln_seq1) * 100) if len(aln_seq1) > 0 else 0,
        "coverage": (matches / max(len(seq1), len(seq2)) * 100)
    }

def copy_button_js(label, text, key):
    safe = text.replace("`", "\\`").replace("\\", "\\\\")
    components.html(f"""
    <button id="btn_{key}" style="background:#ffffff; border:1px solid #d1d5db; border-radius:4px; color:#374151; font-family:'Syne',sans-serif; padding:6px 12px; cursor:pointer; font-size:12px; font-weight:600;">{label}</button>
    <script>
    document.getElementById('btn_{key}').onclick = function() {{
        navigator.clipboard.writeText(`{safe}`).then(() => {{
            this.textContent = 'Copied';
            this.style.color = '#10b981';
            setTimeout(() => {{ this.textContent = '{label}'; this.style.color = '#374151'; }}, 2000);
        }});
    }};
    </script>
    """, height=40)

# ── Main UI ───────────────────────────────────────────────────────────────────

if fetch_button:
    try:
        # Fetch Data
        s1, entity_id = fetch_pdbe_data(pdb_id_in, chain_in)
        s2 = fetch_uniprot_sequence(uniprot_ac)
        p_title = fetch_pdbe_title(pdb_id_in)
        u_name, u_id = fetch_uniprot_metadata(uniprot_ac)
        res = run_alignment(s1, s2)

        # Protein Info Cards
        i1, i2 = st.columns(2)
        with i1:
            pdb_link = f"https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id_in.lower()}?activeTab=macromolecules&entity={entity_id}"
            st.markdown(f"""
                <div class="protein-card">
                    <div class="ptype">Structure (<a href="{pdb_link}" target="_blank">PDB {pdb_id_in.upper()}</a>)</div>
                    <div class="pname">{p_title or "N/A"}</div>
                    <div class="psub">Chain {chain_in} · {len(s1)} AA</div>
                </div>""", unsafe_allow_html=True)
            
        with i2:
            uniprot_link = f"https://www.uniprot.org/uniprotkb/{uniprot_ac}/entry"
            st.markdown(f"""
                <div class="protein-card">
                    <div class="ptype">Protein (<a href="{uniprot_link}" target="_blank">UniProt {uniprot_ac}</a>)</div>
                    <div class="pname">{u_name}</div>
                    <div class="psub">{u_id} · {len(s2)} AA</div>
                </div>""", unsafe_allow_html=True)
            
        c1, c2 = st.columns(2)
        with c1:
            copy_button_js("Copy PDB sequence", s1, "s1")
        with c2:
            copy_button_js("Copy UniProt sequence", s2, "s2")

        # Alignment Stats and Block
        name1, name2 = f"{pdb_id_in.upper()}.{chain_in}", uniprot_ac
        st.markdown(f"""
        <div style="font-size:1.1rem; color:#4b5563; margin-bottom:12px; margin-top:20px;">
            <b>{name1}</b>: residues {res['s1_start']}–{res['s1_end']} aligned to <b>{name2}</b>: residues {res['s2_start']}–{res['s2_end']}
        </div>
        <div style="margin-bottom: 20px;">
            <div class="stat-chip">Identity <span>{res['identity']:.1f}%</span></div>
            <div class="stat-chip">Coverage <span>{res['coverage']:.1f}%</span></div>
            <div class="stat-chip">Overlap <span>{res['length']} AA</span></div>
            <div class="stat-chip">Gap Freq <span>{res['gap_freq']:.1f}%</span></div>
            <div class="stat-chip">Score <span>{res['score']:.1f}</span></div>
        </div>
        """, unsafe_allow_html=True)

        aln_out = []
        p1, p2 = res['s1_start'], res['s2_start']
        for i in range(0, res['length'], 60):
            b1, b2 = res['aln_seq1'][i:i+60], res['aln_seq2'][i:i+60]
            m = "".join("*" if a == b and a != "-" else " " for a, b in zip(b1, b2))
            aln_out.append(f'<div class="aln-line"><span class="aln-name">{name1:<12}</span> <span class="aln-pos">{p1:>4}</span> <span class="aln-seq">{b1}</span></div>')
            aln_out.append(f'<div class="aln-line"><span class="aln-name">{name2:<12}</span> <span class="aln-pos">{p2:>4}</span> <span class="aln-seq">{b2}</span></div>')
            aln_out.append(f'<div class="aln-line">{"":18}<span class="aln-match">{m}</span></div>')
            aln_out.append('<div> </div>')
            p1 += len(b1.replace("-",""))
            p2 += len(b2.replace("-",""))

        st.markdown(f'<div class="aln-block">{"".join(aln_out)}</div>', unsafe_allow_html=True)

    except Exception as e:
        st.error(f"Error: {str(e)}")