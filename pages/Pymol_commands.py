import streamlit as st
import os

st.set_page_config(
    page_title="Haku - Useful PyMOL Commands",
    page_icon="💮",
    initial_sidebar_state="expanded",
)

script_dir = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.join(script_dir, "chain_info")

st.sidebar.markdown("")
st.sidebar.markdown("")
st.sidebar.markdown("### Author")
st.sidebar.markdown(
    """
    <b>Conny WH Yu</b>
    <p style="font-size: 15px; line-height: 1.8;">
    🎓 Data Science and Structural Biology | 🧑🏻‍🔬 10+ years at the bench | 🧑🏻‍💻 Data curation | 🎯 Translate protein structure into functional understanding
    🌐 <a href='https://github.com/connyyu' target='_blank' style="font-size: 14px; color: #666; text-decoration: none;">GitHub profile</a>
    </p>
    """,
    unsafe_allow_html=True
)

st.markdown("##### -- Useful PyMOL Commands --")

st.markdown("**For selection**")

col1, col2 = st.columns([1, 2])
with col1:
    st.write("")
    st.write("_Highlight disulfides_")
with col2:
    st.code("select disulfides, br. CYS/SG and bound_to CYS/SG")

col1, col2 = st.columns([1, 2])
with col1:
    st.write("")
    st.write("_Highlight ligand (e.g. Zn)_")
with col2:
    st.code("sele ZN, hetatm and resn zn")

col1, col2 = st.columns([1, 2])
with col1:
    st.write("")
    st.write("_List the selected residues_")
with col2:
    st.code("iterate sele and name CA, print(f'{chain}:{resi}:{resn}')")

col1, col2 = st.columns([1, 2])
with col1:
    st.write("_Renumbering residues (e.g. obj01)_")
with col2:
    st.code("alter obj01 and chain A, resi=str(int(resi)+12)")

st.markdown("**For visualisation**")

col1, col2 = st.columns([1, 2])
with col1:
    st.write("")
    st.write("_Remove solvent_")
with col2:
    st.code("remove solvent")

col1, col2 = st.columns([1, 2])
with col1:
    st.write("")
    st.write("_Label residue (e.g. 90)_")
with col2:
    st.code("label n. CA and i. 90, '%s %s' % (resn, resi)")

col1, col2 = st.columns([1, 2])
with col1:
    st.write("_Colour pLDDT score using AlphaFold2 colour scheme_")
with col2:
    st.code("spectrum b, rainbow_rev")

st.markdown("**For publication**")

col1, col2 = st.columns([1, 2])
with col1:
    st.write("")
    st.write("_Generate pretty figures_")
    st.write("")
with col2:
    pretty_code_file_path = os.path.join(output_dir, 'pretty_code.txt')
    with open(pretty_code_file_path, 'r') as f:
        st.code(f.read())

st.markdown("##### -- PyMOL Resources --")
st.markdown("""
- PyMOL by Schrödinger: [PyMOL](https://www.pymol.org/)
- Opensource PyMOL: [Github](https://github.com/schrodinger/pymol-open-source)
- PyMOL wiki: [wiki](https://pymolwiki.org/index.php/Main_Page)
""")

# Fix #11: cleaner GitHub icon positioning
st.markdown("""
    <a href='https://github.com/connyyu' target='_blank'>
        <img src='https://upload.wikimedia.org/wikipedia/commons/thumb/9/91/Octicons-mark-github.svg/1024px-Octicons-mark-github.svg.png'
        style='position: fixed; bottom: 20px; left: 20px; width: 30px; height: 30px;'/>
    </a>
""", unsafe_allow_html=True)
