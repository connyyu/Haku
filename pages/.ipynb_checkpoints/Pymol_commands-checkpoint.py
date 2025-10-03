import streamlit as st
import os

st.set_page_config(
    page_title="Haku - Useful PyMOL Commands",
    page_icon="💮", initial_sidebar_state="expanded",
)

script_dir = os.path.dirname(os.path.abspath(__file__)) # location of pages directory
output_dir = os.path.join(script_dir, "chain_info")

st.sidebar.markdown("")
st.sidebar.markdown("")

st.sidebar.markdown("### Author")
st.sidebar.markdown(
    """
    <b>Conny WH Yu</b>
    <p style="font-size: 15px; line-height: 1.8;">
    🎓 Structural biology and bioinformatics | 🧑🏻‍🔬 10+ years at the bench | 🧑🏻‍💻 Data curation | 🎯 Translate protein structure into functional understanding
    🌐 <a href='https://github.com/connyyu' target='_blank' style="font-size: 14px; color: #666; text-decoration: none;">GitHub profile</a>
    </p>
    """, 
    unsafe_allow_html=True
)

st.markdown(
    """
    ##### -- Useful PyMOL Commands --  
"""
)    
col1, col2 = st.columns([1, 2])
with col1:
    st.write("")
    st.write("_Remove solvent_")
    st.write("")
    st.write("")
    st.write("_Highlight disulfides_")
    st.write("")
    st.write("")
    st.write("_Highlight ligand (e.g. Zn)_")
    st.write("")
    st.write("")
    st.write("_Label residue (e.g. 647)_")
    st.write("")
    st.write("_Colour pLDDT score using AlphaFold2 colour scheme_")
    st.write("")
    st.write("")
    st.write("_Generate pretty figures_")
    st.write("")
with col2:
    st.code("remove solvent")
    st.code("select disulfides, br. CYS/SG and bound_to CYS/SG")
    st.code("sele zn, hetatm and resn zn")
    st.code("label n. CA and i. 647, '%s %s' % (resn, resi)")
    st.write("")
    st.code("spectrum b, rainbow_rev")
    pretty_code_file_path = os.path.join(output_dir, 'pretty_code.txt')
    with open(pretty_code_file_path, 'r') as file:
        pretty_code = file.read()
        st.code(pretty_code)

st.markdown(
    """
    ##### -- PyMOL Resources --  
"""
)
st.markdown(
    """    
    - PyMOL by Schrödinger: [PyMOL](https://www.pymol.org/)
    - Opensource PyMOL: [Github](https://github.com/schrodinger/pymol-open-source)
    - PyMOL wiki: [wiki](https://pymolwiki.org/index.php/Main_Page)
    
"""
)
st.markdown("""
    <a href='https://github.com/connyyu' target='_blank'>
        <img src='https://upload.wikimedia.org/wikipedia/commons/thumb/9/91/Octicons-mark-github.svg/1024px-Octicons-mark-github.svg.png' 
        style='position: fixed; bottom: 5%; left: 10%; transform: translateX(-50%); width: 30px; height: 30px;'/>
    </a>
""", unsafe_allow_html=True)