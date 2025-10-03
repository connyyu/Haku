import streamlit as st

st.set_page_config(
    page_title="Haku - Protein structure navigator",
    page_icon="💮", initial_sidebar_state="expanded",
)

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

st.markdown("## 💮 Haku")

st.markdown(
    """
    ##### -- Protein Structure Navigator --  
"""
)    
st.markdown(
    """
    There are so many protein structures.<br>
    The structures are getting larger and more complex than ever.<br>
    <b>Haku</b> is a navigator in the new world of protein structures.<br>
    <i>Bon voyage.</i>""", unsafe_allow_html=True)

st.markdown(
    """    
    #### Features
    - **Find structures**
    -- Find all the structures of a protein using its UniProt AC.
    - **List proteins**
    -- List all the UniProt protein entries in the structure(s).
    - **View interactions**
    -- View all the ligand-protein interactions in the structure(s). 
    - **3D viewer**
    -- Navigate protein structures in a 3D viewer.
    - **Transmembrane annotations**
    -- Visualise transmembrane annotation on a protein structure.
    
    #### Resources
    - Protein database: [UniProt](https://https://www.uniprot.org/)
    - Structure database: [PDBe](https://www.ebi.ac.uk/pdbe/)
    - Interactions in protein structures: [Arpeggio](https://www.sciencedirect.com/science/article/pii/S0022283616305332)
    - Literature database: [Europe PMC](https://europepmc.org/)
    - Structure prediction: [AlphaFold2](https://alphafold.ebi.ac.uk/)
    - Transmembrane topology prediction: [DeepTMHMM](https://dtu.biolib.com/DeepTMHMM)
"""
)
st.markdown("""
    <a href='https://github.com/connyyu' target='_blank'>
        <img src='https://upload.wikimedia.org/wikipedia/commons/thumb/9/91/Octicons-mark-github.svg/1024px-Octicons-mark-github.svg.png' 
        style='position: fixed; bottom: 5%; left: 10%; transform: translateX(-50%); width: 30px; height: 30px;'/>
    </a>
""", unsafe_allow_html=True)