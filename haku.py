import streamlit as st

st.set_page_config(
    page_title="Haku - protein structure navigator",
    page_icon="💮",
)

pages = [
    st.Page("pages/Find_structures.py",      title="Find the structures"),
    st.Page("pages/List_UniProt_entries.py", title="List the proteins"),
    st.Page("pages/View_interactions.py",    title="View the interactions"),
    st.Page("pages/Structure_viewer.py",     title="3D viewer"),
    st.Page("pages/TM_annotations.py",       title="TM annotations"),
    st.Page("pages/Homologs.py",             title="Similar structures"),
    st.Page("pages/Pymol_commands.py",       title="PyMOL commands"),
    st.Page("pages/About_haku.py",           title="About Haku"),
]

pg = st.navigation(pages, position="top", expanded=True)
pg.run()

# Fix #11: use simple fixed positioning without the redundant transform
st.markdown("""
    <a href='https://github.com/connyyu' target='_blank'>
        <img src='https://upload.wikimedia.org/wikipedia/commons/thumb/9/91/Octicons-mark-github.svg/1024px-Octicons-mark-github.svg.png'
        style='position: fixed; bottom: 20px; left: 20px; width: 30px; height: 30px;'/>
    </a>
""", unsafe_allow_html=True)
