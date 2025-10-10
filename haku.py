import streamlit as st

st.set_page_config(
    page_title="Haku - protein structure navigator",
    page_icon="💮",
)

pages = [
        st.Page("pages/Find_structures.py", title="Find the structures"),
        st.Page("pages/List_UniProt_entries.py", title="List the proteins"),
        st.Page("pages/View_interactions.py", title="View the interactions"),
        st.Page("pages/Structure_viewer.py", title="3D viewer"),
        st.Page("pages/TM_annotations.py", title="TM annotations"),
        st.Page("pages/Pymol_commands.py", title="PyMOL commands"),
        st.Page("pages/About_haku.py", title="About Haku"),
]

pg = st.navigation(pages, position="top", expanded=True)
pg.run()




