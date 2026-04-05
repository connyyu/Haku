# Haku, the Protein Structure Navigator.

This repository contains a Streamlit app called **Haku**, a collection of tools to navigate protein structures.

[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://haku3d.streamlit.app/)

## Features

- Find structures — Find all the structures of a protein using its UniProt AC.
- List proteins — List all the UniProt protein entries in the structure(s).
- View interactions — View all the ligand-protein interactions in the structure(s).
- Transmembrane annotations — Visualise transmembrane annotation on a protein structure.
- Similar structures — Find similar proteins with experimental structures.
- 3D viewer — Navigate protein structures in a 3D viewer (py3dmol / Mol*).
  
### Input
- **UniProt accession**
- **PDB code**

## Prerequisites

- **Python 3.x**
- Python libraries: `streamlit`, `pandas`, `requests`, `pybiolib`, `matplotlib`, `aiohttp`, `tenacity`, `py3Dmol`
    
## Author

- **Conny Yu** – [GitHub Profile](https://github.com/connyyu)  
  Haku 2.2 _April 2026_
  Haku 2.0 _October 2025_
