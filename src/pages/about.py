"""Versions, citations, and methodology."""
import importlib.metadata
import shutil
import subprocess
import streamlit as st
from core.databases import load_databases
from predictors.ampeppy import AmPEPpyPredictor

st.title("ℹ️ About")
st.write("This platform combines local peptide databases, NCBI BLAST+, Biopython ProtParam, and optional pretrained AMP classifiers.")
st.subheader("Databases")
try:
    for database in load_databases(enabled_only=False).values():
        st.markdown(f"**{database.name}** — {database.version}  \n{database.description}  \nCitation: {database.citation}")
except Exception as exc:
    st.error(str(exc))
st.subheader("Software and models")
blastp = shutil.which("blastp")
blast_version = subprocess.run([blastp, "-version"], capture_output=True, text=True).stdout.splitlines()[0] if blastp else "Not installed"
try:
    biopython_version = importlib.metadata.version("biopython")
except importlib.metadata.PackageNotFoundError:
    biopython_version = "Not installed"
st.write(f"NCBI BLAST+: {blast_version}  \nBiopython: {biopython_version}")
for predictor in [AmPEPpyPredictor()]:
    available, detail = predictor.availability()
    st.write(f"{predictor.name}: {detail if not available else 'available'}")
st.markdown("""
**Predictor citations**

- Lawrence et al. (2021), *amPEPpy 1.0: a portable and accurate antimicrobial peptide prediction tool*, Bioinformatics 37:2058–2060.
""")
st.subheader("Methodology")
st.markdown("""
- Queries shorter than 30 residues use `blastp-short`; longer queries use `blastp`.
- Database statistics are calculated directly from each registered FASTA file.
- Physicochemical properties use Biopython ProtParam; the displayed net charge is an approximate neutral-pH side-chain count.
- AMP classifications are produced only by installed model wrappers; the platform does not substitute heuristic predictions.
""")
