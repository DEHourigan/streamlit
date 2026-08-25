"""Entry point for the peptide analysis platform."""
from pathlib import Path
import sys
import streamlit as st

SRC_DIR = Path(__file__).resolve().parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

st.set_page_config(page_title="Peptide analysis", page_icon="🧬", layout="wide")
pages = {
    "Overview": [st.Page("pages/home.py", title="Home", icon="🏠", default=True)],
    "Analysis": [
        st.Page("pages/blast_search.py", title="BLAST Search", icon="🔎"),
        st.Page("pages/database_browser.py", title="Database Browser", icon="🗂️"),
        st.Page("pages/amp_prediction.py", title="AMP Prediction", icon="🧪"),
        st.Page("pages/sequence_analysis.py", title="Sequence Analysis", icon="📏"),
    ],
    "Project": [st.Page("pages/about.py", title="About", icon="ℹ️")],
}
st.navigation(pages).run()
