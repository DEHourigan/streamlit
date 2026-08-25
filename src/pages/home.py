"""Platform overview and database summaries."""
import pandas as pd
import streamlit as st
from core.databases import load_databases
from core.fasta import database_stats, read_fasta

@st.cache_data(show_spinner=False)
def cached_summary(path: str, modified: float):
    del modified
    return database_stats(read_fasta(path))

st.title("🧬 Peptide analysis platform")
st.write("Search reference peptide databases, browse their contents, predict antimicrobial activity, and calculate sequence properties.")
try:
    databases = load_databases()
except Exception as exc:
    st.error(str(exc)); st.stop()
if not databases:
    st.warning("No enabled databases are configured in config/databases.toml."); st.stop()

for database in databases.values():
    st.subheader(database.name)
    st.write(database.description)
    if not database.fasta_available:
        st.error(f"FASTA file is missing: {database.fasta}")
        continue
    try:
        stats, lengths = cached_summary(str(database.fasta), database.fasta.stat().st_mtime)
    except Exception as exc:
        st.error(f"Could not read database: {exc}"); continue
    columns = st.columns(5)
    columns[0].metric("Peptides", f"{stats['n_sequences']:,}")
    columns[1].metric("Total residues", f"{stats['total_residues']:,}")
    columns[2].metric("Mean length", f"{stats['mean_length']:.1f} ± {stats['sd_length']:.1f} aa")
    columns[3].metric("Median", f"{stats['median_length']:.0f} aa")
    columns[4].metric("Range", f"{stats['min_length']}–{stats['max_length']} aa")
    counts = pd.Series(lengths).value_counts().sort_index().rename("Number of peptides")
    st.caption("Peptide length distribution")
    st.bar_chart(counts, x_label="Length (aa)", y_label="Number of peptides")

st.subheader("Available tools")
st.markdown("""
- **BLAST Search** finds similar database peptides and displays pairwise alignments.
- **Database Browser** searches and inspects reference sequences.
- **AMP Prediction** provides an extensible interface for optional pretrained classifiers.
- **Sequence Analysis** calculates physicochemical peptide properties.
""")
