"""Search and inspect FASTA database records."""
import pandas as pd
import streamlit as st
from core.databases import load_databases
from core.fasta import read_fasta

@st.cache_data(show_spinner=False)
def cached_records(path: str, modified: float):
    del modified
    return read_fasta(path)

st.title("🗂️ Database browser")
try:
    databases = {key: db for key, db in load_databases().items() if db.fasta_available}
except Exception as exc:
    st.error(str(exc)); st.stop()
if not databases:
    st.error("No configured FASTA database is available."); st.stop()
key = st.selectbox("Database", databases, format_func=lambda item: databases[item].name)
database = databases[key]
try:
    records = cached_records(str(database.fasta), database.fasta.stat().st_mtime)
except Exception as exc:
    st.error(f"Could not read FASTA: {exc}"); st.stop()

maximum = max((len(record.sequence) for record in records), default=1)
query = st.text_input("Search IDs or descriptions", placeholder="Enter a name or identifier")
length_range = st.slider("Sequence length (aa)", 1, maximum, (1, maximum))
needle = query.casefold().strip()
filtered = [record for record in records if length_range[0] <= len(record.sequence) <= length_range[1] and (not needle or needle in record.id.casefold() or needle in record.description.casefold())]
st.caption(f"Showing {len(filtered):,} of {len(records):,} sequences")
table = pd.DataFrame([{"ID": record.id, "Description": record.description, "Length": len(record.sequence), "Sequence": record.sequence} for record in filtered])
st.dataframe(table, use_container_width=True, hide_index=True, column_config={"Sequence": st.column_config.TextColumn(width="large")})
if filtered:
    selected = st.selectbox("Inspect sequence", range(len(filtered)), format_func=lambda index: f"{filtered[index].id} — {filtered[index].description}")
    record = filtered[selected]
    st.subheader(record.id)
    st.write(record.description)
    st.metric("Length", f"{len(record.sequence)} aa")
    st.code(record.sequence, language=None, wrap_lines=True)
    left, right = st.columns(2)
    if left.button("Send to BLAST"):
        st.session_state["blast_sequence"] = record.sequence
        st.switch_page("pages/blast_search.py")
    if right.button("Send to sequence analysis"):
        st.session_state["analysis_sequence"] = record.sequence
        st.switch_page("pages/sequence_analysis.py")
