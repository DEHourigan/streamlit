"""BLAST search UI; execution lives in core.blast."""
import streamlit as st
from core.alignment import format_alignment
from core.blast import BlastError, run_blast
from core.databases import load_databases
from core.sequence import clean_sequence

DISPLAY_COLUMNS = ["Hit", "Identity (%)", "Query coverage (%)", "Alignment length", "Subject length", "E-value", "Bit score", "Description"]

st.title("🔎 Peptide BLAST search")
try:
    databases = load_databases()
except Exception as exc:
    st.error(str(exc)); st.stop()
available = {key: db for key, db in databases.items() if db.blast_available}
if not available:
    st.error("No enabled database has BLAST index files. Check config/databases.toml and build the configured indexes."); st.stop()

database_key = st.selectbox("Database", available, format_func=lambda key: available[key].name)
database = available[database_key]
sequence = st.text_area("Enter amino-acid sequence", value=st.session_state.pop("blast_sequence", ""), height=160, placeholder="MSTK...")
left, right = st.columns(2)
evalue = left.number_input("E-value threshold", min_value=0.0, value=10.0, format="%.2e")
max_hits = right.number_input("Maximum number of hits", min_value=1, max_value=1000, value=50)

if st.button("Search", type="primary"):
    try:
        sequence = clean_sequence(sequence)
        with st.spinner("Running BLAST…"):
            result = run_blast(sequence, database.blast_db, evalue, int(max_hits))
        st.session_state["blast_result"] = result
        st.session_state["blast_database"] = database_key
        st.session_state["blast_query"] = sequence
    except (ValueError, BlastError) as exc:
        st.error(str(exc))

result = st.session_state.get("blast_result")
if result is not None and st.session_state.get("blast_database") == database_key:
    st.info(f"Query length: **{len(st.session_state['blast_query'])} aa** | BLAST task: **{result.task}**")
    hits = result.hits
    if hits.empty:
        st.warning("No hits found.")
    else:
        st.subheader(f"BLAST results ({len(hits)} alignments)")
        st.dataframe(hits[DISPLAY_COLUMNS], use_container_width=True, hide_index=True)
        st.download_button("Download results as CSV", hits.to_csv(index=False).encode(), "blast_results.csv", "text/csv")
        options = [f"{i + 1}. {row['Hit']} — {row['Identity (%)']:.1f}% identity, {row['Query coverage (%)']:.1f}% coverage" for i, row in hits.iterrows()]
        selected = st.selectbox("Select a hit", range(len(hits)), format_func=lambda index: options[index])
        hit = hits.iloc[selected]
        metrics = st.columns(4)
        metrics[0].metric("Identity", f"{hit['Identity (%)']:.1f}%")
        metrics[1].metric("Query coverage", f"{hit['Query coverage (%)']:.1f}%")
        metrics[2].metric("E-value", f"{hit['E-value']:.2e}")
        metrics[3].metric("Bit score", f"{hit['Bit score']:.1f}")
        st.write(f"**Hit:** `{hit['Hit']}`  \n**Description:** {hit['Description']}")
        st.code(format_alignment(hit["Aligned query"], hit["Aligned subject"], hit["Query start"], hit["Subject start"]), language=None)
        st.caption("│ = identical residue     · = substitution")
