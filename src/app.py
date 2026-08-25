import os
import re
import subprocess
import tempfile
from io import StringIO

import pandas as pd
import streamlit as st


# ---------------------------------------------------------
# Configuration
# ---------------------------------------------------------

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

BLAST_DB = os.path.join(
    BASE_DIR,
    "data",
    "blastdb",
    "bagel_db"
)

FASTA_DB = os.path.join(
    BASE_DIR,
    "data",
    "blastdb",
    "bagel_db.faa"
)


# ---------------------------------------------------------
# Database statistics
# ---------------------------------------------------------

@st.cache_data
def get_database_stats(fasta_file):
    """
    Parse FASTA database and calculate sequence length statistics.
    """

    lengths = []
    sequence = []

    with open(fasta_file) as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if sequence:
                    seq = "".join(sequence).replace("*", "")
                    lengths.append(len(seq))
                    sequence = []
            else:
                sequence.append(line)

        # final sequence
        if sequence:
            seq = "".join(sequence).replace("*", "")
            lengths.append(len(seq))

    lengths = pd.Series(lengths)

    stats = {
        "n_sequences": len(lengths),
        "total_residues": int(lengths.sum()),
        "mean_length": lengths.mean(),
        "sd_length": lengths.std(ddof=0),
        "median_length": lengths.median(),
        "min_length": int(lengths.min()),
        "max_length": int(lengths.max()),
    }

    return stats, lengths


# ---------------------------------------------------------
# Alignment formatter
# ---------------------------------------------------------

def format_alignment(qseq, sseq, qstart, sstart, width=60):
    """
    Format BLAST alignment into readable blocks.

    Exact match = │
    mismatch    = ·
    gap         = space
    """

    blocks = []

    q_position = int(qstart)
    s_position = int(sstart)

    for i in range(0, len(qseq), width):

        qblock = qseq[i:i + width]
        sblock = sseq[i:i + width]

        match_line = ""

        for q, s in zip(qblock, sblock):

            if q == "-" or s == "-":
                match_line += " "

            elif q == s:
                match_line += "│"

            else:
                match_line += "·"

        # Count actual residues, ignoring gaps
        q_residues = len(qblock.replace("-", ""))
        s_residues = len(sblock.replace("-", ""))

        q_end = q_position + q_residues - 1
        s_end = s_position + s_residues - 1

        block = (
            f"Query  {q_position:<6} {qblock}  {q_end}\n"
            f"               {match_line}\n"
            f"Hit    {s_position:<6} {sblock}  {s_end}"
        )

        blocks.append(block)

        q_position = q_end + 1
        s_position = s_end + 1

    return "\n\n".join(blocks)


# ---------------------------------------------------------
# Streamlit page
# ---------------------------------------------------------

st.set_page_config(
    page_title="Peptide BLAST",
    page_icon="🧬",
    layout="wide"
)

st.title("Peptide BLAST")

st.write(
    "Search a peptide or protein sequence against the BAGEL peptide database."
)


# ---------------------------------------------------------
# Database overview
# ---------------------------------------------------------

stats, lengths = get_database_stats(FASTA_DB)

st.subheader("Database overview")

col1, col2, col3, col4 = st.columns(4)

with col1:
    st.metric(
        "Peptides",
        f"{stats['n_sequences']:,}"
    )

with col2:
    st.metric(
        "Mean length",
        f"{stats['mean_length']:.1f} ± {stats['sd_length']:.1f} aa"
    )

with col3:
    st.metric(
        "Median length",
        f"{stats['median_length']:.0f} aa"
    )

with col4:
    st.metric(
        "Length range",
        f"{stats['min_length']}–{stats['max_length']} aa"
    )


with st.expander("More database information"):

    col1, col2 = st.columns(2)

    with col1:
        st.write(
            f"""
            **Database:** BAGEL peptide database  
            **Number of peptides:** {stats['n_sequences']:,}  
            **Total amino acids:** {stats['total_residues']:,}  
            """
        )

    with col2:
        st.write(
            f"""
            **Mean length:** {stats['mean_length']:.2f} aa  
            **Standard deviation:** {stats['sd_length']:.2f} aa  
            **Median length:** {stats['median_length']:.0f} aa  
            """
        )

    st.caption("Peptide length distribution")

    length_counts = (
        lengths
        .value_counts()
        .sort_index()
        .rename_axis("Length")
        .rename("Number of peptides")
    )

    st.bar_chart(length_counts)


st.divider()


# ---------------------------------------------------------
# Query input
# ---------------------------------------------------------

st.subheader("Search database")

sequence = st.text_area(
    "Enter amino-acid sequence",
    height=160,
    placeholder="MSTK..."
)

col1, col2 = st.columns(2)

with col1:
    evalue = st.number_input(
        "E-value threshold",
        min_value=0.0,
        value=10.0,
        format="%.2e"
    )

with col2:
    max_hits = st.number_input(
        "Maximum number of hits",
        min_value=1,
        max_value=1000,
        value=50
    )


# ---------------------------------------------------------
# Search
# ---------------------------------------------------------

if st.button("Search", type="primary"):

    # Clean sequence
    sequence = re.sub(r"\s+", "", sequence).upper()

    if sequence.startswith(">"):
        st.error(
            "Please paste the sequence only, without a FASTA header."
        )
        st.stop()

    if not sequence:
        st.error(
            "Please enter an amino-acid sequence."
        )
        st.stop()

    valid_aa = set(
        "ABCDEFGHIKLMNPQRSTVWXYZJUO*"
    )

    invalid = set(sequence) - valid_aa

    if invalid:
        st.error(
            "Invalid characters in sequence: "
            + ", ".join(sorted(invalid))
        )
        st.stop()

    query_length = len(sequence)

    # Short peptide mode
    if query_length < 30:
        task = "blastp-short"
    else:
        task = "blastp"

    st.info(
        f"Query length: **{query_length} aa**  |  "
        f"BLAST task: **{task}**"
    )

    # Write temporary FASTA
    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".fasta",
        delete=False
    ) as tmp:

        tmp.write(">query\n")
        tmp.write(sequence + "\n")

        query_file = tmp.name


    # -----------------------------------------------------
    # BLAST output fields
    # -----------------------------------------------------

    outfmt = (
        "6 "
        "qseqid "
        "sseqid "
        "pident "
        "length "
        "mismatch "
        "gapopen "
        "qstart "
        "qend "
        "sstart "
        "send "
        "evalue "
        "bitscore "
        "qcovs "
        "qlen "
        "slen "
        "qseq "
        "sseq "
        "stitle"
    )

    command = [
        "blastp",
        "-task",
        task,
        "-query",
        query_file,
        "-db",
        BLAST_DB,
        "-evalue",
        str(evalue),
        "-max_target_seqs",
        str(max_hits),
        "-outfmt",
        outfmt
    ]

    result = subprocess.run(
        command,
        capture_output=True,
        text=True
    )

    os.unlink(query_file)

    if result.returncode != 0:
        st.error("BLAST failed.")
        st.code(result.stderr)
        st.stop()

    if not result.stdout.strip():
        st.warning("No hits found.")
        st.stop()


    # -----------------------------------------------------
    # Parse BLAST
    # -----------------------------------------------------

    columns = [
        "Query",
        "Hit",
        "Identity (%)",
        "Alignment length",
        "Mismatches",
        "Gap opens",
        "Query start",
        "Query end",
        "Subject start",
        "Subject end",
        "E-value",
        "Bit score",
        "Query coverage (%)",
        "Query length",
        "Subject length",
        "Aligned query",
        "Aligned subject",
        "Description"
    ]

    df = pd.read_csv(
        StringIO(result.stdout),
        sep="\t",
        names=columns
    )


    # -----------------------------------------------------
    # Results table
    # -----------------------------------------------------

    st.divider()

    st.subheader(
        f"BLAST results ({len(df)} alignments)"
    )

    display_columns = [
        "Hit",
        "Identity (%)",
        "Query coverage (%)",
        "Alignment length",
        "Subject length",
        "E-value",
        "Bit score",
        "Description"
    ]

    st.dataframe(
        df[display_columns],
        use_container_width=True,
        hide_index=True
    )


    # -----------------------------------------------------
    # Alignment viewer
    # -----------------------------------------------------

    st.subheader("Alignment viewer")

    alignment_options = []

    for i, row in df.iterrows():

        alignment_options.append(
            f"{i + 1}. {row['Hit']} "
            f"— {row['Identity (%)']:.1f}% identity, "
            f"{row['Query coverage (%)']:.1f}% coverage"
        )

    selected = st.selectbox(
        "Select a hit",
        options=range(len(df)),
        format_func=lambda x: alignment_options[x]
    )

    hit = df.iloc[selected]


    # -----------------------------------------------------
    # Alignment statistics
    # -----------------------------------------------------

    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric(
            "Identity",
            f"{hit['Identity (%)']:.1f}%"
        )

    with col2:
        st.metric(
            "Query coverage",
            f"{hit['Query coverage (%)']:.1f}%"
        )

    with col3:
        st.metric(
            "E-value",
            f"{hit['E-value']:.2e}"
        )

    with col4:
        st.metric(
            "Bit score",
            f"{hit['Bit score']:.1f}"
        )


    st.write(
        f"**Hit:** `{hit['Hit']}`"
    )

    st.write(
        f"**Description:** {hit['Description']}"
    )


    # -----------------------------------------------------
    # Pretty alignment
    # -----------------------------------------------------

    alignment = format_alignment(
        hit["Aligned query"],
        hit["Aligned subject"],
        hit["Query start"],
        hit["Subject start"]
    )

    st.code(
        alignment,
        language=None
    )

    st.caption(
        "│ = identical residue     · = substitution"
    )


    # -----------------------------------------------------
    # Download results
    # -----------------------------------------------------

    st.divider()

    csv = df.to_csv(index=False)

    st.download_button(
        "Download complete BLAST results",
        csv,
        file_name="blast_results.csv",
        mime="text/csv"
    )