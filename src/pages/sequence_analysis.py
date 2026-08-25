"""Physicochemical sequence analysis page."""
import pandas as pd
import streamlit as st
from core.sequence import calculate_properties

st.title("📏 Sequence analysis")
st.write("Calculate peptide properties using Biopython ProtParam. Net charge is a simple neutral-pH side-chain approximation.")
sequence = st.text_area("Amino-acid sequence", value=st.session_state.pop("analysis_sequence", ""), height=160)
if st.button("Analyse", type="primary"):
    try:
        properties = calculate_properties(sequence)
        st.session_state["sequence_properties"] = properties
    except (ValueError, RuntimeError) as exc:
        st.error(str(exc))
properties = st.session_state.get("sequence_properties")
if properties:
    first = st.columns(4)
    first[0].metric("Length", f"{properties['length']} aa")
    first[1].metric("Molecular weight", f"{properties['molecular_weight']:.1f} Da")
    first[2].metric("Predicted pI", f"{properties['isoelectric_point']:.2f}")
    first[3].metric("Approx. net charge", f"{properties['net_charge']:+.1f}")
    second = st.columns(4)
    second[0].metric("Hydrophobicity (GRAVY)", f"{properties['hydrophobicity']:.3f}")
    second[1].metric("Aromaticity", f"{properties['aromaticity']:.3f}")
    second[2].metric("Instability index", f"{properties['instability_index']:.2f}")
    second[3].metric("Cysteine", f"{properties['cysteine_count']} ({properties['cysteine_fraction']:.1%})")
    composition = pd.DataFrame({"Amino acid": properties["composition"].keys(), "Fraction": properties["composition"].values()}).set_index("Amino acid")
    st.subheader("Amino-acid composition")
    st.bar_chart(composition)
