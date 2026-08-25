"""AMP model comparison page."""
import pandas as pd
import subprocess
import streamlit as st
from core.sequence import clean_sequence
from predictors.ampeppy import AmPEPpyPredictor
from predictors.base import PredictorUnavailable

@st.cache_resource
def predictor_registry():
    return [AmPEPpyPredictor()]

st.title("🧪 AMP prediction")
st.write("Predict antimicrobial peptide activity with the bundled amPEPpy pretrained classifier.")
predictors = predictor_registry()
available_predictors = []
unavailable_details = []
for predictor in predictors:
    available, detail = predictor.availability()
    if available:
        available_predictors.append(predictor)
    else:
        unavailable_details.append(f"{predictor.name}: {detail}")

if available_predictors:
    st.success("Available models: " + ", ".join(predictor.name for predictor in available_predictors))
else:
    st.error("No AMP predictors are installed in this deployment.")
    with st.expander("Deployment diagnostics"):
        for detail in unavailable_details:
            st.write(detail)
    st.stop()

sequence = st.text_area("Amino-acid sequence", value=st.session_state.pop("prediction_sequence", ""), height=160)
if st.button("Predict", type="primary"):
    try:
        sequence = clean_sequence(sequence, standard_only=True)
    except ValueError as exc:
        st.error(str(exc)); st.stop()
    results, errors = [], []
    for predictor in available_predictors:
        try:
            results.append(predictor.predict(sequence))
        except (PredictorUnavailable, subprocess.SubprocessError) as exc:
            errors.append(f"{predictor.name}: {exc}")
    for error in errors:
        st.warning(error)
    if results:
        table = pd.DataFrame([{"Model": result.model, "Prediction": result.prediction, "Probability": result.probability, **result.extra} for result in results])
        st.dataframe(table, use_container_width=True, hide_index=True)
        amps = sum(result.prediction.upper() == "AMP" for result in results)
        st.subheader(f"Consensus: {amps}/{len(results)} AMP")
    elif not errors:
        st.warning("No predictors are configured.")
