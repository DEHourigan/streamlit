"""AMP model comparison page."""
import pandas as pd
import subprocess
import streamlit as st
from core.sequence import clean_sequence
from predictors.ampeppy import AmPEPpyPredictor
from predictors.base import PredictorUnavailable
from predictors.macrel import MacrelPredictor

@st.cache_resource
def predictor_registry():
    return [MacrelPredictor(), AmPEPpyPredictor()]

st.title("🧪 AMP prediction")
st.write("Compare the bundled Macrel and amPEPpy pretrained AMP classifiers. Predictor execution remains isolated behind a common interface.")
predictors = predictor_registry()
for predictor in predictors:
    available, detail = predictor.availability()
    (st.success if available else st.info)(f"{predictor.name}: {'available' if available else detail}")
sequence = st.text_area("Amino-acid sequence", value=st.session_state.pop("prediction_sequence", ""), height=160)
if st.button("Predict", type="primary"):
    try:
        sequence = clean_sequence(sequence, standard_only=True)
    except ValueError as exc:
        st.error(str(exc)); st.stop()
    results, errors = [], []
    for predictor in predictors:
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
