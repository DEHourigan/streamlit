"""Lightweight, self-contained amPEPpy inference adapter.

This computes the same 105 CTD distribution descriptors and runs the bundled
official random-forest model directly, avoiding a legacy Git package build during
Streamlit deployment.
"""
from pathlib import Path
import pickle
import sys

from .base import BasePredictor, PredictionResult, PredictorUnavailable

DEFAULT_MODEL = Path(__file__).resolve().parents[2] / "data" / "models" / "ampeppy" / "amPEP.model"
CTD_GROUPS = {
    "hydrophobicity": ("RKEDQN", "GASTPHY", "CLVIMFW"),
    "normalized.van.der.waals": ("GASTPDC", "NVEQIL", "MHKFRYW"),
    "polarity": ("LIFWCMVY", "PATGS", "HQRKNED"),
    "polarizability": ("GASDT", "CPNVEQIL", "KMHFRYW"),
    "charge": ("KR", "ANCQGHILMFPSTWYV", "DE"),
    "secondary": ("EALMQKRH", "VIYCWFT", "GNPDS"),
    "solvent": ("ALFCGIVW", "RKQEND", "MSPTHY"),
}
PERCENTILES = (0, 25, 50, 75, 100)

def _distribution(positions: list[int]) -> list[float]:
    """Reproduce amPEPpy's 0/25/50/75/100 distribution positions."""
    count = len(positions)
    if count == 0:
        return [0.0] * 5
    if count == 1:
        return [float(positions[0])] * 5
    if count == 2:
        return [float(positions[0])] * 4 + [float(positions[-1])]

    def position_at(fraction: float) -> float:
        one_based_index = round(fraction * count - 0.1)
        return float(positions[one_based_index - 1])

    return [float(positions[0]), position_at(0.25), position_at(0.5),
            position_at(0.75), float(positions[-1])]

def calculate_ctd_features(sequence: str) -> dict[str, float]:
    """Calculate the ordered descriptor vector expected by the pretrained model."""
    features: dict[str, float] = {}
    for property_name, groups in CTD_GROUPS.items():
        for group_number, amino_acids in enumerate(groups, start=1):
            positions = [index for index, residue in enumerate(sequence, start=1) if residue in amino_acids]
            for percentile, value in zip(PERCENTILES, _distribution(positions)):
                features[f"{property_name}.{group_number}.{percentile}"] = value / len(sequence) * 100.0
    return features

class AmPEPpyPredictor(BasePredictor):
    name = "amPEPpy"

    def __init__(self) -> None:
        self._model = None

    def availability(self) -> tuple[bool, str]:
        if not ((3, 10) <= sys.version_info[:2] < (3, 13)):
            return False, (
                f"The bundled model requires Python 3.10–3.12; this deployment uses "
                f"Python {sys.version_info.major}.{sys.version_info.minor}."
            )
        if not DEFAULT_MODEL.is_file():
            return False, f"Bundled model is missing: {DEFAULT_MODEL}"
        try:
            import numpy  # noqa: F401
            import sklearn  # noqa: F401
        except ImportError:
            return False, "Install NumPy and scikit-learn from requirements.txt."
        return True, str(DEFAULT_MODEL)

    def _load_model(self):
        if self._model is None:
            try:
                with DEFAULT_MODEL.open("rb") as handle:
                    self._model = pickle.load(handle)
            except Exception as exc:
                raise PredictorUnavailable(f"Could not load the bundled amPEPpy model: {exc}") from exc
        return self._model

    def predict(self, sequence: str) -> PredictionResult:
        available, detail = self.availability()
        if not available:
            raise PredictorUnavailable(detail)
        import numpy as np

        feature_values = list(calculate_ctd_features(sequence).values())
        frame = np.asarray([feature_values], dtype=float)
        try:
            probabilities = self._load_model().predict_proba(frame)[0]
        except Exception as exc:
            raise PredictorUnavailable(f"amPEPpy inference failed: {exc}") from exc
        non_amp_probability, amp_probability = map(float, probabilities)
        prediction = "AMP" if amp_probability >= non_amp_probability else "Non-AMP"
        return PredictionResult(self.name, prediction, amp_probability,
                                {"non_amp_probability": non_amp_probability})
