"""Dependency-safe wrapper around the amPEPpy command-line predictor."""
import csv
from pathlib import Path
import shutil
import subprocess
import tempfile
from .base import BasePredictor, PredictionResult, PredictorUnavailable

DEFAULT_MODEL = Path(__file__).resolve().parents[2] / "data" / "models" / "ampeppy" / "amPEP.model"

class AmPEPpyPredictor(BasePredictor):
    name = "amPEPpy"
    def availability(self) -> tuple[bool, str]:
        executable = shutil.which("ampep")
        if not executable:
            return False, "Install amPEPpy and expose the `ampep` command on PATH."
        if not DEFAULT_MODEL.is_file():
            return False, f"Place the pretrained model at {DEFAULT_MODEL}."
        return True, executable
    def predict(self, sequence: str) -> PredictionResult:
        available, detail = self.availability()
        if not available:
            raise PredictorUnavailable(detail)
        with tempfile.TemporaryDirectory(prefix="ampeppy-") as directory:
            fasta, output = Path(directory) / "input.fasta", Path(directory) / "prediction.tsv"
            fasta.write_text(f">query\n{sequence}\n", encoding="utf-8")
            completed = subprocess.run([detail, "predict", "-m", str(DEFAULT_MODEL), "-i", str(fasta), "-o", str(output), "--seed", "2012"], capture_output=True, text=True, check=False, timeout=300)
            if completed.returncode or not output.is_file():
                raise PredictorUnavailable(completed.stderr.strip() or "amPEPpy did not produce a prediction file.")
            with output.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
        if not rows:
            raise PredictorUnavailable("amPEPpy returned an empty prediction file.")
        row = rows[0]
        probability = next((float(value) for key, value in row.items() if "prob" in key.lower() and value), None)
        label = next((value for key, value in row.items() if key.lower() in {"prediction", "class", "label"}), None)
        prediction = "AMP" if (label and "amp" in label.lower() and "non" not in label.lower()) or (not label and probability is not None and probability >= 0.5) else "Non-AMP"
        return PredictionResult(self.name, prediction, probability, {"raw": row})
