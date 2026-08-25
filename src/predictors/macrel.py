"""Dependency-safe wrapper around the Macrel command-line predictor."""
import csv
import gzip
from pathlib import Path
import shutil
import subprocess
import tempfile
from .base import BasePredictor, PredictionResult, PredictorUnavailable

class MacrelPredictor(BasePredictor):
    name = "Macrel"
    def availability(self) -> tuple[bool, str]:
        path = shutil.which("macrel")
        return (True, path) if path else (False, "Install Macrel and expose the `macrel` command on PATH.")
    def predict(self, sequence: str) -> PredictionResult:
        available, detail = self.availability()
        if not available:
            raise PredictorUnavailable(detail)
        with tempfile.TemporaryDirectory(prefix="macrel-") as directory:
            fasta, output = Path(directory) / "input.fasta", Path(directory) / "output"
            fasta.write_text(f">query\n{sequence}\n", encoding="utf-8")
            completed = subprocess.run([detail, "peptides", "--fasta", str(fasta), "--output", str(output), "--keep-negatives"], capture_output=True, text=True, check=False, timeout=300)
            prediction_files = list(output.glob("*prediction*.gz")) + list(output.glob("*prediction*")) if output.exists() else []
            if completed.returncode or not prediction_files:
                raise PredictorUnavailable(completed.stderr.strip() or "Macrel did not produce a prediction file.")
            opener = gzip.open if prediction_files[0].suffix == ".gz" else open
            with opener(prediction_files[0], mode="rt", encoding="utf-8") as handle:
                rows = list(csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t"))
        if not rows:
            return PredictionResult(self.name, "Non-AMP", None, {"note": "Macrel returned no AMP candidate."})
        row = rows[0]
        amp_key = next((key for key in row if "amp" in key.lower() and "prob" in key.lower()), None)
        probability = float(row[amp_key]) if amp_key and row[amp_key] else None
        haemolysis = next((row[key] for key in row if "hemol" in key.lower() or "haemol" in key.lower()), None)
        prediction = "AMP" if probability is None or probability >= 0.5 else "Non-AMP"
        return PredictionResult(self.name, prediction, probability, {"haemolysis": haemolysis, "raw": row})
