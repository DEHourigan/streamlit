"""BLAST execution and tabular result parsing."""
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
import shutil
import subprocess
import tempfile
import pandas as pd

BLAST_COLUMNS = ["Query", "Hit", "Identity (%)", "Alignment length", "Mismatches", "Gap opens", "Query start",
    "Query end", "Subject start", "Subject end", "E-value", "Bit score", "Query coverage (%)", "Query length",
    "Subject length", "Aligned query", "Aligned subject", "Description"]
OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen qseq sseq stitle"

class BlastError(RuntimeError):
    pass

@dataclass(frozen=True)
class BlastResult:
    hits: pd.DataFrame
    task: str
    command_version: str

def run_blast(sequence: str, database_prefix: Path, evalue: float = 10.0, max_hits: int = 50) -> BlastResult:
    executable = shutil.which("blastp")
    if not executable:
        raise BlastError("blastp was not found. Install NCBI BLAST+ and ensure blastp is on PATH.")
    task = "blastp-short" if len(sequence) < 30 else "blastp"
    with tempfile.TemporaryDirectory(prefix="peptide-blast-") as directory:
        query = Path(directory) / "query.fasta"
        query.write_text(f">query\n{sequence}\n", encoding="utf-8")
        command = [executable, "-task", task, "-query", str(query), "-db", str(database_prefix), "-evalue",
            str(evalue), "-max_target_seqs", str(max_hits), "-outfmt", OUTFMT]
        try:
            completed = subprocess.run(command, capture_output=True, text=True, timeout=300, check=False)
        except subprocess.TimeoutExpired as exc:
            raise BlastError("BLAST timed out after 5 minutes.") from exc
    if completed.returncode:
        raise BlastError(completed.stderr.strip() or "BLAST exited with an error.")
    hits = pd.read_csv(StringIO(completed.stdout), sep="\t", names=BLAST_COLUMNS) if completed.stdout.strip() else pd.DataFrame(columns=BLAST_COLUMNS)
    version = subprocess.run([executable, "-version"], capture_output=True, text=True, check=False).stdout.splitlines()[0]
    return BlastResult(hits, task, version)
