"""Small dependency-free FASTA reader and database summaries."""
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, median, pstdev

@dataclass(frozen=True)
class FastaRecord:
    id: str
    description: str
    sequence: str

def read_fasta(path: Path) -> list[FastaRecord]:
    records, sequence = [], []
    header = None
    with Path(path).open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(_record(header, sequence))
                header, sequence = line[1:].strip(), []
            else:
                if header is None:
                    raise ValueError(f"Sequence encountered before a FASTA header in {path}")
                sequence.append(line)
    if header is not None:
        records.append(_record(header, sequence))
    return records

def _record(header: str, sequence: list[str]) -> FastaRecord:
    identifier, _, description = header.partition(" ")
    return FastaRecord(identifier, description or header, "".join(sequence).replace("*", "").upper())

def database_stats(records: list[FastaRecord]) -> tuple[dict[str, float | int], list[int]]:
    lengths = [len(record.sequence) for record in records]
    if not lengths:
        raise ValueError("The FASTA database contains no sequences.")
    return ({"n_sequences": len(lengths), "total_residues": sum(lengths), "mean_length": mean(lengths),
        "sd_length": pstdev(lengths), "median_length": median(lengths), "min_length": min(lengths),
        "max_length": max(lengths)}, lengths)
