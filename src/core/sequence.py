"""Sequence validation and peptide-property calculations."""
import re
from collections import Counter

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
BLAST_AA = set("ABCDEFGHIKLMNPQRSTVWXYZJUO*")

def clean_sequence(sequence: str, standard_only: bool = False) -> str:
    sequence = re.sub(r"\s+", "", sequence).upper()
    if sequence.startswith(">"):
        raise ValueError("Paste the sequence only, without a FASTA header.")
    if not sequence:
        raise ValueError("Enter an amino-acid sequence.")
    invalid = set(sequence) - (STANDARD_AA if standard_only else BLAST_AA)
    if invalid:
        raise ValueError("Invalid amino-acid characters: " + ", ".join(sorted(invalid)))
    return sequence.rstrip("*")

def calculate_properties(sequence: str) -> dict:
    sequence = clean_sequence(sequence, standard_only=True)
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
    except ImportError as exc:
        raise RuntimeError("Sequence analysis requires Biopython. Install it with: pip install biopython") from exc
    analysis = ProteinAnalysis(sequence)
    composition = {aa: count / len(sequence) for aa, count in sorted(Counter(sequence).items())}
    net_charge = sequence.count("K") + sequence.count("R") + 0.1 * sequence.count("H") - sequence.count("D") - sequence.count("E")
    return {"length": len(sequence), "molecular_weight": analysis.molecular_weight(),
        "isoelectric_point": analysis.isoelectric_point(), "net_charge": net_charge,
        "hydrophobicity": analysis.gravy(), "aromaticity": analysis.aromaticity(),
        "instability_index": analysis.instability_index(), "cysteine_count": sequence.count("C"),
        "cysteine_fraction": sequence.count("C") / len(sequence), "composition": composition}
