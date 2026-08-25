"""Database registry loading and path validation."""
from dataclasses import dataclass
from pathlib import Path
try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10
    import tomli as tomllib

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_REGISTRY = REPO_ROOT / "config" / "databases.toml"

@dataclass(frozen=True)
class DatabaseConfig:
    key: str
    name: str
    description: str
    fasta: Path
    blast_db: Path
    enabled: bool = True
    version: str = "Not specified"
    citation: str = "Not specified"

    @property
    def fasta_available(self) -> bool:
        return self.fasta.is_file()

    @property
    def blast_available(self) -> bool:
        return any(self.blast_db.parent.glob(f"{self.blast_db.name}.p*"))

def _resolve(path: str) -> Path:
    candidate = Path(path)
    return candidate if candidate.is_absolute() else REPO_ROOT / candidate

def load_databases(registry_path: Path = DEFAULT_REGISTRY, enabled_only: bool = True) -> dict[str, DatabaseConfig]:
    if not registry_path.is_file():
        raise FileNotFoundError(f"Database registry not found: {registry_path}")
    with registry_path.open("rb") as handle:
        raw = tomllib.load(handle).get("databases", {})
    databases = {}
    for key, values in raw.items():
        database = DatabaseConfig(
            key=key, name=values.get("name", key), description=values.get("description", ""),
            fasta=_resolve(values["fasta"]), blast_db=_resolve(values["blast_db"]),
            enabled=values.get("enabled", True), version=values.get("version", "Not specified"),
            citation=values.get("citation", "Not specified"),
        )
        if database.enabled or not enabled_only:
            databases[key] = database
    return databases
