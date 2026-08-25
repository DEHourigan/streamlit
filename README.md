# Peptide analysis platform

A modular Streamlit application for peptide database browsing, BLAST similarity
searches, physicochemical sequence analysis, and optional AMP prediction.

## Run locally

NCBI BLAST+ must be installed and `blastp` must be available on `PATH`.

```bash
python -m pip install -r requirements.txt
streamlit run src/app.py --server.address 127.0.0.1 --server.port 8501
```

The installation command must use the same Python environment as the `streamlit`
command. A reliable alternative is to launch Streamlit through Python:

```bash
python -m pip install -r requirements.txt
python -m streamlit run src/app.py
```

## Deploy on Streamlit Community Cloud

Push the complete repository to GitHub, including `requirements.txt`, `packages.txt`,
`config/databases.toml`, and `data/databases/`. In Streamlit Community Cloud create
an app with `src/app.py` as the main file. The deployment service will install:

- Python packages such as Biopython from `requirements.txt`.
- The NCBI BLAST+ system package from `packages.txt`.

The bundled BAGEL FASTA and BLAST index files are read directly from the repository,
so this database does not need a separate server. Check that Git/GitHub contains all
of the moved files before deployment:

```bash
git status
git add requirements.txt packages.txt config/databases.toml data/databases src
```

Do not commit optional predictor environments or large model dependencies merely to
make the core app deploy. Macrel and amPEPpy remain optional and report their absence
without disabling BLAST, browsing, or sequence analysis.

Macrel and amPEPpy are optional. Their dependencies are isolated from the core UI;
missing predictors are reported without breaking the other tools. See
`requirements-amp.txt` for the optional installation boundary.

Macrel's maintained installation route is Bioconda (`conda install -c bioconda
macrel`). amPEPpy can be installed from its source repository and also needs its
pretrained model file (see below).

## Add a database

Place its FASTA and BLAST+ index files under `data/databases/<key>/`, then add an
enabled `[databases.<key>]` entry to `config/databases.toml`. No Python code changes
are needed. The `blast_db` value is the index prefix, not a filename extension.

```toml
[databases.example]
name = "Example peptides"
description = "Description shown in the application."
fasta = "data/databases/example/example.faa"
blast_db = "data/databases/example/example"
enabled = true
version = "2026-01"
citation = "Relevant citation"
```

Build an index when starting from FASTA:

```bash
makeblastdb -in data/databases/example/example.faa \
  -dbtype prot -out data/databases/example/example
```

## Add a predictor

Implement `BasePredictor` from `src/predictors/base.py`, including `availability()`
and `predict(sequence)`, then register the wrapper in `pages/amp_prediction.py`.
Model imports and command execution belong in the wrapper, not in the page.

Macrel is integrated through its command-line interface. amPEPpy uses its `ampep`
command and expects the pretrained model at `data/models/ampeppy/amPEP.model`.
