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

Select Python 3.12 for the deployment. This has a compatible prebuilt scikit-learn
wheel, and `requirements.txt` pins scikit-learn 1.4.0 to the version used to serialize the
bundled pretrained model. Python cannot be changed on an existing Streamlit Community
Cloud deployment: delete and redeploy the app if it was created with another version.
The NumPy and scikit-learn requirements are guarded for Python 3.10–3.12 so an accidental
newer deployment fails fast instead of attempting a slow incompatible source build.

The bundled BAGEL FASTA and BLAST index files are read directly from the repository,
so this database does not need a separate server. Check that Git/GitHub contains all
of the moved files before deployment:

```bash
git status
git add requirements.txt packages.txt config/databases.toml data/databases src
```

The official amPEPpy pretrained model is stored at
`data/models/ampeppy/amPEP.model`. A lightweight local adapter calculates its CTD
features and runs that model directly, avoiding a legacy Git package build during
deployment. Commit the model file so visitors can use the classifier. Predictor
execution is isolated from the UI, so a model-specific runtime error does not disable
the other tools.

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

The hosted app runs the bundled amPEPpy model directly through the adapter in
`src/predictors/ampeppy.py`. A dormant Macrel wrapper remains available
for installations with a separately managed Macrel environment, but it is not part
of the lightweight Streamlit Cloud deployment.
