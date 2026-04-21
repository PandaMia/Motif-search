# Motif-search

FastAPI service for motif discovery in DNA sequences. It exposes a UI and API to run Greedy, Randomized, Gibbs, Genetic, or EM-based motif search and returns best motifs per gene along with consensus and scores.

Method descriptions: see [METHODS.md](METHODS.md).

## What it does
- Accepts a list of gene sequences and motif length `k`.
- Lets you choose search method and metric.
 - Supports top-level `n_iter` and `n_starts`, plus grouped `genetic_params` for genetic-only controls.
- Returns `best_motifs`, `scores`, and `consensus`.
- UI includes dynamic preset loading from `data/*.json`, where each preset stores metadata plus a `DNA` field with the sequences.

## API
- `POST /find-motif` with JSON payload matching `FindMotif`.
- Response model: `FindMotifResponse`.

### Example request
```json
{
  "genes": ["ACGTACGT", "TACGTACG"],
  "k": 3,
  "method": "em",
  "metric": "hamming",
  "n_iter": 500,
  "n_starts": 20,
  "genetic_params": {
    "n_bots": 100,
    "n_surv": 20,
    "mutation_coef": 0.15
  }
}
```

### Example response
```json
{
  "best_motifs": ["ACG", "ACG"],
  "scores": {"hamming_distance": 0, "entropy_score": 0.0},
  "consensus": "ACG"
}
```

## Run locally
From the repo root:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python app.py
```

The service starts on `http://localhost:8080` by default.

### Custom port
```bash
PORT=8000 python app.py
```

## Project structure
- `app.py` - FastAPI app and routes
- `config/base_models.py` - Pydantic models (`FindMotif`, `FindMotifResponse`)
- `motif_search/` - search algorithm implementations, including an EM/PWM motif finder
- `templates/index.html` - UI
- `static/style.css` - UI styles
- `data/*.json` - preset metadata and `DNA` sequence payloads used by the UI
