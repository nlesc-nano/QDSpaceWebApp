# QDSpace WebApp

A lightweight web application for **building** and **browsing** quantum dots (QDs).  
The project combines a static frontend (served from `/docs`) with a small FastAPI backend (in `/backend`) for structure processing and interactive plotting.

> **Repo:** nlesc-nano/QDSpaceWebApp  
> **License:** MIT

---

## ✨ Features

- **Web UI** for QD building (`/docs/builder/`) and library browsing (`/docs/library/`).
- **Active navigation highlighting** via `docs/assets/nav-active.js` (works under subpaths & GitHub Pages).
- **Quantum‑dot structure library** under `docs/library/` (II–VI, III–V, IV–VI, etc.).
- **Backend APIs** (FastAPI) for generating/transforming structures and interactive plots.
- **Utility scripts** to (re)generate the library index and metadata.

---

## 🗂️ Repository Layout

```
QDSpaceWebApp/
├── LICENSE
├── README.md
├── environment.yml
├── backend/
│   ├── api.py               # API endpoints (builder/library helpers)
│   ├── api_builder.py       # Builder-specific endpoints
│   ├── app.py               # FastAPI app entry (imported by uvicorn)
│   ├── plot_interactive.py  # Exports Plotly HTML
│   ├── latest_version.txt
│   ├── backup.txt
│   └── requirements.txt
├── docs/
│   ├── index.html           # Landing page
│   ├── about.html
│   ├── contact.html
│   ├── assets/
│   │   ├── logos/           # Logos used in the site
│   │   └── nav-active.js    # Robust "active" underline for header nav
│   ├── builder/
│   │   ├── index.html       # QD builder UI
│   │   └── index_old.html
│   └── library/             # Structure explorer UI + data
│       ├── II-VI/
│       ├── II-VI@II-VI/
│       ├── III-V/
│       ├── IV-VI/
│       ├── app.css
│       ├── app.js
│       ├── file_list.js
│       ├── index.html
│       └── metadata.json
├── make_file_list.py        # Rebuild docs/library/file_list.js
└── make_metadata.py         # Rebuild docs/library/metadata.json
```

> Tip: If you also see legacy folders like `docs/II-VI/` at the repository root, they are deprecated. The canonical home of all structures is **`docs/library/`**.

---

## 🚀 Quickstart (Local)

### 1) Create the Python environment

**Option A — Conda (recommended)**

```bash
mamba env create -f environment.yml  # or: conda env create -f environment.yml
conda activate qdspace
```

**Option B — venv + pip**

```bash
python -m venv .venv
source .venv/bin/activate  # (Windows: .venv\Scripts\activate)
pip install -r backend/requirements.txt
```
**Option C — Docker (Local)**
Run everything inside a lightweight local container (no Conda/Python setup needed).

```bash
# From repo root
docker build -f Dockerfile.local -t qdspace-local .

# Run backend + serve docs, mounting local /docs
docker run --rm -p 8000:8000 -v "$PWD/docs:/app/docs" qdspace-local
# On Windows PowerShell, use:
docker run --rm -p 8000:8000 -v "${PWD}\docs:/app/docs" qdspace-local
```

Then open → [**http://127.0.0.1:8000**](http://127.0.0.1:8000)
(from there you can access the **QD Library** and **QD Builder** pages)

### 2) Run the backend (FastAPI)

From the repo root:

```bash
# Preferred (app.py exposes FastAPI instance named "app")
uvicorn backend.app:app --reload --port 8000

# If you wired the app in api.py instead
# uvicorn backend.api:app --reload --port 8000
```

Backend will be available at: `http://127.0.0.1:8000`  
Interactive docs: `http://127.0.0.1:8000/docs`

### 3) Open the frontend

You can simply open the HTML files directly or serve them via a tiny HTTP server:

```bash
# from repo root
python -m http.server 8080
# now visit: http://127.0.0.1:8080/docs/index.html
```

> The builder and library pages expect to reach your backend at `http://127.0.0.1:8000`.  
> If you change ports or run via a tunnel, adjust the fetch URLs accordingly.

---

## 📚 Updating the Library Index

When you add/remove structures under `docs/library/*`, regenerate the list and metadata:

```bash
# builds docs/library/file_list.js
python make_file_list.py

# builds docs/library/metadata.json
python make_metadata.py
```

Commit the changes so they show up in the web UI.

---

## 🌐 Deploying on GitHub Pages

1. Go to **Settings → Pages**.
2. Source: **Deploy from a branch**.  
   Branch: **`main`**, Folder: **`/docs`**.
3. Save. Your site will be published at:  
   `https://<org-or-user>.github.io/QDSpaceWebApp/`

The nav underline script `docs/assets/nav-active.js` is path-agnostic, so highlighting works correctly both locally and on Pages subpaths.

---

## 🧪 API Notes (quick peek)

- `backend/app.py` exposes the FastAPI app: `app = FastAPI(...)`
- `backend/api.py` and `backend/api_builder.py` register route handlers.
- `backend/plot_interactive.py` can export Plotly HTML for interactive PDOS/COOP/fuzzy plots.

See the in-code docstrings and the live **OpenAPI** at `/docs` when the backend is running.

---

## 🤝 Contributing

PRs and issues are welcome. If you’re adding many structures, please:
- Put them under the correct family in `docs/library/`.
- Run the two index scripts (`make_file_list.py`, `make_metadata.py`) before committing.
- Keep filenames/descriptors consistent.

---

## 📄 License

This project is licensed under the **MIT License**. See [`LICENSE`](LICENSE) for details.
