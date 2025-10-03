from pathlib import Path
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from starlette.staticfiles import StaticFiles

# Import your two existing FastAPI apps
from . import api as library_api        # backend/api.py  (miniCAT-based)
from . import api_builder as builder_api  # backend/api_builder.py (QD_Builder-based)

app = FastAPI(title="Quantum Dot Suite")

# Allow CORS (so it works if you open HTML directly from file:// or localhost)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],   # you can restrict later to e.g. ["http://localhost:8000"]
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Mount the Library API at /api
app.mount("/api", library_api.app)

# Mount the Builder API at /builder/api
app.mount("/builder", builder_api.app)

# Serve static frontends (landing page, builder UI, library UI)
docs_dir = Path(__file__).resolve().parent.parent / "docs"
app.mount("/", StaticFiles(directory=str(docs_dir), html=True), name="static")

@app.get("/healthz")
def healthz():
    return {"status": "ok"}

