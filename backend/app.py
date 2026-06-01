# backend/app.py
# Active trigger for uvicorn auto-reload of nanocrystal-builder optimizations
from pathlib import Path
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from starlette.staticfiles import StaticFiles

# Import API modules
from . import api as library_api        # backend/api.py  (miniCAT-based)
from . import api_builder as builder_api  # backend/api_builder.py (QD_Builder-based)

# Initialize FastAPI application
app = FastAPI(title="Quantum Dot Suite")

# Configure CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],   # you can restrict later to e.g. ["http://localhost:8000"]
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Health check endpoint
@app.get("/healthz")
def healthz():
    return {"status": "ok"}

# Mount the Library API at /api
app.mount("/api", library_api.app)

# Mount the Builder API at /builder/api
app.mount("/builder", builder_api.app)

from starlette.responses import FileResponse

# Serve static frontend (Vite compiled output in qd-frontend/dist)
# Only mount when running locally — in Lambda, qd-frontend/ is excluded
# by .dockerignore since CloudFront/S3 serves the frontend.
frontend_dist = Path(__file__).resolve().parent.parent / "qd-frontend" / "dist"

if frontend_dist.exists():
    # Mount /assets specifically so JS/CSS loads fast via StaticFiles
    app.mount("/assets", StaticFiles(directory=str(frontend_dist / "assets")), name="static-assets")

    # SPA Fallback logic for everything else
    @app.get("/{full_path:path}")
    async def serve_spa(full_path: str):
        target_path = frontend_dist / full_path
        if target_path.is_file():
            return FileResponse(target_path)
        return FileResponse(frontend_dist / "index.html")