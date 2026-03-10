import asyncio
import os
from pathlib import Path
from routers import motif_discovery, motif_track, sequence_enrichment, translate
from routers import position_enrichment, recluster
from routers import cluster_msa,cluster_phmm, cluster_list
from routers import preprocess,count,recount,filehandler,progress, cluster, cluster_diversity
from routers import motif_search, distance, mutation_network, data_merge, differential_analysis
from routers import cancel

import state
from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

# Paths that should NOT trigger cancellation (file I/O, progress streaming, health, upload)
_NO_CANCEL_PREFIXES = ('/health', '/api/v1/upload', '/api/v1/clear-files', '/')

class CancelMiddleware:
    """Pure ASGI middleware — passes receive/send through unchanged to avoid
    BaseHTTPMiddleware body-buffering deadlocks on large file uploads."""
    def __init__(self, app):
        self._app = app

    async def __call__(self, scope, receive, send):
        if scope["type"] == "http":
            method = scope.get("method", "")
            path = scope.get("path", "")
            if method == "POST" and not any(path.startswith(p) for p in _NO_CANCEL_PREFIXES):
                # Kill active subprocess (cutadapt, muscle, etc.)
                if state.current_process is not None and state.current_process.poll() is None:
                    state.current_process.kill()
                    state.current_process = None
                # Cancel active asyncio background task (preprocess, etc.)
                if state.current_task is not None and not state.current_task.done():
                    state.current_task.cancel()
                    state.current_task = None
        await self._app(scope, receive, send)

app = FastAPI(
    title="Fastaptamer3",
    description="Fastaptamer3",
    version="1.0.0"
)

# CORS middleware configuration
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
app.add_middleware(CancelMiddleware)

UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))
UPLOAD_DIR.mkdir(exist_ok=True)

# Include routers
app.include_router(recount.router, prefix="/api/v1", tags=["recount"])
app.include_router(count.router, prefix="/api/v1", tags=["count"])
app.include_router(preprocess.router, prefix="/api/v1", tags=["preprocess"])
app.include_router(filehandler.router, prefix="/api/v1", tags=["filehandler"])
app.include_router(progress.router, prefix="/api/v1", tags=["progress"])
app.include_router(cluster.router, prefix="/api/v1", tags=["cluster"])
app.include_router(cluster_diversity.router, prefix="/api/v1", tags=["cluster_diversity"])
app.include_router(cluster_msa.router, prefix="/api/v1", tags=["cluster_msa"])
app.include_router(cluster_phmm.router, prefix="/api/v1", tags=["cluster_phmm"])
app.include_router(cluster_list.router, prefix="/api/v1", tags=["cluster_list"])
app.include_router(recluster.router, prefix="/api/v1", tags=["recluster"])
app.include_router(position_enrichment.router, prefix="/api/v1", tags=["position_enrichment"])
app.include_router(motif_search.router, prefix="/api/v1", tags=["motif_search"])
app.include_router(motif_search.router, prefix="/api/v1", tags=["motif_search"])
app.include_router(motif_discovery.router, prefix="/api/v1", tags=["motif_discovery"])
app.include_router(motif_track.router, prefix="/api/v1", tags=["motif_track"])
app.include_router(sequence_enrichment.router, prefix="/api/v1", tags=["sequence_enrichment"])
app.include_router(translate.router, prefix="/api/v1", tags=["translate"])
app.include_router(distance.router, prefix="/api/v1", tags=["distance"])
app.include_router(mutation_network.router, prefix="/api/v1", tags=["mutation_network"])
app.include_router(data_merge.router, prefix="/api/v1", tags=["data_merge"])
app.include_router(differential_analysis.router, prefix="/api/v1", tags=["differential_analysis"])
app.include_router(cancel.router, prefix="/api/v1", tags=["cancel"])

@app.exception_handler(Exception)
async def general_exception_handler(request: Request, exc: Exception):
    return JSONResponse(
        status_code=500,
        content={"detail": f"{type(exc).__name__}: {str(exc)}"}
    )

@app.get("/")
async def root():
    return {"message": "Welcome to Fastaptamer3 Project"}

@app.get("/health")
async def health_check():
    return {"status": "healthy"}
