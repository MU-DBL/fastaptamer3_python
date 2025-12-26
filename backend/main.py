import os
from pathlib import Path
from routers import motif_discovery, motif_track, sequence_enrichment, translate
from routers import position_enrichment, recluster
from routers import cluster_msa,cluster_phmm, cluster_list
from routers import preprocess,count,recount,filehandler,progress, cluster, cluster_diversity
from routers import motif_search

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

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

@app.get("/")
async def root():
    return {"message": "Welcome to Fastaptamer3 Project"}

@app.get("/health")
async def health_check():
    return {"status": "healthy"}
