import os
from pathlib import Path
from routers import preprocess,count,recount,filehandler
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

@app.get("/")
async def root():
    return {"message": "Welcome to Fastaptamer3 Project"}

@app.get("/health")
async def health_check():
    return {"status": "healthy"}
