import mimetypes
from fastapi import APIRouter, File, UploadFile, HTTPException
from pathlib import Path
import os
from fastapi.responses import FileResponse
import shutil
from datetime import datetime

router = APIRouter()
UPLOAD_DIR = Path(os.getenv("UPLOAD_DIR", "files"))

@router.post("/upload")
async def upload(file: UploadFile = File(...)):

    try:
        # timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Split filename and extension
        file_name = Path(file.filename).stem
        file_extension = Path(file.filename).suffix
        
        # new_filename = f"{file_name}_{timestamp}{file_extension}"
        new_filename = f"{file_name}{file_extension}"
        file_path = UPLOAD_DIR / new_filename

        # Save the file
        with open(file_path, "wb") as buffer:
            shutil.copyfileobj(file.file, buffer)
        
        return {
            "original_filename": file.filename,
            "saved_filename": new_filename,
            "message": "File uploaded successfully",
        }
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"File upload failed: {str(e)}")
    
    finally:
        file.file.close()


@router.get("/download/{filename}")
async def download(filename: str):
    file_path = UPLOAD_DIR / filename
    
    # Check if file exists
    if not file_path.exists():
        raise HTTPException(status_code=404, detail="File not found")
    
    media_type, _ = mimetypes.guess_type(str(file_path))
    if media_type is None:
        media_type = 'application/octet-stream'
    
    # Return the file with download headers
    return FileResponse(
        path=str(file_path),
        filename=filename,
        media_type=media_type,
        headers={
            "Content-Disposition": f'attachment; filename="{filename}"'
        }
    )

def get_media_type(filename: str) -> str:
    """Determine media type based on file extension"""
    extension = filename.lower().split('.')[-1]
    
    media_types = {
        'fastq': 'text/plain',
        'fq': 'text/plain',
        'fasta': 'text/plain',
        'fa': 'text/plain',
        'csv': 'text/csv',
        'tsv': 'text/tab-separated-values',
        'txt': 'text/plain',
        'gz': 'application/gzip',
        'zip': 'application/zip',
        'pdf': 'application/pdf',
    }
    
    return media_types.get(extension, 'application/octet-stream')

@router.get("/files/")
async def list_files():
    files = [f.name for f in UPLOAD_DIR.iterdir() if f.is_file()]
    return {"files": files, "count": len(files)}


@router.delete("/delete/{filename}")
async def delete_file(filename: str):
    file_path = UPLOAD_DIR / filename
    if not file_path.exists():
        raise HTTPException(status_code=404, detail="File not found")
    
    try:
        os.remove(file_path)
        return {"message": f"File '{filename}' deleted successfully"}
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to delete file: {str(e)}")
