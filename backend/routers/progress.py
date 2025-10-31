import time
import json
import asyncio
from fastapi import APIRouter
import threading
import queue
from fastapi.responses import StreamingResponse

router = APIRouter()

from state import progress_queues


@router.get("/progress/{job_id}")
async def progress_stream(job_id: str):
    """SSE endpoint for progress updates"""
    
    async def generate():
        if job_id not in progress_queues:
            yield f"data: {json.dumps({'error': 'Invalid job ID'})}\n\n"
            return
        
        q = progress_queues[job_id]
        
        # Send immediate connection confirmation
        yield f"data: {json.dumps({'stage': 'connected', 'message': 'Connected to job stream'})}\n\n"
        
        while True:
            try:
                # Non-blocking check with small timeout
                data = q.get(timeout=1)  # Changed from 30 to 1
                event = f"data: {json.dumps(data)}\n\n"
                yield event
                
                # Important: Add a small delay to force flush
                await asyncio.sleep(0.01)
                
                # Clean up if job is complete or errored
                if data['stage'] in ['complete', 'error']:
                    threading.Timer(10, lambda: progress_queues.pop(job_id, None)).start()
                    break
                    
            except queue.Empty:
                # Send keepalive more frequently
                yield f": keepalive\n\n"
                await asyncio.sleep(0.1)
    
    return StreamingResponse(
        generate(),
        media_type="text/event-stream",
        headers={
            'Cache-Control': 'no-cache, no-transform',
            'X-Accel-Buffering': 'no',
            'Connection': 'keep-alive',
            'Content-Type': 'text/event-stream; charset=utf-8',
            'Access-Control-Allow-Origin': '*',  # Add CORS
            'Access-Control-Allow-Methods': 'GET',
            'Access-Control-Allow-Headers': 'Content-Type',
        }
    )


def send_progress(job_id, stage, message, progress=None, data=None):
    """Send progress update to SSE stream"""
    if job_id in progress_queues:
        event_data = {
            'stage': stage,
            'message': message,
            'progress': progress,
            'timestamp': time.time(),
            'data': data
        }
        progress_queues[job_id].put(event_data)
        print(f"📤 Progress sent: {job_id} - {stage} - {message}")  # Debug log
