import time
import json
import asyncio
from typing import Dict
from fastapi import APIRouter
from fastapi.responses import StreamingResponse
from state import progress_queues

router = APIRouter()


async def send_progress(job_id, stage, message, progress=None, data=None):
    """Send progress update to SSE stream (async version)"""
    if job_id in progress_queues:
        event_data = {
            "stage": stage,
            "message": message,
            "progress": progress,
            "timestamp": time.time(),
            "data": data,
        }
        await progress_queues[job_id].put(event_data)
        print(f"📤 Progress sent: {job_id} - {stage} - {message}", flush=True)


@router.get("/progress/{job_id}")
async def progress_stream(job_id: str):
    """SSE endpoint that streams progress updates"""
    
    # Create queue if doesn't exist
    if job_id not in progress_queues:
        progress_queues[job_id] = asyncio.Queue()
    
    q = progress_queues[job_id]
    
    async def generate():
        print(f"🎯 [SSE] Starting stream for job {job_id}", flush=True)
        
        yield f"data: {json.dumps({'stage': 'connected', 'message': 'Processing large files may take some time, please be patient.'})}\n\n"

        message_count = 0
        try:
            while True:
                try:
                    data = await asyncio.wait_for(q.get(), timeout=1.0)
                    
                    message_count += 1
                    print(f"📡 [SSE] Message #{message_count}: {data}", flush=True)
                    yield f"data: {json.dumps(data)}\n\n"

                    if data["stage"] in ["complete", "error"]:
                        print(f"✅ [SSE] Job {job_id} finished. Total messages: {message_count}", flush=True)
                        asyncio.create_task(cleanup_queue(job_id))
                        break

                except asyncio.TimeoutError:
                    print(f"⏰ [SSE] Timeout (messages so far: {message_count})", flush=True)
                    yield ": keepalive\n\n"

        except asyncio.CancelledError:
            print(f"🔌 [SSE] Client disconnected from job {job_id}", flush=True)
            raise

    return StreamingResponse(
        generate(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache, no-transform",
            "X-Accel-Buffering": "no",
            "Connection": "keep-alive",
        },
    )


async def cleanup_queue(job_id: str, delay: float = 5.0):
    """Clean up queue after job completes"""
    await asyncio.sleep(delay)
    progress_queues.pop(job_id, None)
    print(f"🗑️ Cleaned up queue for {job_id}", flush=True)
