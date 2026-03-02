import asyncio
from fastapi import APIRouter
import state

router = APIRouter()


@router.post("/cancel")
async def cancel():
    cancelled = []

    if state.current_process is not None and state.current_process.poll() is None:
        state.current_process.kill()
        state.current_process = None
        cancelled.append("process")

    if state.current_task is not None and not state.current_task.done():
        state.current_task.cancel()
        state.current_task = None
        cancelled.append("task")

    return {"status": "ok", "cancelled": cancelled}
