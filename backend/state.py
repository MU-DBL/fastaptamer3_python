import asyncio
import subprocess
from typing import Optional, Union

progress_queues: dict[str, asyncio.Queue] = {}

current_task: Optional[asyncio.Task] = None
current_process: Optional[Union[subprocess.Popen, asyncio.subprocess.Process]] = None


def process_running(proc) -> bool:
    """True while a subprocess.Popen or asyncio subprocess has not exited."""
    if proc is None:
        return False
    if isinstance(proc, subprocess.Popen):
        return proc.poll() is None
    return proc.returncode is None
