import asyncio
import subprocess
from typing import Optional

progress_queues: dict[str, asyncio.Queue] = {}

current_task: Optional[asyncio.Task] = None
current_process: Optional[subprocess.Popen] = None
