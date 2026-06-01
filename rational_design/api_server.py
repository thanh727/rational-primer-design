from __future__ import annotations

import os

import uvicorn


def main() -> None:
    host = os.environ.get("RPD_API_HOST", "127.0.0.1")
    port = int(os.environ.get("RPD_API_PORT", "8000"))
    reload = os.environ.get("RPD_API_RELOAD", "0") == "1"
    uvicorn.run("rational_design.api:app", host=host, port=port, reload=reload)


if __name__ == "__main__":
    main()
