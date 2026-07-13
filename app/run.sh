#!/usr/bin/env bash
# Start the regulatr webapp. No dependencies beyond Python 3.8+ (stdlib only).
# Usage:  ./run.sh          (serves on http://127.0.0.1:8765)
#         PORT=9000 ./run.sh
set -euo pipefail
cd "$(dirname "$0")"
export PORT="${PORT:-8765}"
echo "regulatr → http://127.0.0.1:${PORT}"
exec python3 server.py
