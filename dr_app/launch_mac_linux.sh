#!/bin/bash
# Drip Rate Estimator — macOS / Linux launcher

cd "$(dirname "$0")"

# Check Python
if ! command -v python3 &> /dev/null; then
    echo "ERROR: python3 not found. Please install Python 3.8+."
    exit 1
fi

# Check Flask
if ! python3 -c "import flask" &> /dev/null; then
    echo "Flask not found. Installing..."
    pip3 install flask
fi

echo ""
echo "  Drip Rate Estimator"
echo "  Starting server — open your browser at:  http://localhost:5000"
echo "  Press Ctrl+C to stop."
echo ""

python3 app.py
