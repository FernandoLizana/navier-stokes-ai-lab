#!/usr/bin/env python3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "python"))
exec((ROOT / "scripts" / "verify_c0008_cluster_certificate.py").read_text())
