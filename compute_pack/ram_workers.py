"""Recommend parallel shell workers from available RAM and CPU (no GPU needed)."""

from __future__ import annotations

import os
import sys

# Match repair_full_manifest.py
RAM_GB_PER_WORKER = 1.75
OS_RESERVE_GB = 2.0


def available_ram_gb() -> float:
    if sys.platform == "win32":
        import ctypes

        class MEMORYSTATUSEX(ctypes.Structure):
            _fields_ = [
                ("dwLength", ctypes.c_ulong),
                ("dwMemoryLoad", ctypes.c_ulong),
                ("ullTotalPhys", ctypes.c_ulonglong),
                ("ullAvailPhys", ctypes.c_ulonglong),
                ("ullTotalPageFile", ctypes.c_ulonglong),
                ("ullAvailPageFile", ctypes.c_ulonglong),
                ("ullTotalVirtual", ctypes.c_ulonglong),
                ("ullAvailVirtual", ctypes.c_ulonglong),
                ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
            ]

        stat = MEMORYSTATUSEX()
        stat.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
        if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(stat)):
            return stat.ullAvailPhys / (1024**3)
    else:
        try:
            with open("/proc/meminfo", encoding="utf-8") as f:
                for line in f:
                    if line.startswith("MemAvailable:"):
                        return int(line.split()[1]) / (1024**2)
        except OSError:
            pass
    try:
        import psutil  # type: ignore

        return psutil.virtual_memory().available / (1024**3)
    except Exception:
        return 8.0


def recommend_workers(*, requested_max: int = 32) -> int:
    """Cap by RAM, CPU cores, and requested ceiling."""
    cpu = os.cpu_count() or 4
    ram_cap = max(1, int((available_ram_gb() - OS_RESERVE_GB) // RAM_GB_PER_WORKER))
    cpu_cap = max(1, cpu - 1)
    return max(1, min(requested_max, ram_cap, cpu_cap))


def main() -> None:
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--max", type=int, default=32, help="ceiling (use 32+ on high-RAM boxes)")
    p.add_argument("--json", action="store_true")
    args = p.parse_args()
    ram = available_ram_gb()
    w = recommend_workers(requested_max=args.max)
    if args.json:
        import json

        print(
            json.dumps(
                {
                    "available_ram_gb": round(ram, 2),
                    "cpu_count": os.cpu_count(),
                    "recommended_workers": w,
                    "ram_gb_per_worker": RAM_GB_PER_WORKER,
                }
            )
        )
    else:
        print(w)


if __name__ == "__main__":
    main()
