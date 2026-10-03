"""Run L-0073R cluster certificate honesty audit."""

from __future__ import annotations

import json

from ns_exploration.terminal_weighted.l0073r_audit import audit_full_repair, write_audit_report


def main() -> dict:
    rep = audit_full_repair(prec=128)
    paths = write_audit_report(rep)
    rep["outputs"] = {"json": str(paths[0]), "md": str(paths[1])}
    print(json.dumps(rep, indent=2))
    return rep


if __name__ == "__main__":
    main()
