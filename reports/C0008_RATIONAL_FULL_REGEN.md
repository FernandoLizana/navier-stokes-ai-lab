# Full n=24 rational repair — launch (external PowerShell recommended)

Heavy compute: **87 shells**, checkpointed, **days** with `--workers 4`.

## Launch (recommended — survives Cursor closing)

**Do not** rely on Cursor background terminals; they die when the IDE closes.

```powershell
Set-Location "<REPO_ROOT>\python"
.\scripts\start_c0008_rational_supervisor.ps1
```

Supervisor (`supervise_c0008_rational_regen.py`):
- Spawns regen + inner watchdog in a **detached** process
- **Auto-restarts** on crash or exit before 87/87
- **Kills and restarts** if heartbeat stale >25 min (stuck)
- Log: `experiments/terminal_weighted/supervisor_n24_rational.log`

### Stop

```powershell
.\scripts\stop_c0008_rational_supervisor.ps1
```

### Manual (not recommended in Cursor)

```powershell
python scripts/regenerate_c0008_full_dealias_one_inf.py --n 24 --rational --workers 8 --adaptive --watchdog
```

## Progress

```powershell
python scripts/regenerate_c0008_full_dealias_one_inf.py --n 24 --rational --status
```

## Outputs (separate from float repair)

| File | Purpose |
|------|---------|
| `experiments/terminal_weighted/repair_one_inf_checkpoint_n24_rational.json` | Shell checkpoint |
| `experiments/terminal_weighted/shell_manifest_repair_full_one_inf_n24_rational_partial.json` | Partial manifest |
| `experiments/terminal_weighted/shell_manifest_repair_full_one_inf_n24_rational.json` | Final (when 87/87) |
| `experiments/terminal_weighted/repair_watchdog_state_n24_rational.json` | Watchdog state |
| `experiments/terminal_weighted/shell_one_inf_ck/n24_rational/shell_NNN.json` | Per-shell intra ck |

## Constraints

- **C-0008 remains `exploring`** — rational regen does not claim closure.
- Float manifest `shell_manifest_repair_full_one_inf_n24.json` (I_hi ≈ 170.59) stays **primary Route A** until rational run completes and is evaluated.
- Do not mix float and rational checkpoints.
