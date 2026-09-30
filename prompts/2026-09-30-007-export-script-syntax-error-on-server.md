---
id: 2026-09-30-007-export-script-syntax-error-on-server
timestamp: 2026-09-30T21:03:58+02:00
model: claude-opus-5-5
files_touched:
  - data/kernels/README.md
  - scripts/export_pulse_kernels.py
---

an error: 

  File "export_pulse_kernels.py", line 89
    DYES = [f"dye{i:02d}" for i in range(9)]
                        ^
SyntaxError: invalid syntax
