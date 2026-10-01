---
id: 2026-09-30-003-use-decadal-output-for-continuous-globe
timestamp: 2026-09-30T20:40:50+02:00
model: claude-opus-5-5
files_touched:
  - .gitignore
  - data/kernels/README.md
  - plans/2026-09-30-v2-convolution-explorer-plan.md
  - scripts/export_pulse_kernels.py
---

In step 3 left: globe - it does not need to keep the 1000-yr slices. the dye tracer model provides decadal output, you can use this to compute a quasi continuous surface field.
