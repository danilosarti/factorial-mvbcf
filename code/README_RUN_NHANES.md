# Running the NHANES deep application

This runs on your Mac's real R (it needs live CRAN + internet to the CDC and the
compiled MVBCF package — neither the cloud sandbox nor the device VM has R or can
reach NHANES, so it can't run there).

## One command (from the repo root)

```bash
cd ~/mvbcf_two_treatments
caffeinate -is Rscript real_data/nhanes_application.R
```

First run installs `nhanesA` (and `grf`, `dbarts` if missing) and downloads the
NHANES tables — a few minutes of network, then two MVBCF fits (a few minutes to
~an hour depending on cores). `caffeinate` keeps the Mac awake; run it in its own
Terminal window as you did for the comparison.

## What it produces (in `real_data/outputs/`)

- `nhanes_A_cardio_ATE.csv`, `nhanes_B_neuro_ATE.csv` — average effects (main +
  interaction) per outcome, SD units + natural units, 95% CrI, excludes-zero flag,
  unmeasured-confounding sensitivity value.
- `nhanes_A_cardio_grf.csv`, `nhanes_B_neuro_grf.csv` — grf multi-arm comparison.
- `nhanes_*_forest.pdf` — ATE forest plots.
- `nhanes_B_neuro_vsup_interaction.pdf` — exploratory VSUP interaction map.
- `nhanes_*.rds` — full fit objects.
- `nhanes_run_log.txt` — cell counts, overlap diagnostic, ATE tables (everything
  printed).

## Then

Send me `nhanes_run_log.txt` (or just say it finished — I'll read the committed
outputs). I'll fill the `\TBD{}` slots in `paper/real_data_section.tex` with the
real numbers, build the tables/figures into the manuscript, and recompile.

## The two analyses (all NHANES codes verified against the CDC codebooks)

**A — Cardiometabolic (adults 20+).** Exposures: statin use (from `RXQ_RX`,
drug-name match excluding nystatin) × physical activity (`PAQ605/620/650/665`).
Outcomes: systolic/diastolic BP (`BPXSY*/BPXDI*`), HbA1c (`LBXGH`), total
cholesterol (`LBXTC`). Statin→cholesterol is a positive control.

**B — Neurology/cognition (ages 60+).** Exposures: antihypertensive medication
(`BPQ050A`) × physical activity. Outcomes: CERAD immediate recall
(`CFDCST1+2+3`), CERAD delayed (`CFDCSR`), Animal Fluency (`CFDAST`),
Digit-Symbol (`CFDDS`).

Covariates (both): age, sex, education, race/ethnicity, income-to-poverty, BMI,
diabetes, smoking. Honesty guardrails are enforced in the script and written into
the paper section: ATE is the calibrated claim; interactions only where cells are
supported; CATE maps exploratory; sensitivity value on every ATE.
