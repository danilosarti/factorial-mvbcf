# ACCORD via BioLINCC — data-access checklist + draft request

Goal: obtain the ACCORD trial data to run `real_data/accord_application.R` — the
**randomized 2×2 factorial** confirmation of the interaction estimand (glycemia ×
blood-pressure), the strongest single strengthener for AOAS/JRSS-B. This is the one
item that cannot be executed from here or from your Mac until NHLBI approves access;
below is everything needed to reduce it to a single submission.

## Steps (≈ the whole process)

1. **Create a BioLINCC account** at https://biolincc.nhlbi.nih.gov (login via eRA
   Commons / Login.gov).
2. **Find the study**: search "ACCORD" (Action to Control Cardiovascular Risk in
   Diabetes). Open its Study Datasets page.
3. **IRB / ethics**: most public-use BioLINCC datasets need either IRB approval or a
   documented IRB exemption/non-human-subjects determination from RCSI. Request the
   determination now (it's usually the slow step). A secondary analysis of
   de-identified trial data is typically exempt — RCSI's REC can issue that.
4. **Complete the Data Use Agreement (DUA)**: institutional signing official at RCSI
   signs it (not you personally). Have your RCSI research office ready.
5. **Submit the request package**: DUA + IRB determination + a short research-use
   statement (draft below).
6. **Wait for approval** (typically a few weeks), then **download** the SAS/CSV
   datasets; place them under `~/data/accord/` and fill the `<<< SET >>>` variable
   mappings in `real_data/accord_application.R` from the delivered data dictionary.
7. Run: `cd ~/mvbcf_two_treatments && caffeinate -is Rscript real_data/accord_application.R`.

## Draft research-use statement (paste into the request)

> **Title:** Heterogeneous factorial treatment effects on correlated cardiometabolic
> outcomes in ACCORD.
> **Investigator:** Danilo A. Sarti, Royal College of Surgeons in Ireland (RCSI).
> **Aim:** Apply a factorial multivariate Bayesian causal forest to the randomized
> 2×2 factorial structure of ACCORD (intensive vs standard glycemia crossed with the
> blood-pressure intervention) to estimate the main effects and, in particular, the
> *interaction* of the two interventions on a panel of correlated continuous
> endpoints (e.g., HbA1c, systolic blood pressure, LDL cholesterol), together with
> effect heterogeneity across baseline covariates. Because assignment is randomized,
> the interaction is identified without unconfoundedness assumptions.
> **Outputs:** average main/interaction effects with credible and (debiased)
> confidence intervals, and uncertainty-aware effect maps. Methodological paper;
> no participant re-identification is attempted; data handled per the DUA.
> **Data security:** stored on an RCSI-managed encrypted device, not redistributed,
> destroyed at project end per the DUA.

## Alternative if ACCORD access stalls
The Women's Health Initiative (WHI) hormone-therapy × dietary-modification ×
calcium/vitamin-D design (also BioLINCC) gives a randomized multi-factor structure;
the same script applies with the WHI arm variables. For a fully-open (no-DUA) fallback,
the paper already ships the NHANES observational application with the debiased
estimator — the randomized version is an upgrade, not a prerequisite for submission.

## Live links (verified Sep 2026)
- ACCORD study page: https://biolincc.nhlbi.nih.gov/studies/accord/
- Specimen & Data Request form: https://biolincc.nhlbi.nih.gov/requests/specimen-and-data-request/form/
- BioLINCC FAQ (access process): https://www.biolincc.nhlbi.nih.gov/faq/
