# Real datasets to illustrate the multi-treatment MVBCF

The model needs, at minimum: a multivariate (or single) outcome, **two or more
crossed binary treatments**, unit covariates that can moderate the effects, and —
for the interaction to be identified — **overlap across the treatment cells**
(ideally a genuine factorial design, so all $2^T$ combinations occur). Below are
real datasets grouped by how directly they map onto the model, starting with the
one you can run in five minutes.

## A. Ready-to-run in R (no download)

**`datasets::npk`** — the canonical worked example. A classic agricultural
factorial: nitrogen (N), phosphate (P) and potassium (K), each applied or not, on
pea yield, in 6 blocks (24 plots). This is *exactly* a $T=3$ binary-treatment
factorial. Map `N,P,K → Z1,Z2,Z3`, `yield → response`, `block → covariate`, fit
with `order = 2` (or `3`), and read off the main effects and the N×P / N×K / P×K
interactions. Small, but perfect for a first end-to-end demonstration and unit
test.

**`MASS::oats`** — a split-plot with **variety × nitrogen (4 rates)**. Recode the
nitrogen rates into ordered binary contrasts (or dummy the levels) to get several
crossed "treatments"; `Block` and `Variety` are covariates/moderators. Good for
showing the ordinal-treatment expansion described in the theory (Remark on
non-binary treatments).

## B. Agricultural multi-environment factorials (closest to the drought pipeline)

**`agridat` package** (CRAN) — a large curated collection of *real* published
agronomic experiments, many of them factorial and multi-environment. Directly
relevant families:

- **Irrigation × fertilizer / nitrogen trials** (e.g. water regime × N rate across
  sites) — the natural two-treatment MET analogue of your drought pipeline: `Z1 =
  water/drought regime`, `Z2 = N (or inoculant)`, `env = site-year`, soil/climate
  covariates as moderators, multiple traits as the multivariate outcome.
- **Split-plot and factorial yield trials** (`agridat::gomez.*`,
  `agridat::yates.*`, several `.splitplot` / `.strip` datasets) with two crossed
  factors and blocking — ideal to exercise `fit_mvbcf2` with real, messy,
  unbalanced data.
- Many entries already carry genotype × environment structure, so the
  gen/env grid + counterfactual prediction of never-grown cells transfers directly.

**FACE / AGFACE experiments** (Free-Air CO₂ Enrichment; data via the hosting
institutions, e.g. AGFACE wheat) — factorial **CO₂ × irrigation × N** on wheat
traits. A real $2\times2$ (or higher) climate-treatment factorial with strong
expected interactions (CO₂ benefit depends on water/N), which is precisely what
$\tau_{12}$ is meant to capture.

**CIMMYT / ICARDA drought-management trials** (public breeding-trial repositories)
— drought stress × agronomic management across many environments, multi-trait.
The scale and unbalance match your `innovar_drought_bcf` setting.

## C. Clinical 2×2 factorial trials (canonical two-treatment causal examples)

These are textbook *randomised* factorials, so overlap is guaranteed and the
interaction estimand is clean — excellent for validating the causal machinery
against a known randomised benchmark.

- **ISIS-2** — aspirin × streptokinase, $2\times2$ factorial in acute myocardial
  infarction. The classic demonstration that two treatments and their combination
  can be estimated jointly.
- **Physicians' Health Study** — aspirin × beta-carotene, $2\times2$ factorial.
- **Women's Health Study** — aspirin × vitamin E, $2\times2$ factorial.

For these, treatment effects on the (possibly multivariate) clinical endpoints map
to $\tau_1,\tau_2$ and the drug–drug interaction to $\tau_{12}$; baseline patient
covariates moderate the CATEs. Access is via the trial data-sharing portals
(e.g. NHLBI BioLINCC) subject to their data-use terms.

## D. Social / marketing multi-treatment experiments

- **Hillstrom "MineThatData" e-mail challenge** — a three-arm marketing experiment
  (men's e-mail / women's e-mail / no e-mail) with rich customer covariates and
  spend/visit/conversion outcomes. Recode the arms into two binary treatment
  contrasts to study heterogeneous and interacted effects. Widely used in the
  uplift/causal-ML literature and freely downloadable.
- **Gerber & Green get-out-the-vote experiments** — factorial mobilisation
  treatments (mail × phone × canvassing) with turnout outcomes; multi-arm designs
  that fit the order-$r$ framework.
- **Microcredit / development multi-arm RCTs** (several public replication
  packages) — multiple crossed interventions with household covariates.

## Practical checklist when bringing your own data

1. **Overlap first.** Cross-tabulate the treatment cells within covariate strata.
   $\tau_S$ is only identified where $\prod_{t\in S}Z_t$ varies given $\bx$;
   otherwise the interaction is prior-driven (the half-scale prior shrinks it to 0,
   which is the safe default but should be reported as *unsupported*, not *null*).
2. **Randomised vs observational.** If randomised (B–D), the propensity controls
   are harmless and effects are unbiased by design. If observational (much of B),
   the per-treatment propensity scores in $\mu$ are doing real work — check
   covariate balance and positivity.
3. **Multivariate outcomes** share strength through $\Sigma$; include correlated
   traits (yield, grain weight, height; or multiple clinical endpoints) rather than
   collapsing to one.
4. **Start at `order = 2`.** Go to higher order only where the design populates the
   higher cells.
