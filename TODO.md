# TODO — episurv

> **Status: Milestone stub.** Detailed atomization to follow after `pyrbmi` v0.3.0.

---

## Phase 0 — Setup (same pattern as pyrbmi)
- [x] 0.1 uv init, pyproject.toml, ruff, mypy
- [x] 0.2 CI/CD: lint, test, R parity weekly run
- [x] 0.3 License: Apache-2.0
- [x] 0.4 MkDocs documentation
- [x] 0.5 .gitignore
- [x] 0.6 .pre-commit-config.yaml

## v0.1.0 — Instantaneous Rt (EpiEstim parity)
- [ ] `Incidence` data object (date-indexed, right-censoring support)
- [ ] `SerialInterval`: parametric (Gamma, LogNormal) + empirical
- [ ] `estimate_rt_instant()`: renewal equation, sliding window, Gamma posterior
- [ ] `RtResult`: per-window R, lower/upper CI, posterior mean
- [ ] Validation: compare to `EpiEstim::estimate_R()` on `flu1918` dataset (tol=1e-4)

## v0.2.0 — Uncertain SI + Compartmental Models
- [ ] `estimate_rt_uncertain_si()`: bootstrap SI uncertainty (matches EpiEstim `method="uncertain_si"`)
- [ ] `SIR`, `SEIR`, `SEIRD` compartmental models
- [ ] PyMC Bayesian calibration for compartmental models
- [ ] Contact matrix loading (POLYMOD Prem et al. 2017 data)

## v0.3.0 — Prospective Surveillance (R `surveillance` parity)
- [ ] `CUSUMDetector`: CUSUM control chart for count time series
- [ ] `EARSDetector`: CDC EARS C1/C2/C3 methods
- [ ] `OutbreakAlert`: unified alert output object

## v0.4.0 — EpiNow2 Analog (Rt + Nowcasting)
- [ ] Delay distribution estimation (symptom onset → report date)
- [ ] Nowcasting: right-truncation adjustment
- [ ] Full Bayesian Rt pipeline (PyMC)
- [ ] Validation against EpiNow2 on COVID-19 England data

## v0.5.0 — Scan Statistics + Retrospective Detection
- [ ] Temporal scan statistic (Kulldorff)
- [ ] Spatial scan statistic (requires coordinates)

## v0.6.0 — Contact Matrices + Age-Stratified Models
- [ ] `ContactMatrix` from POLYMOD, CoMix data
- [ ] Age-stratified SEIR

## v0.7.0 — R Parity Validation Suite
- [ ] rpy2-based comparison: EpiEstim, EpiNow2, surveillance
- [ ] Auto-generated validation report

## v1.0.0 — Stable + JOSS Paper
