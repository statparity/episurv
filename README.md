# episurv

**Python library for infectious disease surveillance — EpiEstim / EpiNow2 / surveillance parity.**

[![License](https://img.shields.io/badge/license-Apache--2.0-blue)](LICENSE)

---

## Motivation

The canonical Python tools for infectious disease surveillance are incomplete relative to R's `EpiEstim` (Rt estimation), `EpiNow2` (Rt + nowcasting + reporting delays), and `surveillance` (prospective/retrospective aberration detection). `episurv` closes this gap with direct feature parity, enabling Python-native public health workflows.

---

## Scope (Target Parity)

| Capability | R package | `episurv` | Target |
|---|---|---|---|
| Instantaneous Rt (renewal equation) | `EpiEstim` | `episurv.rt.instant` | v0.1.0 |
| Rt with uncertain serial interval | `EpiEstim` | `episurv.rt.uncertain_si` | v0.2.0 |
| Rt + nowcasting + reporting delays | `EpiNow2` | `episurv.rt.epinow` | v0.4.0 |
| Prospective aberration detection (CUSUM) | `surveillance` | `episurv.detect.cusum` | v0.3.0 |
| Prospective aberration detection (EARS) | `surveillance` | `episurv.detect.ears` | v0.3.0 |
| Retrospective outbreak detection (scan) | `surveillance` | `episurv.detect.scan` | v0.5.0 |
| SIR/SEIR compartmental models | `EpiModel` | `episurv.models.compartmental` | v0.2.0 |
| Age-stratified contact matrices | `socialmixr` | `episurv.contacts` | v0.6.0 |
| R parity validation suite | — | `tests/test_vs_r/` | v0.7.0 |

---

## License

Apache-2.0 — see [LICENSE](LICENSE).
