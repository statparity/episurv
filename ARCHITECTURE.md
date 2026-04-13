# ARCHITECTURE — episurv

> **Status: Architecture stub.** Full design to be elaborated after `pyrbmi` v0.3.0.

---

## 1. Package Structure (Planned)

```
src/episurv/
├── rt/
│   ├── instant.py       # EpiEstim: instantaneous Rt via renewal equation
│   ├── uncertain_si.py  # Rt with bootstrapped serial interval uncertainty
│   └── epinow.py        # EpiNow2 analog: Rt + nowcasting + delay adjustment
│
├── detect/
│   ├── cusum.py         # CUSUM prospective surveillance
│   ├── ears.py          # EARS (CDC Early Aberration Reporting System)
│   └── scan.py          # Scan statistic (Kulldorff spatial/temporal)
│
├── models/
│   ├── compartmental.py # SIR, SEIR, SEIRD with Bayesian calibration (PyMC)
│   └── network.py       # Network-based transmission (EoN analog)
│
├── contacts/
│   └── matrices.py      # Age-stratified contact matrices (POLYMOD/CoMix data)
│
└── data/
    ├── incidence.py     # Incidence object (daily/weekly, right-censored)
    └── linelist.py      # Individual-level case data (onset, report dates)
```

---

## 2. Core Algorithm: Instantaneous Rt (EpiEstim)

The instantaneous reproduction number R(t) is estimated via the renewal equation:

```
I(t) = R(t) × Σ_{s=1}^{t} I(t-s) × w(s)
```

Where `w(s)` is the serial interval distribution (discrete probability mass).

Assuming `R(t)` is constant over a sliding window `[t-τ+1, t]`, the posterior is Gamma-distributed (analytical solution when prior is Gamma):

```
R(t) | data ~ Gamma(α + Σ I(t), β + Σ Λ(t))
```

Where `Λ(t) = Σ I(t-s) w(s)` is the total infectiousness.

This matches EpiEstim exactly and enables direct output parity.

---

## 3. EpiNow2 Analog Design

More complex: requires:
1. Delay distribution estimation (symptom onset → report)
2. Nowcasting (adjust right-truncated case counts)
3. Renewal equation Rt estimation on nowcast-corrected incidence
4. Stan/PyMC backend for full Bayesian inference

Primary backend: PyMC v5. Secondary: CmdStanPy (Stan) for exact EpiNow2 parity on NUTS sampler.

---

## 4. Key Dependencies

| Package | Role |
|---|---|
| `numpy`, `scipy` | Numerical core, distributions |
| `pandas` | Incidence time series |
| `pymc ≥ 5.16` | Bayesian Rt inference, compartmental calibration |
| `arviz ≥ 0.20` | MCMC diagnostics |
| `statsmodels` | CUSUM, time series components |
| `rpy2` (dev) | EpiEstim / EpiNow2 parity validation |
