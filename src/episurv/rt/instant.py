"""Instantaneous Rt estimation via renewal equation (EpiEstim parity)."""

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from episurv.data.incidence import Incidence
from episurv.data.serial_interval import SerialInterval


@dataclass
class RtResult:
    """Result of instantaneous Rt estimation.

    Attributes:
        dates: DatetimeIndex of estimation dates (end of each window).
        r_mean: Posterior mean R(t) for each window.
        r_median: Posterior median R(t).
        r_lower: Lower CI bound (default 95%).
        r_upper: Upper CI bound (default 95%).
        window_size: Sliding window size in days.
        si_mean: Mean serial interval used.
        si_sd: SD of serial interval used.

    Example:
        >>> result = estimate_rt_instant(incidence, si, window_size=7)
        >>> print(result.r_mean)
    """

    dates: pd.DatetimeIndex
    r_mean: np.ndarray
    r_median: np.ndarray
    r_lower: np.ndarray
    r_upper: np.ndarray
    window_size: int
    ci_width: float
    si_mean: float
    si_sd: float

    def to_dataframe(self) -> pd.DataFrame:
        """Convert results to DataFrame."""
        return pd.DataFrame(
            {
                "date": self.dates,
                "r_mean": self.r_mean,
                "r_median": self.r_median,
                f"r_lower_{int(self.ci_width*100)}": self.r_lower,
                f"r_upper_{int(self.ci_width*100)}": self.r_upper,
            }
        )

    def __len__(self) -> int:
        """Number of Rt estimates."""
        return len(self.dates)


def estimate_rt_instant(
    incidence: Incidence,
    si: SerialInterval,
    window_size: int = 7,
    prior_mean: float = 5.0,
    prior_sd: float = 5.0,
    ci_width: float = 0.95,
    min_cases: int = 12,
) -> RtResult:
    """Estimate instantaneous Rt via renewal equation with Gamma posterior.

    Implements the EpiEstim method: assumes R(t) constant over sliding windows,
    computes analytical Gamma posterior given Gamma prior.

    Renewal equation:
        I(t) = R(t) × Σ_{s=1}^{t} I(t-s) × w(s)

    Where w(s) is the serial interval PMF.

    Args:
        incidence: Incidence time series.
        si: Serial interval distribution.
        window_size: Sliding window size in days (default 7).
        prior_mean: Prior mean for R (default 5.0).
        prior_sd: Prior SD for R (default 5.0).
        ci_width: Credible interval width (default 0.95).
        min_cases: Minimum total cases in window for estimation (default 12).

    Returns:
        RtResult with posterior summaries per window.

    Reference:
        Cori et al. (2013) A new framework and software to estimate time-varying
        reproduction numbers during epidemics. AJE.
    """
    # Input validations
    if not 0 < ci_width < 1:
        raise ValueError("ci_width must be between 0 and 1")
    if window_size <= 0:
        raise ValueError("window_size must be positive")
    if prior_mean <= 0:
        raise ValueError("prior_mean must be positive")
    if prior_sd <= 0:
        raise ValueError("prior_sd must be positive")

    n = incidence.n
    if window_size > n:
        raise ValueError(f"window_size ({window_size}) cannot exceed data length ({n})")

    values = incidence.values.values
    pmf = si.get_pmf(max_len=n)

    # Compute total infectiousness Λ(t) = Σ I(t-s) w(s)
    # This is a convolution: Λ = I * w (reversed)
    infectiousness = np.zeros(n)
    for t in range(1, n):
        # Sum over possible serial intervals
        max_s = min(t, len(pmf))
        for s in range(1, max_s + 1):
            infectiousness[t] += values[t - s] * pmf[s - 1]

    # Prior parameters (Gamma: shape=α, scale=β)
    # mean = αβ, var = αβ²
    prior_var = prior_sd**2
    prior_scale = prior_var / prior_mean
    prior_shape = prior_mean / prior_scale

    # Sliding window estimation
    results = []
    dates = []

    for t_end in range(window_size - 1, n):
        t_start = t_end - window_size + 1

        # Extract window
        window_cases = values[t_start : t_end + 1]
        window_lambda = infectiousness[t_start : t_end + 1]

        # Skip if insufficient cases
        total_cases = window_cases.sum()
        if total_cases < min_cases:
            results.append((np.nan, np.nan, np.nan, np.nan))
            dates.append(incidence.dates[t_end])
            continue

        # Posterior parameters (conjugate Gamma)
        # α_post = α_prior + Σ I(t)
        # β_post = β_prior + Σ Λ(t)
        post_shape = prior_shape + total_cases
        post_scale = 1.0 / (1.0 / prior_scale + window_lambda.sum())

        # Posterior summaries
        mean_r = post_shape * post_scale
        median_r = stats.gamma.ppf(0.5, post_shape, scale=post_scale)

        alpha = (1 - ci_width) / 2
        lower_r = stats.gamma.ppf(alpha, post_shape, scale=post_scale)
        upper_r = stats.gamma.ppf(1 - alpha, post_shape, scale=post_scale)

        results.append((mean_r, median_r, lower_r, upper_r))
        dates.append(incidence.dates[t_end])

    # Unpack results
    results_array = np.array(results)

    return RtResult(
        dates=pd.DatetimeIndex(dates),
        r_mean=results_array[:, 0],
        r_median=results_array[:, 1],
        r_lower=results_array[:, 2],
        r_upper=results_array[:, 3],
        window_size=window_size,
        ci_width=ci_width,
        si_mean=si.mean(),
        si_sd=si.sd(),
    )
