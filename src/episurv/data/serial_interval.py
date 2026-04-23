"""Serial interval distributions: parametric and empirical."""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd
from scipy import stats
from typing_extensions import Self


@dataclass
class SerialInterval:
    """Serial interval distribution for infectious disease transmission.

    Supports parametric distributions (Gamma, LogNormal) and empirical PMF.

    Attributes:
        dist_type: Distribution type ("gamma", "lognormal", or "empirical").
        params: Distribution parameters (for parametric) or None.
        pmf: Discrete probability mass function (for empirical or computed).
        max_si: Maximum serial interval days (truncation point).

    Example:
        >>> # Gamma distribution with mean=5.2, sd=1.5
        >>> si = SerialInterval.gamma(mean=5.2, sd=1.5)
        >>> # Or from empirical data
        >>> si = SerialInterval.empirical([2, 3, 3, 4, 5, 5, 5, 6, 7])
    """

    dist_type: Literal["gamma", "lognormal", "empirical"]
    params: dict[str, float] | None = None
    pmf: np.ndarray | None = None
    max_si: int = 30

    def __post_init__(self) -> None:
        """Compute PMF if parametric distribution specified."""
        if self.dist_type in ("gamma", "lognormal") and self.pmf is None:
            self.pmf = self._compute_pmf()

        pmf = self.pmf
        if pmf is not None:
            # Ensure PMF sums to 1
            pmf_norm = pmf / pmf.sum()
            self.pmf = pmf_norm
            # Ensure max_si matches PMF length
            self.max_si = len(pmf_norm)

    def _compute_pmf(self) -> np.ndarray:
        """Compute discrete PMF from parametric distribution."""
        if self.dist_type == "gamma":
            return self._gamma_pmf()
        elif self.dist_type == "lognormal":
            return self._lognormal_pmf()
        else:
            raise ValueError(f"Cannot compute PMF for {self.dist_type}")

    def _gamma_pmf(self) -> np.ndarray:
        """Compute Gamma PMF using EpiEstim's discr_si formula (shifted Gamma, shift=1).

        Matches EpiEstim::discr_si exactly:
            a = ((mu-1)/sigma)^2
            b = sigma^2 / (mu-1)
            w[k] = k*F(k,a,b) + (k-2)*F(k-2,a,b) - 2*(k-1)*F(k-1,a,b)
                   + a*b*(2*F(k-1,a+1,b) - F(k-2,a+1,b) - F(k,a+1,b))
            w[0] = 0 always.
        """
        if self.params is None:
            raise ValueError("Gamma params required")

        mu = self.params["mean"]
        sigma = self.params["sd"]

        if mu <= 1:
            raise ValueError("Gamma SI mean must be > 1 (EpiEstim shifted Gamma convention)")
        if sigma <= 0:
            raise ValueError("Gamma SI sd must be > 0")

        # EpiEstim parameterization (shape/scale of unshifted Gamma)
        a = ((mu - 1) / sigma) ** 2
        b = sigma**2 / (mu - 1)

        k = np.arange(0, self.max_si + 1, dtype=float)

        def cdf(x: np.ndarray, shape: float, scale: float) -> np.ndarray:
            result: np.ndarray = stats.gamma.cdf(x, shape, scale=scale)
            return result

        # EpiEstim formula
        w = (
            k * cdf(k, a, b)
            + (k - 2) * cdf(k - 2, a, b)
            - 2 * (k - 1) * cdf(k - 1, a, b)
            + a * b * (2 * cdf(k - 1, a + 1, b) - cdf(k - 2, a + 1, b) - cdf(k, a + 1, b))
        )
        # Clamp negatives (numerical noise)
        w = np.maximum(w, 0.0)
        # w[0] is always 0 by construction; enforce explicitly
        w[0] = 0.0

        # Normalize
        total = w.sum()
        if total == 0:
            raise ValueError("Gamma PMF sums to zero — check mean/sd parameters")
        result: np.ndarray = w / total
        return result

    def _lognormal_pmf(self) -> np.ndarray:
        """Compute LogNormal distribution PMF (naive discretization).

        Note: EpiEstim does not provide a discr_si equivalent for LogNormal;
        this uses P(SI=k) = F(k+0.5) - F(k-0.5) (continuity-corrected).
        w[0] = 0 is enforced.
        """
        if self.params is None:
            raise ValueError("LogNormal params required")

        mu = self.params["mu"]
        sigma = self.params["sigma"]

        k = np.arange(0, self.max_si + 1, dtype=float)
        # Continuity-corrected CDF differences
        cdf_vals: np.ndarray = stats.lognorm.cdf(k + 0.5, sigma, scale=np.exp(mu))
        w = np.diff(cdf_vals, prepend=stats.lognorm.cdf(0.5, sigma, scale=np.exp(mu)))
        w[0] = 0.0
        w = np.maximum(w, 0.0)

        total = w.sum()
        if total == 0:
            raise ValueError("LogNormal PMF sums to zero — check mean/sd parameters")
        result: np.ndarray = w / total
        return result

    @classmethod
    def gamma(cls, mean: float, sd: float, max_si: int = 30) -> Self:
        """Create Gamma-distributed serial interval using EpiEstim's discr_si convention.

        Uses a shifted Gamma (shift=1) matching EpiEstim::discr_si exactly.
        Requires mean > 1.

        Args:
            mean: Mean serial interval in days (must be > 1).
            sd: Standard deviation (must be > 0).
            max_si: Maximum serial interval length (truncation).

        Returns:
            SerialInterval with Gamma distribution (w[0]=0).
        """
        if mean <= 1:
            raise ValueError("mean must be > 1 for EpiEstim shifted Gamma (shift=1)")
        if sd <= 0:
            raise ValueError("sd must be > 0")

        return cls(
            dist_type="gamma",
            params={"mean": mean, "sd": sd},
            max_si=max_si,
        )

    @classmethod
    def lognormal(cls, mean: float, sd: float, max_si: int = 30) -> Self:
        """Create LogNormal-distributed serial interval.

        Args:
            mean: Mean serial interval in days (must be > 0).
            sd: Standard deviation (must be > 0).
            max_si: Maximum serial interval length (truncation).

        Returns:
            SerialInterval with LogNormal distribution (w[0]=0).
        """
        if mean <= 0:
            raise ValueError("mean must be > 0")
        if sd <= 0:
            raise ValueError("sd must be > 0")

        # Convert mean/sd to lognormal mu/sigma
        var = sd**2
        mu = np.log(mean**2 / np.sqrt(var + mean**2))
        sigma = np.sqrt(np.log(1 + var / mean**2))

        return cls(
            dist_type="lognormal",
            params={"mu": mu, "sigma": sigma, "mean": mean, "sd": sd},
            max_si=max_si,
        )

    @classmethod
    def empirical(cls, si_data: np.ndarray | list, max_si: int = 30) -> Self:
        """Create empirical serial interval from observed data.

        Args:
            si_data: Array of observed serial intervals (days).
            max_si: Maximum serial interval to consider.

        Returns:
            SerialInterval with empirical PMF.
        """
        si_array = np.asarray(si_data)

        # Warn if any values fall outside [0, max_si]
        n_out = int(np.sum((si_array < 0) | (si_array > max_si)))
        if n_out > 0:
            warnings.warn(
                f"{n_out} serial interval value(s) outside [0, {max_si}] were excluded.",
                UserWarning,
                stacklevel=2,
            )

        # Compute empirical PMF via histogram
        counts, _ = np.histogram(si_array, bins=np.arange(-0.5, max_si + 1.5))
        # Enforce w[0] = 0 (zero-day SI is epidemiologically impossible)
        counts[0] = 0
        total = counts.sum()
        if total == 0:
            raise ValueError("No serial interval data in valid range")
        pmf = counts / total

        return cls(
            dist_type="empirical",
            params={"mean": float(np.mean(si_array)), "sd": float(np.std(si_array))},
            pmf=pmf,
            max_si=max_si,
        )

    def get_pmf(self, max_len: int | None = None) -> np.ndarray:
        """Get PMF truncated to max_len.

        Args:
            max_len: Maximum length (defaults to max_si).

        Returns:
            Probability mass function array.
        """
        if self.pmf is None:
            raise ValueError("PMF not computed")

        max_len = self.max_si if max_len is None else max_len
        pmf = self.pmf[:max_len]

        # Renormalize if truncated
        if len(pmf) < len(self.pmf):
            pmf_sum = pmf.sum()
            if pmf_sum == 0:
                raise ValueError("Truncated PMF sums to zero")
            pmf = pmf / pmf_sum

        return pmf

    def mean(self) -> float:
        """Mean serial interval."""
        if self.params and "mean" in self.params:
            return self.params["mean"]
        if self.pmf is not None:
            return float(np.sum(np.arange(len(self.pmf)) * self.pmf))
        raise ValueError("Cannot compute mean")

    def sd(self) -> float:
        """Standard deviation of serial interval."""
        if self.params and "sd" in self.params:
            return self.params["sd"]
        if self.pmf is not None:
            mean = self.mean()
            variance = np.sum((np.arange(len(self.pmf)) - mean) ** 2 * self.pmf)
            return float(np.sqrt(variance))
        raise ValueError("Cannot compute sd")

    def to_dataframe(self) -> pd.DataFrame:
        """Convert to DataFrame with day/probability columns."""
        if self.pmf is None:
            raise ValueError("PMF not available")

        return pd.DataFrame({"day": np.arange(len(self.pmf)), "probability": self.pmf})
