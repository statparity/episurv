"""Serial interval distributions: parametric and empirical."""

from __future__ import annotations

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
        """Compute Gamma distribution PMF (discretized)."""
        if self.params is None:
            raise ValueError("Gamma params required")

        shape = self.params["shape"]
        scale = self.params["scale"]

        # Discretized Gamma: P(SI = k) = F(k+1) - F(k)
        k_vals = np.arange(0, self.max_si + 1)
        cdf_vals: np.ndarray = stats.gamma.cdf(k_vals + 1, shape, scale=scale)
        pmf = np.diff(cdf_vals, prepend=0)

        # Truncate and normalize
        pmf_truncated: np.ndarray = pmf[: self.max_si]
        result: np.ndarray = pmf_truncated / pmf_truncated.sum()
        return result

    def _lognormal_pmf(self) -> np.ndarray:
        """Compute LogNormal distribution PMF (discretized)."""
        if self.params is None:
            raise ValueError("LogNormal params required")

        mu = self.params["mu"]
        sigma = self.params["sigma"]

        # Discretized LogNormal
        k_vals = np.arange(0, self.max_si + 1)
        cdf_vals: np.ndarray = stats.lognorm.cdf(k_vals + 1, sigma, scale=np.exp(mu))
        pmf = np.diff(cdf_vals, prepend=0)

        pmf_truncated: np.ndarray = pmf[: self.max_si]
        result: np.ndarray = pmf_truncated / pmf_truncated.sum()
        return result

    @classmethod
    def gamma(cls, mean: float, sd: float, max_si: int = 30) -> Self:
        """Create Gamma-distributed serial interval.

        Args:
            mean: Mean serial interval in days.
            sd: Standard deviation.
            max_si: Maximum serial interval length (truncation).

        Returns:
            SerialInterval with Gamma distribution.
        """
        # Gamma parameterization: shape = mean^2 / var, scale = var / mean
        var = sd**2
        shape = mean**2 / var
        scale = var / mean

        return cls(
            dist_type="gamma",
            params={"shape": shape, "scale": scale, "mean": mean, "sd": sd},
            max_si=max_si,
        )

    @classmethod
    def lognormal(cls, mean: float, sd: float, max_si: int = 30) -> Self:
        """Create LogNormal-distributed serial interval.

        Args:
            mean: Mean serial interval in days.
            sd: Standard deviation.
            max_si: Maximum serial interval length (truncation).

        Returns:
            SerialInterval with LogNormal distribution.
        """
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

        # Compute empirical PMF via histogram
        counts, _ = np.histogram(si_array, bins=np.arange(-0.5, max_si + 1.5))
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

        max_len = max_len or self.max_si
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
