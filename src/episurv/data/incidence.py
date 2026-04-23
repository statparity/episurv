"""Incidence data object with date indexing and right-censoring support."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd
from typing_extensions import Self


@dataclass
class Incidence:
    """Time-indexed incidence data with optional right-censoring.

    Attributes:
        dates: DatetimeIndex of observation dates.
        values: Array-like of case counts.
        is_censored: Boolean mask for right-censored values (default all False).

    Example:
        >>> import pandas as pd
        >>> dates = pd.date_range("2020-03-01", periods=10)
        >>> inc = Incidence(dates, [1, 2, 5, 3, 8, 10, 12, 15, 18, 20])
    """

    dates: pd.DatetimeIndex
    values: pd.Series
    is_censored: pd.Series | None = None

    def __post_init__(self) -> None:
        """Validate inputs and align types."""
        # Ensure dates is DatetimeIndex
        if not isinstance(self.dates, pd.DatetimeIndex):
            self.dates = pd.DatetimeIndex(self.dates)

        # Ensure values is Series with same index as dates
        if not isinstance(self.values, pd.Series):
            self.values = pd.Series(self.values, index=self.dates)
        else:
            self.values = self.values.copy()
            self.values.index = self.dates

        # Default is_censored to all False
        if self.is_censored is None:
            censored = pd.Series(False, index=self.dates)
        elif not isinstance(self.is_censored, pd.Series):
            censored = pd.Series(self.is_censored, index=self.dates)
        else:
            censored = self.is_censored.copy()
            censored.index = self.dates
        self.is_censored = censored

        # Validate lengths match
        if not (len(self.dates) == len(self.values) == len(self.is_censored)):
            raise ValueError("dates, values, and is_censored must have same length")

        # Validate values are non-negative
        if (self.values < 0).any():
            raise ValueError("incidence values must be non-negative")

    @property
    def n(self) -> int:
        """Number of time points."""
        return len(self.dates)

    @property
    def total(self) -> int:
        """Total case count (including censored)."""
        return int(self.values.sum())

    @property
    def uncensored_values(self) -> pd.Series:
        """Values with censored entries set to NaN."""
        censored = self.is_censored
        assert censored is not None  # Type narrowing
        result: pd.Series = self.values.where(~censored)
        return result

    def window(self, start_idx: int, end_idx: int) -> Self:
        """Extract a time window [start_idx, end_idx] (inclusive)."""
        censored = self.is_censored
        assert censored is not None  # Type narrowing
        return self.__class__(
            dates=self.dates[start_idx : end_idx + 1],
            values=self.values.iloc[start_idx : end_idx + 1],
            is_censored=censored.iloc[start_idx : end_idx + 1],
        )

    def to_dataframe(self) -> pd.DataFrame:
        """Convert to DataFrame with dates as index."""
        return pd.DataFrame(
            {"incidence": self.values, "is_censored": self.is_censored},
            index=self.dates,
        )

    @classmethod
    def from_dataframe(
        cls,
        df: pd.DataFrame,
        date_col: str | None = None,
        value_col: str = "incidence",
        censored_col: str | None = None,
    ) -> Self:
        """Create Incidence from DataFrame.

        Args:
            df: Input DataFrame.
            date_col: Column with dates (default uses index).
            value_col: Column with incidence values.
            censored_col: Column with censoring flags (optional).
        """
        if date_col is not None:
            if date_col not in df.columns:
                raise KeyError(f"Column '{date_col}' not found in DataFrame")
            dates = pd.DatetimeIndex(df[date_col])
        else:
            dates = pd.DatetimeIndex(df.index)

        if value_col not in df.columns:
            raise KeyError(f"Column '{value_col}' not found in DataFrame")
        values = df[value_col]

        is_censored = None
        if censored_col is not None:
            if censored_col not in df.columns:
                raise KeyError(f"Column '{censored_col}' not found in DataFrame")
            is_censored = df[censored_col]

        return cls(dates=dates, values=values, is_censored=is_censored)
