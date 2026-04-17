"""episurv - Python library for infectious disease surveillance.

This library provides Python-native implementations of infectious disease
surveillance methods with parity to R packages EpiEstim, EpiNow2, and surveillance.
"""

from episurv.data import Incidence, SerialInterval
from episurv.rt import RtResult, estimate_rt_instant

__version__ = "0.1.0"
__all__ = [
    "__version__",
    "Incidence",
    "SerialInterval",
    "RtResult",
    "estimate_rt_instant",
]
