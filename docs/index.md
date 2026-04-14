# episurv

Python library for infectious disease surveillance — EpiEstim / EpiNow2 / surveillance parity.

## Installation

```bash
pip install episurv
```

## Quick Start

```python
import episurv
```

## Features

- **Rt Estimation**: Instantaneous reproduction number via renewal equation (EpiEstim parity)
- **Nowcasting**: Right-truncation adjustment for reporting delays (EpiNow2 parity)
- **Outbreak Detection**: CUSUM, EARS, and scan statistics (surveillance parity)
- **R Validation**: Automated parity testing against reference R implementations
