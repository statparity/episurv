"""Basic tests for episurv package."""

import episurv


def test_version() -> None:
    """Test that version is accessible."""
    assert episurv.__version__ == "0.1.0"
