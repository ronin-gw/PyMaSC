"""Type stubs for readlen module."""
from typing import Literal


def estimate_readlen(path: str, esttype: Literal['MEAN', 'MEDIAN', 'MODE', 'MIN', 'MAX'], mapq_criteria: int) -> int:
    """Estimate read length from BAM file using specified method.

    Args:
        path: Path to BAM file
        esttype: Estimation method ('MEAN', 'MEDIAN', 'MODE', 'MIN', 'MAX')
        mapq_criteria: Minimum mapping quality threshold

    Returns:
        Estimated read length in base pairs

    Raises:
        ValueError: If estimation method is invalid or insufficient data
    """
    ...
