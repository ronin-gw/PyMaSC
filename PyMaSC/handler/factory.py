"""Factory pattern implementations for PyMaSC calculators and workers.

Provides factory functions for creating calculator and worker instances
with appropriate configurations.

Key components:
- CompositeCalculator: Combines NCC and MSCC calculations
- create_calculator: Factory function for calculator creation
- create_worker: Factory function for worker creation
"""
import logging
from typing import Optional, Any, Tuple, cast
from multiprocessing.synchronize import Lock

from PyMaSC.result import BothChromResult, BothGenomeWideResult

from PyMaSC.interfaces.config import (
    PyMaSCConfig, Algorithm, CalculationTarget, mappability_configured
)
from PyMaSC.interfaces.calculator import (
    CrossCorrelationCalculator,
    NCCCalculatorModel, MSCCCalculatorModel, BothCalculatorModel
)

from PyMaSC.core.successive.ncc import NaiveCCCalculator
from PyMaSC.core.successive.mscc import MSCCCalculator
from PyMaSC.core.bitarray.mscc import CCBitArrayCalculator

from PyMaSC.reader.bigwig import BigWigReader
from PyMaSC.core.mappability import BWFeederWithMappableRegionSum

logger = logging.getLogger(__name__)


class CompositeCalculator(BothCalculatorModel):
    """Composite calculator that runs both NCC and MSCC calculations.

    Used when SUCCESSIVE algorithm is combined with mappability data.
    """

    _ncc_calculator: NCCCalculatorModel
    _mscc_calculator: MSCCCalculatorModel

    def __init__(self, ncc_calculator: NCCCalculatorModel, mscc_calculator: MSCCCalculatorModel):
        """Initialize composite calculator.

        Args:
            ncc_calculator: NCC calculator instance
            mscc_calculator: MSCC calculator instance
        """
        self._ncc_calculator = ncc_calculator
        self._mscc_calculator = mscc_calculator

    def feed_forward_read(self, chrom: str, pos: int, readlen: int) -> None:
        """Process a forward strand read in both calculators."""
        self._ncc_calculator.feed_forward_read(chrom, pos, readlen)
        self._mscc_calculator.feed_forward_read(chrom, pos, readlen)

    def feed_reverse_read(self, chrom: str, pos: int, readlen: int) -> None:
        """Process a reverse strand read in both calculators."""
        self._ncc_calculator.feed_reverse_read(chrom, pos, readlen)
        self._mscc_calculator.feed_reverse_read(chrom, pos, readlen)

    def finishup_calculation(self) -> None:
        """Complete calculation in both calculators."""
        self._ncc_calculator.finishup_calculation()
        self._mscc_calculator.finishup_calculation()

    def flush(self, chrom: Optional[str] = None) -> None:
        """Flush both calculators."""
        self._ncc_calculator.flush(chrom)
        self._mscc_calculator.flush(chrom)

    def get_result(self, chrom: str) -> BothChromResult:
        """Get combined results for a specific chromosome.

        Args:
            chrom: Chromosome name

        Returns:
            BothChromResult containing both NCC and MSCC results
        """
        return BothChromResult(
            chrom=self._ncc_calculator.get_result(chrom),
            mappable_chrom=self._mscc_calculator.get_result(chrom)
        )

    def get_whole_result(self) -> BothGenomeWideResult:
        """Get combined genome-wide results from both calculators.

        Returns:
            BothGenomeWideResult containing aggregated NCC and MSCC data
        """
        nccresult = self._ncc_calculator.get_whole_result()
        msccresult = self._mscc_calculator.get_whole_result()

        return BothGenomeWideResult(
            genomelen=nccresult.genomelen,
            forward_sum=nccresult.forward_sum,
            reverse_sum=nccresult.reverse_sum,
            chroms=nccresult.chroms.copy(),
            mappable_chroms=msccresult.chroms.copy(),
            forward_read_len_sum=msccresult.forward_read_len_sum,
            reverse_read_len_sum=msccresult.reverse_read_len_sum
        )


def create_calculator(
    config: PyMaSCConfig,
    logger_lock: Optional[Lock] = None,
    progress_hook: Optional[Any] = None
) -> Tuple[CrossCorrelationCalculator, Optional[BigWigReader]]:
    """Create a cross-correlation calculator based on configuration.

    Args:
        config: PyMaSC configuration containing implementation and target settings
        logger_lock: Optional lock for thread-safe logging
        progress_hook: Optional progress reporting hook for BITARRAY implementation

    Returns:
        Tuple of (calculator instance, optional BigWig reader)

    Raises:
        ValueError: If unsupported implementation type is specified
    """
    # Create calculator based on implementation and target
    if config.implementation is Algorithm.SUCCESSIVE:
        return _create_successive_calculator(
            config, logger_lock
        )
    elif config.implementation is Algorithm.BITARRAY:
        return _create_bitarray_calculator(
            config, logger_lock, progress_hook
        )
    else:
        raise ValueError(f"Unsupported implementation: {config.implementation.value}")


def _create_successive_calculator(
    config: PyMaSCConfig,
    logger_lock: Optional[Lock]
) -> Tuple[CrossCorrelationCalculator, Optional[BigWigReader]]:
    """Create calculator using SUCCESSIVE algorithm.

    Args:
        config: Configuration for calculator setup
        logger_lock: Optional lock for thread-safe logging

    Returns:
        Tuple of (SUCCESSIVE calculator, optional BigWig reader)

    Raises:
        ValueError: If unsupported calculation target is specified
    """
    if config.target is CalculationTarget.NCC:
        # Create NCC calculator
        return NaiveCCCalculator(
            config.max_shift,
            config.references,
            config.lengths,
            logger_lock
        ), None

    elif config.target is CalculationTarget.MSCC:
        # Create MSCC calculator
        return _create_mscc_calculator(
            config, logger_lock
        )

    elif config.target is CalculationTarget.BOTH:
        # Create composite calculator for both NCC and MSCC
        ncc_calc = NaiveCCCalculator(
            config.max_shift,
            config.references,
            config.lengths,
            logger_lock
        )

        mscc_calc, bwfeeder = _create_mscc_calculator(config, logger_lock)

        return CompositeCalculator(ncc_calc, mscc_calc), bwfeeder

    else:
        raise ValueError(f"Unsupported target: {config.target.value}")


def _create_mscc_calculator(
    config: PyMaSCConfig,
    logger_lock: Optional[Lock]
) -> Tuple[MSCCCalculator, BWFeederWithMappableRegionSum]:
    """Create pure MSCC calculator with mappability correction.

    Args:
        config: Configuration with mappability_path and read_length set
        logger_lock: Optional lock for thread-safe logging

    Returns:
        Tuple of (MSCC calculator, BigWig feeder for mappability data)

    Raises:
        AssertionError: If required mappability configuration is missing
    """

    assert config.mappability_path is not None, "Mappability path must be set for MSCC calculation"
    assert config.read_length is not None, "Read length must be set for MSCC calculation"

    bwfeeder = BWFeederWithMappableRegionSum(
        str(config.mappability_path),
        config.max_shift
    )

    calculator = MSCCCalculator(
        config.max_shift,
        config.read_length,
        list(config.references),
        config.lengths,
        bwfeeder,
        logger_lock
    )

    return calculator, bwfeeder


def _create_bitarray_calculator(
    config: PyMaSCConfig,
    logger_lock: Optional[Lock],
    progress_hook: Optional[Any]
) -> Tuple[CrossCorrelationCalculator, Optional[BigWigReader]]:
    """Create calculator using BITARRAY algorithm.

    Args:
        config: Configuration for calculator setup
        logger_lock: Optional lock for thread-safe logging
        progress_hook: Optional progress reporting hook

    Returns:
        Tuple of (BITARRAY calculator, optional BigWig reader)
    """

    bwfeeder = None
    if mappability_configured(config):
        bwfeeder = BigWigReader(config.mappability_path)

    config = cast(PyMaSCConfig, config)
    assert config.read_length is not None, "Read length must be set for BitArray calculator"

    calculator = CCBitArrayCalculator(
        config.max_shift,
        config.read_length,
        list(config.references),
        config.lengths,
        bwfeeder,
        config.skip_ncc,
        logger_lock,
        progress_hook
    )

    return calculator, bwfeeder
