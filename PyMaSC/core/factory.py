"""Factory pattern implementations for PyMaSC calculators and workers.

This module provides factory classes that encapsulate the creation logic
for various calculator and worker implementations. The factories handle
the complexity of different constructor signatures and dependencies,
providing a unified interface for object creation.

Key factories:
- CalculatorFactory: Creates cross-correlation calculator instances
- WorkerFactory: Creates worker instances for multiprocessing
- CalculatorAdapter: Adapts existing calculators to new interface
"""
import logging
from typing import Optional, Any
from multiprocessing.synchronize import Lock

from PyMaSC.core.result import BothChromResult, BothGenomeWideResult

from .interfaces.calculator import CrossCorrelationCalculator, BothCalculatorModel
from .models import (
    CalculationConfig, MappabilityConfig,
    CalculationTarget, ImplementationAlgorithm
)

# Import existing calculator implementations
# These will be wrapped to conform to the new interface
from PyMaSC.core.ncc import NaiveCCCalculator as _NativeCCCalculator
from PyMaSC.core.mscc import MSCCCalculator as _NativeMSCCCalculator
from PyMaSC.bacore.mscc import CCBitArrayCalculator as _NativeBitArrayCalculator

# Import BigWig reader for mappability
from PyMaSC.reader.bigwig import BigWigReader
from PyMaSC.core.mappability import BWFeederWithMappableRegionSum

logger = logging.getLogger(__name__)


class CompositeCalculator(BothCalculatorModel):
    """Composite calculator that runs both NCC and MSCC calculations.

    This calculator is used when SUCCESSIVE algorithm is combined with
    mappability data, matching the behavior of the original implementation.
    """

    _ncc_calculator: _NativeCCCalculator
    _mscc_calculator: _NativeMSCCCalculator

    def __init__(self, ncc_calculator: Any, mscc_calculator: Any):
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
        return BothChromResult(
            chrom=self._ncc_calculator.get_result(chrom),
            mappable_chrom=self._mscc_calculator.get_result(chrom)
        )

    def get_whole_result(self) -> BothGenomeWideResult:
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


class CalculatorFactory:
    """Factory for creating cross-correlation calculator instances.

    This factory encapsulates the complex logic of creating different
    calculator types with their varying constructor signatures and
    dependencies. It provides a unified interface for calculator creation
    regardless of the underlying implementation.
    """

    @staticmethod
    def create_calculator(
        target: CalculationTarget,
        implementation: ImplementationAlgorithm,
        config: CalculationConfig,
        mappability_config: Optional[MappabilityConfig] = None,
        logger_lock: Optional[Lock] = None,
        progress_hook: Optional[Any] = None
    ) -> CrossCorrelationCalculator:
        """Create a calculator instance based on target and implementation.

        Args:
            target: What type of cross-correlation to calculate
            implementation: How to implement the calculation
            config: Calculation configuration
            mappability_config: Optional mappability configuration
            logger_lock: Optional thread lock for logging
            progress_hook: Optional progress reporting hook

        Returns:
            Calculator instance wrapped in adapter

        Raises:
            ValueError: If combination is not supported or configuration is invalid
        """
        # Validate mappability requirements
        if target in [CalculationTarget.MSCC, CalculationTarget.BOTH]:
            if not mappability_config or not mappability_config.is_enabled():
                raise ValueError(
                    f"Target {target.value} requires mappability configuration"
                )

        # Create calculator based on implementation and target
        if implementation == ImplementationAlgorithm.SUCCESSIVE:
            return CalculatorFactory._create_successive_calculator(
                target, config, mappability_config, logger_lock
            )
        elif implementation == ImplementationAlgorithm.BITARRAY:
            return CalculatorFactory._create_bitarray_calculator(
                target, config, mappability_config, logger_lock, progress_hook
            )
        else:
            raise ValueError(f"Unsupported implementation: {implementation.value}")

    @staticmethod
    def _create_successive_calculator(
        target: CalculationTarget,
        config: CalculationConfig,
        mappability_config: Optional[MappabilityConfig],
        logger_lock: Optional[Lock]
    ) -> CrossCorrelationCalculator:
        """Create calculator using SUCCESSIVE algorithm."""
        if target == CalculationTarget.NCC:
            # Create NCC calculator
            return _NativeCCCalculator(
                config.max_shift,
                config.references,
                config.lengths,
                logger_lock
            )

        elif target == CalculationTarget.MSCC:
            # Create MSCC calculator
            return CalculatorFactory._create_mscc_calculator(
                config, mappability_config, logger_lock
            )

        elif target == CalculationTarget.BOTH:
            # Create composite calculator for both NCC and MSCC
            ncc_calc = _NativeCCCalculator(
                config.max_shift,
                config.references,
                config.lengths,
                logger_lock
            )

            # Create BWFeeder for MSCC
            bwfeeder = BWFeederWithMappableRegionSum(
                str(mappability_config.mappability_path),
                config.max_shift
            )

            read_len = config.read_length or mappability_config.read_len
            if read_len is None:
                raise ValueError("Read length must be specified for MSCC")

            mscc_calc = _NativeMSCCCalculator(
                config.max_shift,
                read_len,
                config.references,
                config.lengths,
                bwfeeder,
                logger_lock
            )

            return CompositeCalculator(ncc_calc, mscc_calc)

    @staticmethod
    def _create_bitarray_calculator(
        target: CalculationTarget,
        config: CalculationConfig,
        mappability_config: Optional[MappabilityConfig],
        logger_lock: Optional[Lock],
        progress_hook: Optional[Any]
    ) -> CrossCorrelationCalculator:
        """Create calculator using BITARRAY algorithm."""
        # BitArray handles all targets in a single calculator
        bwfeeder = None
        if mappability_config and mappability_config.is_enabled():
            assert mappability_config.mappability_path is not None
            bwfeeder = BigWigReader(mappability_config.mappability_path)

        read_len = config.read_length
        if mappability_config and mappability_config.read_len:
            read_len = mappability_config.read_len

        if read_len is None:
            raise ValueError("Read length must be specified for BitArray")

        # Determine skip_ncc based on target
        skip_ncc = target == CalculationTarget.MSCC or config.skip_ncc

        return _NativeBitArrayCalculator(
            config.max_shift,
            read_len,
            config.references,
            config.lengths,
            bwfeeder,
            skip_ncc,
            logger_lock,
            progress_hook
        )

    @staticmethod
    def _create_mscc_calculator(
        config: CalculationConfig,
        mappability_config: MappabilityConfig,
        logger_lock: Optional[Lock]
    ) -> CrossCorrelationCalculator:
        """Create pure MSCC calculator."""
        bwfeeder = BWFeederWithMappableRegionSum(
            str(mappability_config.mappability_path),
            config.max_shift
        )

        read_len = config.read_length or mappability_config.read_len
        if read_len is None:
            raise ValueError("Read length must be specified for MSCC")

        return _NativeMSCCCalculator(
            config.max_shift,
            read_len,
            config.references,
            config.lengths,
            bwfeeder,
            logger_lock
        )
