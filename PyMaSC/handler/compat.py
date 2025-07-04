"""Compatibility layer for legacy handler usage.

This module provides wrapper classes that maintain backward compatibility
with the existing CCCalcHandler and BACalcHandler interfaces while
delegating to the new UnifiedCalcHandler implementation.

These wrappers ensure that existing code continues to work without
modification during the transition to the unified architecture.
"""
import logging
from typing import Optional, Any

from .unified import UnifiedCalcHandler
from PyMaSC.core.models import (
    CalculationConfig, MappabilityConfig, ExecutionConfig,
    AlgorithmType, ExecutionMode
)

logger = logging.getLogger(__name__)


class CCCalcHandler(UnifiedCalcHandler):
    """Compatibility wrapper for CCCalcHandler.
    
    This class maintains the original CCCalcHandler interface while
    using UnifiedCalcHandler internally. It converts the legacy
    parameter format to the new configuration objects.
    """
    
    def __init__(self, path: str, esttype: str, max_shift: int, 
                 mapq_criteria: int, nworker: int = 1, 
                 skip_ncc: bool = False, chromfilter: Optional[Any] = None):
        """Initialize with legacy CCCalcHandler parameters.
        
        Args:
            path: Path to input BAM file
            esttype: Read length estimation method
            max_shift: Maximum shift distance
            mapq_criteria: Minimum mapping quality
            nworker: Number of worker processes
            skip_ncc: Whether to skip NCC calculation
            chromfilter: Chromosome filtering patterns
        """
        # Create configuration objects from legacy parameters
        calc_config = CalculationConfig(
            algorithm=AlgorithmType.SUCCESSIVE,  # CCCalcHandler uses successive
            max_shift=max_shift,
            mapq_criteria=mapq_criteria,
            skip_ncc=skip_ncc
        )
        
        # Store esttype as attribute for backward compatibility
        calc_config.esttype = esttype
        calc_config.chromfilter = chromfilter
        
        exec_config = ExecutionConfig(
            mode=ExecutionMode.MULTI_PROCESS if nworker > 1 else ExecutionMode.SINGLE_PROCESS,
            worker_count=nworker
        )
        
        # Initialize unified handler
        super().__init__(path, calc_config, exec_config)


class BACalcHandler(UnifiedCalcHandler):
    """Compatibility wrapper for BACalcHandler.
    
    This class maintains the original BACalcHandler interface while
    using UnifiedCalcHandler internally. It converts the legacy
    parameter format to the new configuration objects.
    """
    
    def __init__(self, path: str, esttype: str, max_shift: int, 
                 mapq_criteria: int, nworker: int = 1, 
                 skip_ncc: bool = False, chromfilter: Optional[Any] = None):
        """Initialize with legacy BACalcHandler parameters.
        
        Args:
            path: Path to input BAM file
            esttype: Read length estimation method
            max_shift: Maximum shift distance
            mapq_criteria: Minimum mapping quality
            nworker: Number of worker processes
            skip_ncc: Whether to skip NCC calculation
            chromfilter: Chromosome filtering patterns
        """
        # Create configuration objects from legacy parameters
        calc_config = CalculationConfig(
            algorithm=AlgorithmType.BITARRAY,  # BACalcHandler uses BitArray
            max_shift=max_shift,
            mapq_criteria=mapq_criteria,
            skip_ncc=skip_ncc
        )
        
        # Store esttype as attribute for backward compatibility
        calc_config.esttype = esttype
        calc_config.chromfilter = chromfilter
        
        exec_config = ExecutionConfig(
            mode=ExecutionMode.MULTI_PROCESS if nworker > 1 else ExecutionMode.SINGLE_PROCESS,
            worker_count=nworker
        )
        
        # Initialize unified handler
        super().__init__(path, calc_config, exec_config)