"""Validation service for PyMaSC input validation.

This module provides the ValidationService that centralizes all
validation logic, removing it from handlers and ensuring consistent
validation across the application.

Key features:
- BAM file validation
- Configuration validation
- Input parameter validation
- Consistent error reporting
"""
import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Dict, Any, Tuple

from pysam import AlignmentFile

from PyMaSC.core.models import (
    CalculationConfig, MappabilityConfig, ExecutionConfig,
    CalculationTarget, ImplementationAlgorithm
)

logger = logging.getLogger(__name__)


@dataclass
class ValidationResult:
    """Result of validation operation.

    Attributes:
        is_valid: Whether validation passed
        errors: List of error messages
        warnings: List of warning messages
        metadata: Additional validation metadata
    """
    is_valid: bool
    errors: List[str]
    warnings: List[str]
    metadata: Optional[Dict[str, Any]] = None

    def add_error(self, message: str) -> None:
        """Add an error message."""
        self.errors.append(message)
        self.is_valid = False

    def add_warning(self, message: str) -> None:
        """Add a warning message."""
        self.warnings.append(message)

    def merge(self, other: 'ValidationResult') -> None:
        """Merge another validation result into this one."""
        self.errors.extend(other.errors)
        self.warnings.extend(other.warnings)
        self.is_valid = self.is_valid and other.is_valid

        # Merge metadata
        if other.metadata:
            if self.metadata is None:
                self.metadata = {}
            self.metadata.update(other.metadata)


class ValidationService(ABC):
    """Abstract base for validation services.

    Defines the interface for all validation operations,
    enabling different validation strategies.
    """

    @abstractmethod
    def validate_bam_file(self, path: str) -> ValidationResult:
        """Validate a BAM file.

        Args:
            path: Path to BAM file

        Returns:
            Validation result
        """
        pass

    @abstractmethod
    def validate_mappability_file(self, path: str) -> ValidationResult:
        """Validate a mappability file.

        Args:
            path: Path to mappability file

        Returns:
            Validation result
        """
        pass

    @abstractmethod
    def validate_calculation_config(self,
                                  config: CalculationConfig,
                                  bam_info: Optional[Dict[str, Any]] = None) -> ValidationResult:
        """Validate calculation configuration.

        Args:
            config: Calculation configuration
            bam_info: Optional BAM file information

        Returns:
            Validation result
        """
        pass

    @abstractmethod
    def validate_workflow_request(self,
                                bam_path: str,
                                output_prefix: str,
                                calculation_config: CalculationConfig,
                                mappability_config: Optional[MappabilityConfig] = None,
                                execution_config: Optional[ExecutionConfig] = None) -> ValidationResult:
        """Validate a complete workflow request.

        Args:
            bam_path: Path to BAM file
            output_prefix: Output file prefix
            calculation_config: Calculation configuration
            mappability_config: Optional mappability configuration
            execution_config: Optional execution configuration

        Returns:
            Validation result
        """
        pass


class StandardValidationService(ValidationService):
    """Standard implementation of validation service.

    Provides comprehensive validation for all PyMaSC inputs
    with detailed error messages and metadata.
    """

    def __init__(self, strict_mode: bool = False):
        """Initialize validation service.

        Args:
            strict_mode: Whether to use strict validation rules
        """
        self.strict_mode = strict_mode

    def validate_bam_file(self, path: str) -> ValidationResult:
        """Validate a BAM file for PyMaSC requirements."""
        result = ValidationResult(is_valid=True, errors=[], warnings=[])

        # Check file exists
        if not Path(path).exists():
            result.add_error(f"BAM file not found: {path}")
            return result

        try:
            with AlignmentFile(path) as bamfile:
                # Check if file is valid
                if bamfile.closed:
                    result.add_error("Cannot open BAM file")
                    return result

                # Check for references
                if not bamfile.references:
                    result.add_error("BAM file has no reference sequences")
                    return result

                # Check if sorted
                try:
                    header_dict = bamfile.header.to_dict()
                    if 'HD' in header_dict:
                        hd = header_dict['HD']
                        if 'SO' in hd and hd['SO'] != 'coordinate':
                            result.add_error("BAM file must be sorted by coordinate")
                    else:
                        result.add_warning("Cannot determine if BAM file is sorted")
                except (KeyError, AttributeError):
                    result.add_warning("Cannot determine if BAM file is sorted")

                # Check for index
                has_index = self._check_bam_index(bamfile)
                if not has_index:
                    result.add_warning("BAM file has no index - multiprocessing will be disabled")

                # Add metadata
                result.metadata = {
                    'n_references': bamfile.nreferences,
                    'references': list(bamfile.references),
                    'lengths': list(bamfile.lengths),
                    'has_index': has_index
                }

        except Exception as e:
            result.add_error(f"Error reading BAM file: {str(e)}")

        return result

    def _check_bam_index(self, bamfile: AlignmentFile) -> bool:
        """Check if BAM file has an index."""
        try:
            # Try to use index
            if bamfile.references:
                # Attempt a small fetch operation
                _ = list(bamfile.fetch(bamfile.references[0], 1, 2))
                return True
        except (OSError, ValueError):
            return False
        return False

    def validate_mappability_file(self, path: str) -> ValidationResult:
        """Validate a mappability file."""
        result = ValidationResult(is_valid=True, errors=[], warnings=[])

        # Check file exists
        if not Path(path).exists():
            result.add_error(f"Mappability file not found: {path}")
            return result

        # Check file extension
        path_obj = Path(path)
        if path_obj.suffix.lower() not in ['.bw', '.bigwig', '.json']:
            result.add_warning(f"Unexpected mappability file extension: {path_obj.suffix}")

        # Try to open and validate based on type
        if path_obj.suffix.lower() in ['.bw', '.bigwig']:
            result.merge(self._validate_bigwig_file(path))
        elif path_obj.suffix.lower() == '.json':
            result.merge(self._validate_json_mappability(path))

        return result

    def _validate_bigwig_file(self, path: str) -> ValidationResult:
        """Validate a BigWig mappability file."""
        result = ValidationResult(is_valid=True, errors=[], warnings=[])

        try:
            # Import BigWigReader to validate
            from PyMaSC.reader.bigwig import BigWigReader

            with BigWigReader(path) as reader:
                # Check if we can read chromosomes
                chroms = reader.chroms()
                if not chroms:
                    result.add_error("BigWig file has no chromosomes")

                result.metadata = {
                    'n_chromosomes': len(chroms),
                    'chromosomes': list(chroms.keys())
                }

        except ImportError:
            result.add_error("pyBigWig not installed - cannot validate BigWig files")
        except Exception as e:
            result.add_error(f"Error reading BigWig file: {str(e)}")

        return result

    def _validate_json_mappability(self, path: str) -> ValidationResult:
        """Validate a JSON mappability file."""
        result = ValidationResult(is_valid=True, errors=[], warnings=[])

        try:
            import json
            with open(path) as f:
                data = json.load(f)

            # Check expected structure
            if not isinstance(data, dict):
                result.add_error("JSON mappability file must contain a dictionary")
                return result

            # Check for expected keys
            if 'stat' not in data:
                result.add_error("JSON mappability file missing 'stat' key")

            if 'chromosomes' in data:
                result.metadata = {
                    'n_chromosomes': len(data['chromosomes']),
                    'chromosomes': list(data['chromosomes'].keys())
                }

        except json.JSONDecodeError as e:
            result.add_error(f"Invalid JSON format: {str(e)}")
        except Exception as e:
            result.add_error(f"Error reading JSON file: {str(e)}")

        return result

    def validate_calculation_config(self,
                                  config: CalculationConfig,
                                  bam_info: Optional[Dict[str, Any]] = None) -> ValidationResult:
        """Validate calculation configuration."""
        result = ValidationResult(is_valid=True, errors=[], warnings=[])

        # Validate target and implementation
        if config.target not in CalculationTarget:
            result.add_error(f"Invalid calculation target: {config.target}")
        if config.implementation not in ImplementationAlgorithm:
            result.add_error(f"Invalid implementation algorithm: {config.implementation}")

        # Validate max_shift
        if config.max_shift <= 0:
            result.add_error("max_shift must be positive")
        elif config.max_shift > 10000:
            result.add_warning("Very large max_shift may impact performance")

        # Validate mapq_criteria
        if config.mapq_criteria < 0:
            result.add_error("mapq_criteria cannot be negative")
        elif config.mapq_criteria > 60:
            result.add_warning("Very high mapq_criteria may filter out too many reads")

        # Validate read_length if provided
        if config.read_length is not None:
            if config.read_length <= 0:
                result.add_error("read_length must be positive")
            elif config.read_length > 1000:
                result.add_warning("Unusually long read_length specified")

        # Validate references against BAM info if provided
        if bam_info and 'references' in bam_info:
            bam_refs = set(bam_info['references'])

            if config.references:
                for ref in config.references:
                    if ref not in bam_refs:
                        result.add_error(f"Reference '{ref}' not found in BAM file")

        # Target-specific validation
        if config.target == CalculationTarget.MSCC and config.skip_ncc:
            result.add_warning("MSCC target with skip_ncc may produce limited results")

        return result

    def validate_workflow_request(self,
                                bam_path: str,
                                output_prefix: str,
                                calculation_config: CalculationConfig,
                                mappability_config: Optional[MappabilityConfig] = None,
                                execution_config: Optional[ExecutionConfig] = None) -> ValidationResult:
        """Validate a complete workflow request."""
        result = ValidationResult(is_valid=True, errors=[], warnings=[])

        # Validate BAM file
        bam_result = self.validate_bam_file(bam_path)
        result.merge(bam_result)

        # Validate output path
        if not output_prefix:
            result.add_error("Output prefix cannot be empty")
        else:
            output_dir = Path(output_prefix).parent
            if output_dir and not output_dir.exists():
                result.add_warning(f"Output directory does not exist: {output_dir}")

        # Validate calculation config with BAM info
        config_result = self.validate_calculation_config(
            calculation_config,
            bam_result.metadata
        )
        result.merge(config_result)

        # Validate mappability if needed
        if calculation_config.target in [CalculationTarget.MSCC, CalculationTarget.BOTH]:
            if mappability_config is None:
                result.add_error(f"{calculation_config.target.value} target requires mappability configuration")
            else:
                if mappability_config.mappability_path is not None:
                    map_result = self.validate_mappability_file(str(mappability_config.mappability_path))
                else:
                    map_result = ValidationResult(is_valid=False, errors=[], warnings=[])
                    map_result.add_error("Mappability path is required but not provided")
                result.merge(map_result)

                # Check chromosome compatibility
                if (bam_result.metadata and map_result.metadata and
                    'chromosomes' in bam_result.metadata and 'chromosomes' in map_result.metadata):

                    bam_chroms = set(bam_result.metadata['chromosomes'])
                    map_chroms = set(map_result.metadata['chromosomes'])

                    missing_in_map = bam_chroms - map_chroms
                    if missing_in_map:
                        result.add_warning(
                            f"Chromosomes in BAM but not in mappability: {missing_in_map}"
                        )

        # Validate execution config if provided
        if execution_config:
            if (execution_config.worker_count > 1 and 
                bam_result.metadata is not None and 
                not bam_result.metadata.get('has_index', False)):
                result.add_warning("Multiprocessing requested but BAM file has no index")

        return result


# Factory function for service creation
def create_validation_service(strict: bool = False) -> ValidationService:
    """Create validation service.

    Args:
        strict: Whether to use strict validation

    Returns:
        Validation service instance
    """
    return StandardValidationService(strict_mode=strict)