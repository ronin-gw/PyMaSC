"""Statistics file reading utilities.

Reads and parses PyMaSC statistics files containing analysis results.
"""
import logging
from typing import Dict, Union, Literal
from pathlib import Path

from PyMaSC.utils.output import catch_IOError
from PyMaSC.interfaces.output import (
    StatLabels,
    SummaryItems,
    CorrelationStats, CorrelationItems, NCC_LABELS, MSCC_LABELS, SUMMARY_LABELS,
    OutputStats, nanint, nanfloat
)

logger = logging.getLogger(__name__)


def _parse_stats_section(
    raw_data: Dict[str, str],
    labels: StatLabels,
    dataclass_type: type,
    section_name: str = ""
) -> Dict[str, Union[str, int, float]]:
    """Parse a section of stats data using label mappings.

    Args:
        raw_data: Raw data from the stats file
        labels: StatLabels instance for this section
        dataclass_type: The dataclass type to parse into (e.g., SummaryItems, CorrelationItems)
        section_name: Name of the section for error messages (e.g., "NCC", "MSCC")

    Returns:
        Dictionary of parsed field values

    Raises:
        ValueError: If required fields are missing
    """
    # Create label-to-field mapping
    mapping = {}
    for field_name, label_value in labels.__dict__.items():
        mapping[label_value] = field_name

    # Parse data
    parsed_data = {}
    for label, field_name in mapping.items():
        if label in raw_data:
            field_type = _get_field_type(field_name, dataclass_type)
            parsed_data[field_name] = _parse_value(raw_data[label], field_name, field_type)
        else:
            # Handle missing fields
            if _is_nanfield(field_name, dataclass_type):
                prefix = f"{section_name} " if section_name else ""
                logger.warning(f"Missing {prefix}field '{field_name}', setting to 'nan'")
                parsed_data[field_name] = 'nan'
            else:
                prefix = f"{section_name} " if section_name else ""
                raise ValueError(f"Required {prefix}field '{field_name}' is missing from statistics file")

    return parsed_data


def _parse_value(value_str: str, field_name: str, target_type: type) -> Union[str, int, float]:
    """Parse string value to appropriate type with nan handling."""
    if value_str == 'nan':
        return 'nan'

    try:
        return target_type(value_str)
    except (ValueError, TypeError) as e:
        logger.warning(f"Failed to parse value '{value_str}' for field '{field_name}': {e}")
        return 'nan'


def _get_field_type(field_name: str, dataclass_type: type) -> type:
    """Get the expected type for a field from dataclass annotations."""
    annotations = dataclass_type.__annotations__
    if field_name in annotations:
        annotation = annotations[field_name]
        # Handle Union types like nanint, nanfloat
        if hasattr(annotation, '__origin__') and annotation.__origin__ is Union:
            # For nanint/nanfloat, get the non-literal type (int/float)
            for arg in annotation.__args__:
                if arg is not type(None) and not (hasattr(arg, '__origin__') and arg.__origin__ is Literal):
                    return arg
        return annotation
    return str


def _is_nanfield(field_name: str, dataclass_type: type) -> bool:
    """Check if a field is annotated as nanint or nanfloat."""
    annotations = dataclass_type.__annotations__
    if field_name in annotations:
        annotation = annotations[field_name]
        return annotation == nanint or annotation == nanfloat
    return False


def _read_stats(path: Union[str, Path]) -> Dict[str, str]:
    """WIP
    """
    logger.info("Load statistics from '{}'.".format(path))

    with open(path) as f:
        return dict(line.strip().split('\t', 1) for line in f)


@catch_IOError(logger)
def load_stats(path: Union[str, Path]) -> OutputStats:
    """Load statistics from tab-delimited file and create OutputStats object.

    Args:
        path: Path to the statistics file

    Returns:
        OutputStats object with loaded data

    Raises:
        ValueError: If required fields are missing
    """
    # Read the raw data
    raw_data = _read_stats(path)

    # Parse each section using the integrated function
    summary_data = _parse_stats_section(raw_data, SUMMARY_LABELS, SummaryItems)
    ncc_data = _parse_stats_section(raw_data, NCC_LABELS, CorrelationItems, "NCC")
    mscc_data = _parse_stats_section(raw_data, MSCC_LABELS, CorrelationItems, "MSCC")

    # Construct data objects
    summary_items = SummaryItems(**summary_data)  # type: ignore[arg-type]
    ncc_items = CorrelationItems(**ncc_data)  # type: ignore[arg-type]
    mscc_items = CorrelationItems(**mscc_data)  # type: ignore[arg-type]

    # Create correlation stats with labels
    ncc_stats = CorrelationStats(stats=ncc_items, labels=NCC_LABELS)
    mscc_stats = CorrelationStats(stats=mscc_items, labels=MSCC_LABELS)

    # Create and return OutputStats
    return OutputStats(
        stats=summary_items,
        ncc_stats=ncc_stats,
        mscc_stats=mscc_stats
    )
