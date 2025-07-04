"""Mappability service for handling mappability data operations.

This module provides the MappabilityService that centralizes all
mappability-related operations, removing them from handlers and
ensuring consistent mappability handling across the application.

Key features:
- BigWig/JSON mappability data reading
- Mappability statistics calculation
- Precalculation support
- Caching for performance
"""
import json
import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

import numpy as np

from PyMaSC.reader.bigwig import BigWigReader
from PyMaSC.core.models import MappabilityConfig

logger = logging.getLogger(__name__)


@dataclass
class MappabilityStats:
    """Statistics for mappability data.
    
    Attributes:
        chromosome: Chromosome name
        mappable_length: Length of mappable regions
        total_length: Total chromosome length
        mappable_fraction: Fraction of mappable bases
        mean_mappability: Mean mappability value
        median_mappability: Median mappability value
    """
    chromosome: str
    mappable_length: int
    total_length: int
    mappable_fraction: float
    mean_mappability: float
    median_mappability: float


@dataclass
class GenomeWideMappabilityStats:
    """Genome-wide mappability statistics.
    
    Attributes:
        chromosome_stats: Per-chromosome statistics
        total_mappable_length: Total mappable bases genome-wide
        total_genome_length: Total genome length
        genome_mappable_fraction: Overall mappable fraction
    """
    chromosome_stats: Dict[str, MappabilityStats]
    total_mappable_length: int
    total_genome_length: int
    genome_mappable_fraction: float


class MappabilityService(ABC):
    """Abstract base for mappability services.
    
    Defines the interface for all mappability operations,
    enabling different implementations and strategies.
    """
    
    @abstractmethod
    def get_mappability_values(self,
                             chromosome: str,
                             start: int,
                             end: int) -> np.ndarray:
        """Get mappability values for a genomic region.
        
        Args:
            chromosome: Chromosome name
            start: Start position (1-based)
            end: End position (1-based, inclusive)
            
        Returns:
            Array of mappability values
        """
        pass
    
    @abstractmethod
    def calculate_mappability_stats(self,
                                  chromosome: str,
                                  length: int) -> MappabilityStats:
        """Calculate mappability statistics for a chromosome.
        
        Args:
            chromosome: Chromosome name
            length: Chromosome length
            
        Returns:
            Mappability statistics
        """
        pass
    
    @abstractmethod
    def calculate_genome_wide_stats(self,
                                  chromosome_lengths: Dict[str, int]) -> GenomeWideMappabilityStats:
        """Calculate genome-wide mappability statistics.
        
        Args:
            chromosome_lengths: Dictionary of chromosome lengths
            
        Returns:
            Genome-wide statistics
        """
        pass
    
    @abstractmethod
    def is_position_mappable(self,
                           chromosome: str,
                           position: int,
                           threshold: float = 0.5) -> bool:
        """Check if a position is mappable.
        
        Args:
            chromosome: Chromosome name
            position: Genomic position (1-based)
            threshold: Mappability threshold
            
        Returns:
            True if position is mappable
        """
        pass


class BigWigMappabilityService(MappabilityService):
    """Mappability service using BigWig files.
    
    Provides efficient access to mappability data stored
    in BigWig format with caching for performance.
    """
    
    def __init__(self, 
                 bigwig_path: str,
                 config: Optional[MappabilityConfig] = None,
                 cache_size: int = 100):
        """Initialize BigWig mappability service.
        
        Args:
            bigwig_path: Path to BigWig file
            config: Optional mappability configuration
            cache_size: Number of regions to cache
        """
        self.bigwig_path = bigwig_path
        self.config = config or MappabilityConfig(
            mappability_path=bigwig_path,
            read_len=50
        )
        self.cache_size = cache_size
        
        # Initialize reader
        self._reader = None
        self._cache = {}  # Simple LRU-style cache
        self._cache_order = []
        
    def _get_reader(self) -> BigWigReader:
        """Get or create BigWig reader."""
        if self._reader is None:
            self._reader = BigWigReader(self.bigwig_path)
        return self._reader
    
    def get_mappability_values(self,
                             chromosome: str,
                             start: int,
                             end: int) -> np.ndarray:
        """Get mappability values for a genomic region."""
        # Check cache
        cache_key = (chromosome, start, end)
        if cache_key in self._cache:
            # Move to end (LRU)
            self._cache_order.remove(cache_key)
            self._cache_order.append(cache_key)
            return self._cache[cache_key].copy()
        
        # Read from BigWig
        reader = self._get_reader()
        
        try:
            # BigWigReader expects 0-based coordinates
            values = reader.get_as_array(chromosome, start - 1, end)
            
            # Convert None values to 0.0
            result = np.array([v if v is not None else 0.0 for v in values])
            
            # Update cache
            self._update_cache(cache_key, result)
            
            return result
            
        except Exception as e:
            logger.error(f"Error reading mappability for {chromosome}:{start}-{end}: {e}")
            # Return zeros on error
            return np.zeros(end - start)
    
    def _update_cache(self, key: Tuple, value: np.ndarray) -> None:
        """Update cache with LRU eviction."""
        if len(self._cache) >= self.cache_size:
            # Evict oldest
            oldest = self._cache_order.pop(0)
            del self._cache[oldest]
        
        self._cache[key] = value.copy()
        self._cache_order.append(key)
    
    def calculate_mappability_stats(self,
                                  chromosome: str,
                                  length: int) -> MappabilityStats:
        """Calculate mappability statistics for a chromosome."""
        reader = self._get_reader()
        
        try:
            # Get chromosome info
            chroms = reader.chroms()
            if chromosome not in chroms:
                logger.warning(f"Chromosome {chromosome} not found in mappability file")
                return MappabilityStats(
                    chromosome=chromosome,
                    mappable_length=0,
                    total_length=length,
                    mappable_fraction=0.0,
                    mean_mappability=0.0,
                    median_mappability=0.0
                )
            
            # Calculate stats in chunks to avoid memory issues
            chunk_size = 1000000  # 1 Mb chunks
            mappable_bases = 0
            all_values = []
            
            for start in range(0, length, chunk_size):
                end = min(start + chunk_size, length)
                values = self.get_mappability_values(chromosome, start + 1, end + 1)
                
                # Count mappable bases (using read length threshold)
                threshold = 1.0 / (2 * self.config.read_len)
                mappable_bases += np.sum(values > threshold)
                
                # Collect non-zero values for statistics
                non_zero = values[values > 0]
                if len(non_zero) > 0:
                    all_values.extend(non_zero)
            
            # Calculate statistics
            if all_values:
                mean_map = np.mean(all_values)
                median_map = np.median(all_values)
            else:
                mean_map = 0.0
                median_map = 0.0
            
            return MappabilityStats(
                chromosome=chromosome,
                mappable_length=int(mappable_bases),
                total_length=length,
                mappable_fraction=mappable_bases / length if length > 0 else 0.0,
                mean_mappability=float(mean_map),
                median_mappability=float(median_map)
            )
            
        except Exception as e:
            logger.error(f"Error calculating mappability stats for {chromosome}: {e}")
            return MappabilityStats(
                chromosome=chromosome,
                mappable_length=0,
                total_length=length,
                mappable_fraction=0.0,
                mean_mappability=0.0,
                median_mappability=0.0
            )
    
    def calculate_genome_wide_stats(self,
                                  chromosome_lengths: Dict[str, int]) -> GenomeWideMappabilityStats:
        """Calculate genome-wide mappability statistics."""
        chromosome_stats = {}
        total_mappable = 0
        total_length = 0
        
        for chrom, length in chromosome_lengths.items():
            stats = self.calculate_mappability_stats(chrom, length)
            chromosome_stats[chrom] = stats
            total_mappable += stats.mappable_length
            total_length += stats.total_length
        
        return GenomeWideMappabilityStats(
            chromosome_stats=chromosome_stats,
            total_mappable_length=total_mappable,
            total_genome_length=total_length,
            genome_mappable_fraction=total_mappable / total_length if total_length > 0 else 0.0
        )
    
    def is_position_mappable(self,
                           chromosome: str,
                           position: int,
                           threshold: float = 0.5) -> bool:
        """Check if a position is mappable."""
        values = self.get_mappability_values(chromosome, position, position + 1)
        return len(values) > 0 and values[0] > threshold
    
    def close(self) -> None:
        """Close resources."""
        if self._reader is not None:
            self._reader.close()
            self._reader = None
        self._cache.clear()
        self._cache_order.clear()
    
    def __enter__(self):
        """Context manager entry."""
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()


class JSONMappabilityService(MappabilityService):
    """Mappability service using pre-calculated JSON files.
    
    Provides access to pre-calculated mappability statistics
    stored in JSON format.
    """
    
    def __init__(self, json_path: str):
        """Initialize JSON mappability service.
        
        Args:
            json_path: Path to JSON file
        """
        self.json_path = json_path
        self._data = None
        self._load_data()
    
    def _load_data(self) -> None:
        """Load JSON data."""
        try:
            with open(self.json_path) as f:
                self._data = json.load(f)
        except Exception as e:
            logger.error(f"Error loading JSON mappability file: {e}")
            self._data = {'stat': {}, 'chromosomes': {}}
    
    def get_mappability_values(self,
                             chromosome: str,
                             start: int,
                             end: int) -> np.ndarray:
        """Get mappability values (not supported for JSON)."""
        logger.warning("get_mappability_values not supported for JSON format")
        return np.ones(end - start)  # Assume all mappable
    
    def calculate_mappability_stats(self,
                                  chromosome: str,
                                  length: int) -> MappabilityStats:
        """Get pre-calculated mappability statistics."""
        chrom_data = self._data.get('chromosomes', {}).get(chromosome, {})
        
        if chrom_data:
            return MappabilityStats(
                chromosome=chromosome,
                mappable_length=chrom_data.get('mappable_length', 0),
                total_length=length,
                mappable_fraction=chrom_data.get('mappable_fraction', 0.0),
                mean_mappability=chrom_data.get('mean_mappability', 0.0),
                median_mappability=chrom_data.get('median_mappability', 0.0)
            )
        else:
            return MappabilityStats(
                chromosome=chromosome,
                mappable_length=0,
                total_length=length,
                mappable_fraction=0.0,
                mean_mappability=0.0,
                median_mappability=0.0
            )
    
    def calculate_genome_wide_stats(self,
                                  chromosome_lengths: Dict[str, int]) -> GenomeWideMappabilityStats:
        """Get pre-calculated genome-wide statistics."""
        # Use pre-calculated stats if available
        if 'stat' in self._data:
            stat = self._data['stat']
            return GenomeWideMappabilityStats(
                chromosome_stats={
                    chrom: self.calculate_mappability_stats(chrom, length)
                    for chrom, length in chromosome_lengths.items()
                },
                total_mappable_length=stat.get('total_mappable_length', 0),
                total_genome_length=stat.get('total_genome_length', 0),
                genome_mappable_fraction=stat.get('genome_mappable_fraction', 0.0)
            )
        else:
            # Calculate from chromosome data
            chromosome_stats = {}
            total_mappable = 0
            total_length = 0
            
            for chrom, length in chromosome_lengths.items():
                stats = self.calculate_mappability_stats(chrom, length)
                chromosome_stats[chrom] = stats
                total_mappable += stats.mappable_length
                total_length += stats.total_length
            
            return GenomeWideMappabilityStats(
                chromosome_stats=chromosome_stats,
                total_mappable_length=total_mappable,
                total_genome_length=total_length,
                genome_mappable_fraction=total_mappable / total_length if total_length > 0 else 0.0
            )
    
    def is_position_mappable(self,
                           chromosome: str,
                           position: int,
                           threshold: float = 0.5) -> bool:
        """Check if position is mappable (always True for JSON)."""
        return True  # JSON format doesn't have position-level data


# Factory function for service creation
def create_mappability_service(path: str,
                             config: Optional[MappabilityConfig] = None) -> MappabilityService:
    """Create appropriate mappability service based on file type.
    
    Args:
        path: Path to mappability file
        config: Optional mappability configuration
        
    Returns:
        Mappability service instance
    """
    path_obj = Path(path)
    
    if path_obj.suffix.lower() in ['.bw', '.bigwig']:
        return BigWigMappabilityService(path, config)
    elif path_obj.suffix.lower() == '.json':
        return JSONMappabilityService(path)
    else:
        raise ValueError(f"Unsupported mappability file format: {path_obj.suffix}")