"""Generate minimal test data for PyMaSC unit tests."""

import os
import tempfile
import numpy as np
from typing import List, Tuple, Dict


def create_mock_read_data(n_reads: int = 100, max_shift: int = 200) -> Tuple[List[int], List[int]]:
    """
    Create mock read position data for testing cross-correlation algorithms.
    
    Args:
        n_reads: Number of reads to generate
        max_shift: Maximum shift distance for cross-correlation
        
    Returns:
        Tuple of (forward_positions, reverse_positions)
    """
    np.random.seed(42)  # Reproducible results
    
    # Generate reads with realistic fragment length distribution
    fragment_length = 150  # Expected fragment length
    read_length = 50
    
    forward_positions = []
    reverse_positions = []
    
    for i in range(n_reads):
        start_pos = np.random.randint(1000, 10000)
        
        # Forward read at start position
        forward_positions.append(start_pos)
        
        # Reverse read at end position (start + fragment_length - read_length)
        reverse_pos = start_pos + fragment_length - read_length
        reverse_positions.append(reverse_pos)
        
        # Add some noise reads
        if i % 10 == 0:
            noise_pos = np.random.randint(1000, 10000)
            if np.random.random() > 0.5:
                forward_positions.append(noise_pos)
            else:
                reverse_positions.append(noise_pos)
    
    return sorted(forward_positions), sorted(reverse_positions)


def create_mock_reference_data() -> Tuple[List[str], List[int]]:
    """
    Create mock reference chromosome data.
    
    Returns:
        Tuple of (chromosome_names, chromosome_lengths)
    """
    references = ['chr1', 'chr2', 'chr3']
    lengths = [20000, 15000, 10000]
    return references, lengths


def create_mock_mappability_data(chrom_size: int = 20000, 
                                 mappable_fraction: float = 0.8) -> np.ndarray:
    """
    Create mock mappability data for testing.
    
    Args:
        chrom_size: Size of chromosome
        mappable_fraction: Fraction of positions that are mappable
        
    Returns:
        Array of mappability values (0.0 to 1.0)
    """
    np.random.seed(42)
    
    # Create basic mappability pattern
    mappability = np.ones(chrom_size, dtype=np.float32)
    
    # Add unmappable regions (repetitive sequences, etc.)
    n_unmappable = int(chrom_size * (1 - mappable_fraction))
    unmappable_positions = np.random.choice(chrom_size, n_unmappable, replace=False)
    mappability[unmappable_positions] = 0.0
    
    # Add partially mappable regions
    n_partial = int(chrom_size * 0.1)
    partial_positions = np.random.choice(chrom_size, n_partial, replace=False)
    mappability[partial_positions] = np.random.uniform(0.3, 0.7, n_partial)
    
    return mappability


def create_expected_cc_result(max_shift: int = 200, 
                              fragment_length: int = 150) -> Dict[int, float]:
    """
    Create expected cross-correlation results for validation.
    
    Args:
        max_shift: Maximum shift distance
        fragment_length: Expected fragment length
        
    Returns:
        Dictionary mapping shift distance to CC value
    """
    cc_values = {}
    
    for shift in range(max_shift + 1):
        if shift == fragment_length:
            # Peak at fragment length
            cc_values[shift] = 1.0
        elif abs(shift - fragment_length) < 20:
            # Values near peak
            distance = abs(shift - fragment_length)
            cc_values[shift] = 1.0 - (distance / 20.0) * 0.8
        else:
            # Background noise
            cc_values[shift] = 0.1 + np.random.random() * 0.1
    
    return cc_values


class MockBAMData:
    """Mock BAM-like data structure for testing."""
    
    def __init__(self, references: List[str], lengths: List[int]):
        self.references = references
        self.lengths = lengths
        self.reads = self._generate_reads()
    
    def _generate_reads(self) -> List[Dict]:
        """Generate mock read data."""
        reads = []
        
        for ref_idx, (ref_name, ref_length) in enumerate(zip(self.references, self.lengths)):
            n_reads = min(100, ref_length // 100)  # Reasonable read density
            
            for i in range(n_reads):
                pos = np.random.randint(1, ref_length - 50)
                
                read = {
                    'reference_id': ref_idx,
                    'reference_name': ref_name,
                    'reference_start': pos,
                    'reference_end': pos + 50,  # 50bp read length
                    'is_reverse': np.random.random() > 0.5,
                    'mapping_quality': np.random.randint(20, 60),
                    'query_length': 50
                }
                reads.append(read)
        
        return reads
    
    def fetch(self, reference: str = None):
        """Mock fetch method like pysam.AlignmentFile."""
        if reference:
            return [read for read in self.reads if read['reference_name'] == reference]
        return self.reads


def create_temp_test_file(content: bytes, suffix: str = '.tmp') -> str:
    """
    Create a temporary file with given content.
    
    Args:
        content: Binary content to write
        suffix: File suffix
        
    Returns:
        Path to temporary file
    """
    fd, path = tempfile.mkstemp(suffix=suffix)
    try:
        with os.fdopen(fd, 'wb') as f:
            f.write(content)
    except:
        os.close(fd)
        raise
    
    return path


def cleanup_temp_file(path: str):
    """Clean up temporary file."""
    try:
        os.unlink(path)
    except OSError:
        pass