# Test Data Directory

This directory contains test data files for PyMaSC numerical validation tests.

## ENCODE Test Data

### ChIP-seq Data
- `ENCFF000RMB-test.bam` - Small subset of ENCODE ChIP-seq reads (chr1:700k-850k region)
- `ENCFF000RMB-test.bam.bai` - BAM index file
- `ENCFF000RMB-test.sam` - Human-readable SAM format (reference)

### Mappability Data
- `hg19_36mer-test.bigwig` - 36-mer mappability values for hg19
- `hg19_36mer-test.bedGraph` - Human-readable BedGraph format (reference)

## Test Data Properties

### BAM File Statistics
- Total reads: 2,501
- chr1 reads: 2,501 (concentrated in 700k-850k region)
- Forward reads: ~1,245
- Reverse reads: ~1,256

### Usage in Tests
These files are used in:
- `tests/integration/test_numerical_validation.py`
- `tests/integration/test_expected_results.py`

## Notes
- These are small subsets of real ENCODE data suitable for testing
- The limited region ensures fast test execution while maintaining realistic data properties
- All test data files should remain in this directory for consistency