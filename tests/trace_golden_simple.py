"""
Simple approach to trace PyMaSC golden calculation by modifying Python wrapper code.
"""

import os
import sys
import tempfile
import subprocess
from pathlib import Path
import shutil


def create_traced_pymasc_script():
    """Create a modified PyMaSC script that logs calculation details."""
    
    script_content = '''#!/usr/bin/env python3
"""
Modified PyMaSC script with calculation tracing.
"""

import sys
import logging
from PyMaSC.pymasc import main
from PyMaSC.core.ncc import NaiveCCCalculator
from PyMaSC.core.mscc import MSCCCalculator

# Set up detailed logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(levelname)s:%(name)s:%(message)s',
    handlers=[
        logging.FileHandler('/tmp/pymasc_detailed.log'),
        logging.StreamHandler(sys.stdout)
    ]
)

# Create a wrapper to capture calculation details
class TracedCalculation:
    def __init__(self):
        self.read_count = 0
        self.forward_positions = []
        self.reverse_positions = []
        
    def log_read(self, strand, ref, pos, readlen):
        self.read_count += 1
        if strand == "forward":
            self.forward_positions.append((ref, pos, readlen))
        else:
            self.reverse_positions.append((ref, pos, readlen))
            
        if self.read_count % 100 == 0:
            print(f"TRACE: Processed {self.read_count} reads")
            
        # Log first 50 reads for detailed analysis
        if self.read_count <= 50:
            print(f"TRACE_READ: {strand} {ref} {pos} {readlen}")

tracer = TracedCalculation()

# Monkey patch the worker class to add tracing
import PyMaSC.handler.masc_worker

original_feed_read = PyMaSC.handler.masc_worker.CalcWorkerBase._feed_read

def traced_feed_read(self, read):
    """Traced version of _feed_read method."""
    if (read.is_read2 or read.mapping_quality < self.mapq_criteria or
            read.is_unmapped or read.is_duplicate):
        return
    
    readlen = read.infer_query_length()
    pos = read.reference_start + 1
    ref = read.reference_name
    
    if read.is_reverse:
        tracer.log_read("reverse", ref, pos, readlen)
    else:
        tracer.log_read("forward", ref, pos, readlen)
    
    # Call original method
    return original_feed_read(self, read)

# Apply the patch
PyMaSC.handler.masc_worker.CalcWorkerBase._feed_read = traced_feed_read

# Run main
if __name__ == "__main__":
    try:
        main()
        print(f"TRACE_SUMMARY: Processed {tracer.read_count} total reads")
        print(f"TRACE_SUMMARY: {len(tracer.forward_positions)} forward reads")
        print(f"TRACE_SUMMARY: {len(tracer.reverse_positions)} reverse reads")
        
        # Save detailed positions for unit tests
        import json
        with open('/tmp/pymasc_read_positions.json', 'w') as f:
            json.dump({
                'forward_positions': tracer.forward_positions[:1000],  # First 1000
                'reverse_positions': tracer.reverse_positions[:1000],   # First 1000
                'total_reads': tracer.read_count
            }, f, indent=2)
        print("TRACE: Saved read positions to /tmp/pymasc_read_positions.json")
        
    except Exception as e:
        print(f"TRACE_ERROR: {e}")
        raise
'''
    
    script_path = "/tmp/traced_pymasc.py"
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    os.chmod(script_path, 0o755)
    return script_path


def run_traced_calculation():
    """Run PyMaSC with tracing enabled."""
    
    script_path = create_traced_pymasc_script()
    
    # Test data paths
    bam_file = "tests/data/ENCFF000RMB-test.bam"
    bigwig_file = "tests/data/hg19_36mer-test.bigwig"
    
    with tempfile.TemporaryDirectory() as temp_dir:
        cmd = [
            'python', script_path,
            bam_file,
            '-m', bigwig_file,
            '-o', temp_dir,
            '-d', '300',
            '-q', '10', 
            '-r', '36'
        ]
        
        print(f"Running traced PyMaSC: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
        
        print("STDOUT:")
        print(result.stdout)
        
        if result.stderr:
            print("STDERR:")
            print(result.stderr)
        
        if result.returncode == 0:
            print("Traced calculation completed successfully")
            
            # Copy trace files
            trace_dir = Path("tests/traces/golden_calculation")
            trace_dir.mkdir(parents=True, exist_ok=True)
            
            if Path("/tmp/pymasc_detailed.log").exists():
                shutil.copy("/tmp/pymasc_detailed.log", trace_dir)
            if Path("/tmp/pymasc_read_positions.json").exists():
                shutil.copy("/tmp/pymasc_read_positions.json", trace_dir)
            
            # Copy output files
            for file_pattern in ["*_stats.tab", "*_cc.tab", "*_mscc.tab"]:
                for result_file in Path(temp_dir).glob(file_pattern):
                    shutil.copy(result_file, trace_dir)
            
            print(f"Trace files saved to: {trace_dir}")
            return True
        else:
            print(f"Traced calculation failed with return code: {result.returncode}")
            return False


if __name__ == "__main__":
    run_traced_calculation()