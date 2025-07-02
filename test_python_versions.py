#!/usr/bin/env python
"""Test PyMaSC compatibility across Python versions"""
import sys
import subprocess
import os

PYTHON_VERSIONS = ["3.8.9", "3.9.21", "3.10.15", "3.11.1", "3.12.11", "3.13.3"]

def test_version(version):
    """Test PyMaSC with a specific Python version"""
    print(f"\n{'='*60}")
    print(f"Testing Python {version}")
    print('='*60)
    
    env = os.environ.copy()
    env['PYENV_VERSION'] = version
    
    # Test basic import
    cmd = ["python", "test_basic_imports.py"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, env=env)
        if result.returncode == 0:
            print(result.stdout)
            print(f"✅ Python {version}: Basic imports passed")
            return True
        else:
            print(result.stdout)
            print(result.stderr)
            print(f"❌ Python {version}: Basic imports failed")
            return False
    except Exception as e:
        print(f"❌ Python {version}: Error running test: {e}")
        return False

def main():
    results = {}
    for version in PYTHON_VERSIONS:
        # Check if version is installed
        check_cmd = ["pyenv", "versions", "--bare"]
        installed = subprocess.run(check_cmd, capture_output=True, text=True)
        if version not in installed.stdout:
            print(f"⚠️  Python {version} not installed, skipping...")
            continue
        
        results[version] = test_version(version)
    
    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print('='*60)
    for version, passed in results.items():
        status = "✅ PASSED" if passed else "❌ FAILED"
        print(f"Python {version}: {status}")

if __name__ == "__main__":
    main()