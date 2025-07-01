"""
Instrumentation for tracing PyMaSC function calls during golden test execution.
"""

import os
import json
import pickle
import functools
import numpy as np
from pathlib import Path


class FunctionTracer:
    """Records function calls with inputs and outputs."""
    
    def __init__(self, output_dir="/tmp/pymasc_traces"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.call_count = 0
        
    def trace_function(self, func_name):
        """Decorator to trace function calls."""
        def decorator(func):
            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                self.call_count += 1
                call_id = f"{func_name}_{self.call_count}"
                
                # Serialize inputs
                inputs = {
                    'args': self._serialize_args(args),
                    'kwargs': self._serialize_kwargs(kwargs)
                }
                
                # Execute function
                result = func(*args, **kwargs)
                
                # Serialize outputs
                outputs = self._serialize_data(result)
                
                # Save trace
                trace_data = {
                    'function': func_name,
                    'call_id': call_id,
                    'inputs': inputs,
                    'outputs': outputs,
                    'call_order': self.call_count
                }
                
                trace_file = self.output_dir / f"{call_id}.json"
                with open(trace_file, 'w') as f:
                    json.dump(trace_data, f, indent=2)
                
                print(f"TRACE: {func_name} -> {trace_file}")
                return result
            return wrapper
        return decorator
    
    def _serialize_args(self, args):
        """Serialize function arguments."""
        serialized = []
        for i, arg in enumerate(args):
            if i == 0:  # Skip 'self' parameter
                serialized.append("<self>")
            else:
                serialized.append(self._serialize_data(arg))
        return serialized
    
    def _serialize_kwargs(self, kwargs):
        """Serialize keyword arguments."""
        return {k: self._serialize_data(v) for k, v in kwargs.items()}
    
    def _serialize_data(self, data):
        """Convert data to JSON-serializable format."""
        if isinstance(data, np.ndarray):
            return {
                'type': 'numpy_array',
                'shape': data.shape,
                'dtype': str(data.dtype),
                'data': data.tolist() if data.size < 1000 else f"<array_size_{data.size}>"
            }
        elif isinstance(data, dict):
            return {k: self._serialize_data(v) for k, v in data.items()}
        elif isinstance(data, (list, tuple)):
            if len(data) < 100:  # Only serialize small lists
                return [self._serialize_data(item) for item in data]
            else:
                return f"<{type(data).__name__}_length_{len(data)}>"
        elif isinstance(data, (int, float, str, bool, type(None))):
            return data
        else:
            return f"<{type(data).__name__}_{str(data)[:100]}>"


# Global tracer instance
global_tracer = FunctionTracer()