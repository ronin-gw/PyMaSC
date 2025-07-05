# Handler Migration Guide

This guide explains how to migrate existing PyMaSC handlers to the new simplified, service-oriented architecture.

## Overview

The new architecture separates concerns through dependency injection and service delegation:

- **ValidationService**: All input validation logic
- **IOService**: All file operations
- **MappabilityService**: Mappability-specific operations
- **CalculationService**: Pure calculation logic
- **WorkflowService**: Process orchestration

## Migration Steps

### 1. Replace Direct File Access

**Before:**
```python
class OldHandler:
    def __init__(self, bam_path):
        self.bamfile = AlignmentFile(bam_path)
        # Direct file access throughout
```

**After:**
```python
class NewHandler:
    def __init__(self, bam_path, io_service=None):
        self.io_service = io_service or create_io_service()
        self.bam_info = self.io_service.get_bam_info(bam_path)
```

### 2. Move Validation to Service

**Before:**
```python
def validate_inputs(self):
    if not Path(self.bam_path).exists():
        raise ValueError("File not found")
    # Complex validation logic mixed with handler
```

**After:**
```python
def validate(self):
    result = self.validation_service.validate_workflow_request(
        self.bam_path, self.output_prefix, self.config
    )
    return result.is_valid
```

### 3. Use Dependency Injection

**Before:**
```python
class OldHandler:
    def __init__(self, path, max_shift, mapq):
        # Create all dependencies internally
        self.calculator = NaiveCCCalculator(...)
        self.mappability_handler = MappabilityHandler(...)
```

**After:**
```python
class NewHandler:
    def __init__(self, path, config, dependencies=None):
        self.deps = dependencies or HandlerDependencies()
        self.deps.ensure_services()
```

### 4. Simplify Initialization with Builder

**Before:**
```python
handler = ComplexHandler(
    path=bam_path,
    esttype='mean',
    max_shift=500,
    mapq_criteria=30,
    mappability_path=map_path,
    nworker=4,
    # Many more parameters...
)
```

**After:**
```python
handler = HandlerBuilder() \
    .with_bam_file(bam_path) \
    .with_max_shift(500) \
    .with_mappability(map_path) \
    .with_multiprocessing(4) \
    .build()
```

### 5. Delegate Operations to Services

**Before:**
```python
def run_calculation(self):
    # Complex calculation logic in handler
    for chrom in self.chromosomes:
        reads = self._process_reads(chrom)
        result = self._calculate_correlation(reads)
        # More complex logic...
```

**After:**
```python
def run_calculation(self):
    request = WorkflowRequest(
        bam_path=self.bam_path,
        calculation_config=self.config
    )
    result = self.workflow_service.execute(request)
    return result.calculation_result
```

## Complete Example

### Old Handler (Before)
```python
class CCCalcHandler(BaseCalcHandler):
    def __init__(self, path, esttype, max_shift, mapq_criteria, 
                 nworker=1, mappability_path=None):
        # Complex initialization
        super().__init__(path, mapq_criteria, max_shift, nworker)
        self.esttype = esttype

        # Direct file validation
        if not os.path.exists(path):
            raise ValueError("BAM file not found")

        # Create internal dependencies
        if mappability_path:
            self.mappability_handler = MappabilityHandler(
                mappability_path, max_shift, self.read_len, nworker
            )

        # Complex setup logic
        self._setup_calculators()
        self._validate_chromosomes()

    def run_calculation(self):
        # Complex workflow mixed with business logic
        if self.nworker > 1:
            return self._run_multiprocess()
        else:
            return self._run_singleprocess()
```

### New Handler (After)
```python
class SimplifiedCalcHandler:
    def __init__(self, bam_path, config, dependencies=None):
        self.bam_path = bam_path
        self.config = config
        self.deps = dependencies or HandlerDependencies()
        self.deps.ensure_services()

    def run_calculation(self):
        # Validate
        if not self.validate():
            raise ValueError("Validation failed")

        # Execute via service
        request = WorkflowRequest(
            bam_path=self.bam_path,
            calculation_config=self.config
        )
        result = self.deps.workflow_service.execute(request)

        if not result.is_successful:
            raise RuntimeError(f"Calculation failed: {result.error}")

        return result.calculation_result
```

## Testing Benefits

### Old Testing Approach
```python
def test_handler():
    # Need real files
    with tempfile.NamedTemporaryFile() as f:
        # Complex setup
        handler = OldHandler(f.name, ...)
        # Hard to mock internal dependencies
```

### New Testing Approach
```python
def test_handler():
    # Use mock services
    deps = HandlerDependencies(
        io_service=InMemoryIOService(),
        validation_service=Mock(),
        workflow_service=Mock()
    )

    # Easy to test
    handler = SimplifiedCalcHandler("/fake/path", config, deps)
    result = handler.run_calculation()
```

## Migration Checklist

- [ ] Identify all direct file access → Move to IOService
- [ ] Extract validation logic → Move to ValidationService  
- [ ] Identify mappability operations → Move to MappabilityService
- [ ] Replace internal dependency creation → Use dependency injection
- [ ] Simplify constructor → Use builder pattern
- [ ] Replace complex workflow → Delegate to WorkflowService
- [ ] Update tests → Use mock services
- [ ] Remove obsolete code → Clean up old implementations

## Gradual Migration

You can migrate gradually using the adapter pattern:

```python
# Use ServiceBasedCalcHandler as a bridge
handler = ServiceBasedCalcHandler(
    path=bam_path,
    esttype='mean',
    max_shift=500,
    # Uses services internally but maintains old interface
)
```

This allows existing code to work while you transition to the new architecture.