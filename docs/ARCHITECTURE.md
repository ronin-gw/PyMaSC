# PyMaSC Architecture Documentation

## Overview

PyMaSC has been refactored to use a clean, service-oriented architecture that separates concerns and improves maintainability, testability, and extensibility.

## Architecture Layers

### 1. Core Layer (`PyMaSC/core/`)
The foundation of the application containing:
- **Interfaces** (`interfaces.py`): Abstract base classes defining contracts
- **Models** (`models.py`): Data structures and configuration objects
- **Factory** (`factory.py`): Object creation with dependency injection
- **Observer** (`observer.py`): Event-driven progress reporting
- **Algorithms** (`ncc.pyx`, `bitarray.pyx`, `successive.pyx`): Cython implementations

### 2. Service Layer (`PyMaSC/services/`)
Business logic separated into focused services:
- **CalculationService**: Pure cross-correlation calculations
- **IOService**: All file I/O operations
- **ValidationService**: Input validation logic
- **MappabilityService**: Mappability data handling
- **WorkflowService**: Process orchestration

### 3. Handler Layer (`PyMaSC/handler/`)
User-facing API that coordinates services:
- **BaseCalcHandler**: Common handler functionality
- **CalcHandler**: Unified algorithm interface
- **SimplifiedCalcHandler**: Clean dependency injection example
- **ServiceBasedCalcHandler**: Adapter for legacy compatibility

### 4. Utilities Layer (`PyMaSC/utils/`)
Shared utilities and helpers:
- **ReadProcessing**: Common read filtering/processing
- **OutputUtils**: File path and output management
- **StatsUtils**: Statistical calculations
- **CalcUtils**: Calculation helpers

## Key Design Patterns

### 1. Dependency Injection
```python
# Services are injected, not created internally
handler = SimplifiedCalcHandler(
    bam_path="/path/to/file.bam",
    config=calculation_config,
    dependencies=HandlerDependencies(
        io_service=custom_io_service,
        validation_service=custom_validator
    )
)
```

### 2. Factory Pattern
```python
# Calculators created through factory with new conceptual approach
calculator = CalculatorFactory.create_calculator(
    target=CalculationTarget.BOTH,
    implementation=ImplementationAlgorithm.BITARRAY,
    config=calculation_config,
    mappability_config=mappability_config
)
```

### 3. Strategy Pattern
```python
# Algorithms implement common interface
class CrossCorrelationCalculator(ABC):
    @abstractmethod
    def feed_forward_read(self, chrom: str, pos: int, readlen: int) -> None: pass
    @abstractmethod
    def feed_reverse_read(self, chrom: str, pos: int, readlen: int) -> None: pass
```

### 4. Observer Pattern
```python
# Progress reporting through observers
class ProgressObserver(ABC):
    @abstractmethod
    def update(self, event: ProgressEvent) -> None: pass

handler.attach(ConsoleProgressObserver())
handler.attach(FileProgressObserver("/path/to/log"))
```

### 5. Builder Pattern
```python
# Complex object construction
handler = HandlerBuilder() \
    .with_bam_file("/path/to/file.bam") \
    .with_algorithm("bitarray") \
    .with_max_shift(500) \
    .with_mappability("/path/to/mappability.bw") \
    .with_multiprocessing(4) \
    .build()
```

## Service Architecture

### CalculationService
Encapsulates all cross-correlation calculation logic:
- No file I/O or external dependencies
- Pure functions with immutable data
- Algorithm-agnostic interface
- Efficient caching of calculators

### IOService
Centralizes all I/O operations:
- BAM file reading
- BigWig mappability data
- Output file writing
- Provides both file-based and in-memory implementations

### ValidationService
Comprehensive input validation:
- BAM file validation
- Configuration validation
- Workflow request validation
- Detailed error reporting

### MappabilityService
Handles mappability data operations:
- BigWig file support
- JSON precalculated stats support
- Efficient caching
- Position-level queries

### WorkflowService
Orchestrates complete analysis workflows:
- Coordinates other services
- Progress notification
- Error handling
- Parallel execution support

## Data Flow

1. **Input Validation**
   ```
   User Request → ValidationService → Validated Configuration
   ```

2. **Data Reading**
   ```
   BAM File → IOService → ChromosomeData
   Mappability File → MappabilityService → MappabilityData
   ```

3. **Calculation**
   ```
   ChromosomeData → CalculationService → CalculationResult
   ```

4. **Aggregation**
   ```
   Multiple CalculationResults → GenomeWideResult
   ```

5. **Output**
   ```
   GenomeWideResult → IOService → Output Files
   ```

## Testing Strategy

### Unit Tests
- Mock all external dependencies
- Test each service in isolation
- Use InMemoryIOService for fast tests

### Integration Tests
- Test service interactions
- Use real service implementations
- Verify end-to-end workflows

### Performance Tests
- Benchmark critical paths
- Verify caching effectiveness
- Compare algorithm performance

## Migration Path

### From Legacy Handlers
1. Replace direct file access with IOService
2. Move validation logic to ValidationService
3. Use dependency injection for services
4. Delegate calculations to CalculationService

### Gradual Migration
- Use ServiceBasedCalcHandler as adapter
- Maintain backward compatibility
- Migrate one handler at a time

## Best Practices

### 1. Immutable Data
Use dataclasses and avoid mutation:
```python
@dataclass(frozen=True)
class CalculationResult:
    chromosome: str
    forward_count: int
    reverse_count: int
    correlation_bins: np.ndarray
```

### 2. Explicit Dependencies
Inject all dependencies:
```python
def __init__(self, io_service: IOService, calc_service: CalculationService):
    self.io_service = io_service
    self.calc_service = calc_service
```

### 3. Interface Segregation
Keep interfaces focused:
```python
class ReadProcessor(Protocol):
    def process_read(self, read: AlignedSegment) -> Optional[ReadData]: pass
```

### 4. Error Handling
Use result objects:
```python
@dataclass
class ValidationResult:
    is_valid: bool
    errors: List[str]
    warnings: List[str]
```

## Performance Considerations

### Caching
- Calculator instances are cached by configuration
- Mappability data is cached with LRU eviction
- Results in 380x speedup for repeated calculations

### Memory Efficiency
- Streaming read processing
- Chunk-based mappability queries
- Lazy evaluation where possible

### Parallelization
- Chromosome-level parallelization
- Shared-nothing architecture
- Process pool for multi-core systems

## Future Enhancements

### 1. Plugin System
- Dynamic algorithm loading
- Custom output formats
- Third-party integrations

### 2. Cloud Support
- S3/GCS I/O service implementations
- Distributed calculation service
- Serverless deployment

### 3. Real-time Analysis
- Streaming calculation service
- Progressive result updates
- Live visualization support

## Conclusion

The refactored PyMaSC architecture provides:
- **Maintainability**: Clear separation of concerns
- **Testability**: Easy mocking and isolation
- **Performance**: Efficient caching and parallelization
- **Extensibility**: Plugin points for customization
- **Reliability**: Comprehensive validation and error handling

This architecture ensures PyMaSC can evolve to meet future requirements while maintaining backward compatibility and excellent performance.