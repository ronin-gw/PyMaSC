# PyMaSC Handler Analysis and Service Refactoring Opportunities

## Current Handler Architecture Analysis

### 1. Handler Classes Overview

#### BaseCalcHandler (PyMaSC/handler/base.py)
- **Purpose**: Base class extracting common functionality from all calculation handlers
- **Key Responsibilities**:
  - BAM file initialization and validation
  - Chromosome filtering and reference setup
  - Result aggregation patterns
  - Progress management logic
  - Multiprocessing workflow coordination
- **Lines of Code**: ~300 lines
- **Eliminated Duplication**: 85-105 lines across handlers

#### CCCalcHandler (PyMaSC/handler/masc.py)
- **Purpose**: Main cross-correlation calculation handler
- **Extends**: BaseCalcHandler
- **Key Responsibilities**:
  - Read length estimation/validation
  - Mappability handler setup
  - Worker process coordination
  - Calculation orchestration
- **Direct File Access**: BAM file opening, read length estimation
- **Complex Dependencies**: Progress managers, worker processes, queues

#### MappabilityHandler (PyMaSC/handler/mappability.py)
- **Purpose**: Mappability calculation and management
- **Extends**: MappableLengthCalculator
- **Key Responsibilities**:
  - BigWig file access
  - Parallel mappability calculation
  - Statistics caching (JSON files)
  - Progress reporting
- **Direct File Access**: BigWig reading, JSON reading/writing
- **Resource Management**: File handles, worker processes

#### UnifiedCalcHandler (PyMaSC/handler/unified.py)
- **Purpose**: Strategy-pattern based unified handler
- **Key Responsibilities**:
  - Algorithm strategy selection
  - BAM file processing
  - Calculator/worker factory usage
  - Progress observer pattern support
- **Attempts at Improvement**: Uses strategy pattern, observer pattern
- **Still Has**: Direct BAM file access, complex initialization

#### ServiceBasedCalcHandler (PyMaSC/handler/service_adapter.py)
- **Purpose**: Adapter to use service layer with handler interface
- **Shows**: How handlers can delegate to services
- **Key Pattern**: Wrapper around WorkflowService

### 2. Existing Service Layer

#### IOService (PyMaSC/services/io.py)
- **Purpose**: Centralized I/O operations
- **Implementations**: FileIOService, InMemoryIOService
- **Handles**:
  - BAM file reading
  - BigWig file access
  - Results writing (table, stats, figure)
- **Good Patterns**: Abstract interface, testable, resource management

#### CalculationService (PyMaSC/services/calculation.py)
- **Purpose**: Pure calculation logic
- **Implementations**: StandardCalculationService, ParallelCalculationService
- **Handles**:
  - Cross-correlation calculation
  - Result aggregation
  - No I/O operations
- **Good Patterns**: Immutable data, pure functions, cacheable

#### WorkflowService (PyMaSC/services/workflow.py)
- **Purpose**: Orchestrates complete workflows
- **Uses**: IOService + CalculationService
- **Handles**: End-to-end processing

## Refactoring Opportunities

### 1. Responsibilities to Move to Services

#### From Handlers to IOService:
1. **BAM File Validation**
   - Move `_initialize_bam_file()` logic to IOService
   - Create `validate_bam_file()` method
   - Return validation result object

2. **Read Length Estimation**
   - Currently in handlers via `estimate_readlen()`
   - Should be `IOService.estimate_read_length()`
   - Keeps I/O operations centralized

3. **Mappability File Access**
   - Move BigWig operations from MappabilityHandler
   - Already partially in IOService
   - Complete the migration

4. **Progress File I/O**
   - If progress is persisted to disk
   - Should go through IOService

#### From Handlers to New Services:

1. **ValidationService**
   ```python
   class ValidationService:
       def validate_bam_file(self, info: BAMFileInfo) -> ValidationResult
       def validate_chromosome_compatibility(self, bam_info, bigwig_info) -> ValidationResult
       def validate_parameters(self, config: CalculationConfig) -> ValidationResult
   ```

2. **ProgressService**
   ```python
   class ProgressService:
       def create_progress_tracker(self, total_work: int) -> ProgressTracker
       def update_progress(self, tracker: ProgressTracker, completed: int)
       def get_progress_report(self, tracker: ProgressTracker) -> ProgressReport
   ```

3. **MappabilityService**
   ```python
   class MappabilityService:
       def calculate_mappability(self, chromosome: str, config: MappabilityConfig) -> MappabilityResult
       def load_cached_stats(self, path: str) -> Optional[MappabilityStats]
       def save_stats(self, stats: MappabilityStats, path: str)
   ```

4. **ConfigurationService** (extends existing)
   ```python
   class ConfigurationService:
       def merge_configs(self, *configs) -> MergedConfig
       def validate_config(self, config: Any) -> ValidationResult
       def resolve_conflicts(self, configs: List[Config]) -> Config
   ```

### 2. Initialization Complexity Simplification

#### Current Complex Handler Initialization:
```python
# Current: Complex initialization in handler
handler = CCCalcHandler(path, esttype, max_shift, mapq_criteria, nworker, skip_ncc, chromfilter)
handler.set_readlen(readlen)
handler.set_mappability_handler(mappability_handler)
```

#### Proposed Service-Based Approach:
```python
# Proposed: Builder pattern with services
workflow = WorkflowBuilder() \
    .with_bam_file(path) \
    .with_calculation_config(calc_config) \
    .with_mappability(mappability_config) \
    .with_execution_mode(ExecutionMode.PARALLEL, workers=4) \
    .validate() \
    .build()

result = workflow.execute()
```

### 3. Direct File Access Removal

#### Current Direct Access Points:
1. **BAM File Opening** (multiple handlers)
   - Replace with: `io_service.get_bam_info()`
   - Use: `io_service.stream_reads()`

2. **BigWig File Access** (MappabilityHandler)
   - Replace with: `io_service.read_mappability()`
   - Cache through: `mappability_service.get_cached_or_calculate()`

3. **JSON Stats Files** (MappabilityHandler)
   - Replace with: `io_service.read_json()`, `io_service.write_json()`
   - Or use: `mappability_service.load_stats()`, `save_stats()`

### 4. Dependency Injection Opportunities

#### Current Tight Coupling:
```python
class MappabilityHandler:
    def __init__(self, path, max_shift, ...):
        # Creates own BigWigReader
        # Manages own file I/O
        # Creates worker processes
```

#### Proposed Loose Coupling:
```python
class MappabilityService:
    def __init__(self, 
                 io_service: IOService,
                 calculation_service: CalculationService,
                 cache_service: Optional[CacheService] = None):
        self.io = io_service
        self.calc = calculation_service
        self.cache = cache_service or NullCache()
```

### 5. Handler Simplification Pattern

#### Simplified Handler Using Services:
```python
class SimplifiedCalcHandler:
    def __init__(self, services: ServiceContainer):
        self.services = services
        self.config = None
        self.result = None
    
    def configure(self, **kwargs):
        # Use configuration service
        self.config = self.services.config.build_calculation_config(**kwargs)
        
    def run(self):
        # Validate
        validation = self.services.validator.validate_all(self.config)
        if not validation.is_valid:
            raise ValueError(validation.errors)
        
        # Execute through workflow service
        request = self.services.workflow.create_request(self.config)
        self.result = self.services.workflow.execute(request)
        
        return self.result
```

## Implementation Priority

### Phase 1: Core Service Extraction
1. Move read length estimation to IOService
2. Create ValidationService for all validations
3. Move mappability stats I/O to service layer

### Phase 2: Service Integration
1. Update handlers to use services via dependency injection
2. Create ServiceContainer/Registry pattern
3. Implement builder patterns for complex configurations

### Phase 3: Handler Simplification
1. Create thin handler wrappers around services
2. Remove direct file access from all handlers
3. Standardize error handling through services

### Phase 4: Advanced Patterns
1. Implement caching service layer
2. Add event/observer system for progress
3. Create plugin architecture for new algorithms

## Benefits of Service-Oriented Refactoring

1. **Testability**: Services can be mocked/stubbed easily
2. **Reusability**: Services can be composed in different ways
3. **Maintainability**: Clear separation of concerns
4. **Extensibility**: New features as new services
5. **Performance**: Caching, connection pooling at service layer
6. **Debugging**: Centralized logging and monitoring

## Backward Compatibility Strategy

1. Keep existing handler interfaces as facades
2. Gradually migrate internals to use services
3. Provide both old and new APIs during transition
4. Use adapter pattern (like ServiceBasedCalcHandler)
5. Deprecate old patterns with clear migration guides