# PyMaSC Legacy Code Cleanup Plan

## Overview

This document outlines the plan for removing legacy code after the architectural refactoring. The cleanup should be done in phases to ensure backward compatibility is maintained for users during the transition period.

## Cleanup Phases

### Phase 1: Compatibility Layer Removal (After 3-month deprecation period)

These files exist solely for backward compatibility and should be removed first:

1. **`PyMaSC/core/worker_compat.py`**
   - Purpose: Backward compatibility for old worker classes
   - Dependencies: Check for imports in user scripts
   - Action: Remove after deprecation warnings have been shown for 3 months

2. **`PyMaSC/core/progress_migration.py`**
   - Purpose: Migration utilities for progress system
   - Dependencies: Used by legacy handlers
   - Action: Remove once all handlers use observer pattern

3. **`PyMaSC/core/progress_adapter.py`**
   - Purpose: Bridge between old and new progress systems
   - Dependencies: Used by compatibility layers
   - Action: Remove with other compatibility code

4. **`PyMaSC/handler/compat.py`**
   - Purpose: Legacy handler interface compatibility
   - Dependencies: May be imported by user scripts
   - Action: Remove after deprecation period

### Phase 2: Deprecated Worker Implementations

These workers have been replaced by UnifiedWorker:

1. **`PyMaSC/handler/masc_worker.py`**
   - Status: Already marked as deprecated
   - Replacement: UnifiedWorker
   - Action: Remove and update any remaining references

2. **`PyMaSC/handler/masc_noindex_worker.py`**
   - Purpose: Single-process for non-indexed BAMs
   - Replacement: Unified architecture handles this
   - Action: Verify functionality is covered, then remove

### Phase 3: Transitional Components

1. **`PyMaSC/handler/service_adapter.py`**
   - Purpose: Adapter for service layer integration
   - Future: Not needed once all code uses services directly
   - Action: Keep until full service migration

### Phase 4: Legacy Handlers (Requires careful review)

1. **`PyMaSC/handler/masc.py`** (CCCalcHandler)
2. **`PyMaSC/handler/bamasc.py`** (BACalcHandler)
   - Status: Functionality moved to UnifiedCalcHandler
   - Risk: May still be used by external scripts
   - Action: Add deprecation warnings first, remove later

## Pre-Cleanup Checklist

Before removing any file:

- [ ] Search entire codebase for imports
- [ ] Check documentation for references
- [ ] Verify replacement functionality works
- [ ] Add deprecation warnings (if not present)
- [ ] Update migration guide
- [ ] Test with common use cases

## Code Search Commands

```bash
# Find imports of a module
grep -r "from PyMaSC.core.worker_compat import" .
grep -r "import PyMaSC.core.worker_compat" .

# Find class usage
grep -r "NaiveCCCalcWorker" . --exclude-dir=.git
grep -r "CCCalcHandler" . --exclude-dir=.git

# Find test dependencies
find tests/ -name "*.py" -exec grep -l "worker_compat\|progress_migration" {} \;
```

## Migration Messages

Add clear deprecation messages before removal:

```python
import warnings

warnings.warn(
    "CCCalcHandler is deprecated and will be removed in PyMaSC 1.0. "
    "Please use UnifiedCalcHandler instead. "
    "See migration guide: https://github.com/bfichera/PyMaSC/docs/migration.md",
    DeprecationWarning,
    stacklevel=2
)
```

## Post-Cleanup Tasks

After removing legacy code:

1. **Update Documentation**
   - Remove references to deprecated classes
   - Update API documentation
   - Update examples

2. **Update Tests**
   - Remove tests for deleted code
   - Ensure coverage remains high
   - Update integration tests

3. **Update Dependencies**
   - Remove unused imports
   - Update setup.py if needed
   - Clean requirements files

4. **Communication**
   - Announce removals in release notes
   - Update migration guide
   - Notify major users

## Risk Mitigation

1. **Gradual Removal**: Remove in small batches
2. **Version Tags**: Tag version before major removals
3. **Rollback Plan**: Keep removed code in a branch
4. **User Testing**: Beta release before final removal

## Timeline Recommendation

- **Month 1-3**: Add deprecation warnings
- **Month 4-6**: Remove Phase 1 compatibility layers
- **Month 7-9**: Remove Phase 2 deprecated workers
- **Month 10-12**: Remove remaining legacy code

## Conclusion

The cleanup should be done gradually with clear communication to users. The new architecture is fully functional, but rushing the removal of legacy code could break existing workflows. A measured approach ensures a smooth transition for all users.