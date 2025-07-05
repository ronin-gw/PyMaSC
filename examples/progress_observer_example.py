#!/usr/bin/env python
"""Example of using the new observer-based progress system with PyMaSC.

This example demonstrates various ways to use the observer pattern
for progress monitoring while maintaining backward compatibility.
"""
import sys
import os
import tempfile
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from PyMaSC.core.progress_migration import (
    enable_progress_migration, 
    ProgressMigrationContext,
    with_progress_observer,
    enable_file_logging,
    create_handler_with_observers
)
from PyMaSC.core.observer import (
    TerminalProgressObserver,
    FileProgressObserver,
    AggregateProgressObserver
)
from PyMaSC.handler.unified import UnifiedCalcHandler
from PyMaSC.core.models import CalculationConfig, AlgorithmType


def example_1_global_migration():
    """Example 1: Enable observer pattern globally."""
    print("\n=== Example 1: Global Migration ===")

    # Enable migration globally - all progress bars will use observers
    enable_progress_migration(True)

    # Now any code using ProgressBar will automatically use observers
    from PyMaSC.utils.progress import ProgressBar

    progress = ProgressBar()
    progress.set("Example 1", 100)
    for i in range(0, 101, 20):
        progress.update(i)
    progress.clean()

    # Disable migration
    enable_progress_migration(False)
    print("Migration disabled - back to normal progress bars")


def example_2_context_manager():
    """Example 2: Use context manager for temporary migration."""
    print("\n=== Example 2: Context Manager ===")

    # Create observer to track events
    observer = TerminalProgressObserver()

    # Use context for temporary migration
    with ProgressMigrationContext() as ctx:
        ctx.attach_observer(observer)

        from PyMaSC.utils.progress import ProgressBar
        progress = ProgressBar()
        progress.set("Context Example", 50)
        progress.update(25)
        progress.clean()

    # Outside context - normal behavior
    print("Outside context - migration disabled")


def example_3_file_logging():
    """Example 3: Log progress to file."""
    print("\n=== Example 3: File Logging ===")

    with tempfile.NamedTemporaryFile(mode='w', suffix='.log', delete=False) as f:
        log_path = f.name

    try:
        # Enable file logging for all progress
        observer = enable_file_logging(log_path, format_type='text')

        from PyMaSC.utils.progress import ProgressBar
        progress = ProgressBar()
        progress.set("File Logging Test", 30)
        for i in range(0, 31, 10):
            progress.update(i)
        progress.clean()

        # Read and display log
        print(f"\nProgress log written to: {log_path}")
        with open(log_path, 'r') as f:
            print("Log contents:")
            print(f.read())

    finally:
        os.unlink(log_path)
        enable_progress_migration(False)


def example_4_decorator():
    """Example 4: Use decorator for automatic observer attachment."""
    print("\n=== Example 4: Decorator ===")

    @with_progress_observer(TerminalProgressObserver, use_progress_bar=False)
    def process_data():
        """Function that uses progress bars."""
        from PyMaSC.utils.progress import ProgressBar

        progress = ProgressBar()
        progress.set("Decorated Function", 40)
        for i in range(0, 41, 10):
            progress.update(i)
        progress.clean()

        return "Processing complete"

    result = process_data()
    print(f"Result: {result}")


def example_5_handler_with_observers():
    """Example 5: Create handler with multiple observers."""
    print("\n=== Example 5: Handler with Observers ===")

    # This would normally use real BAM files
    print("Demonstrating handler observer attachment...")

    # Create configuration
    config = CalculationConfig(
        algorithm=AlgorithmType.SUCCESSIVE,
        max_shift=200,
        mapq_criteria=20
    )

    # Create observers
    terminal_obs = TerminalProgressObserver(use_progress_bar=False)
    aggregate_obs = AggregateProgressObserver()

    # Mock the handler creation (would use real BAM in practice)
    print("Would create handler with observers:")
    print(f"  - Terminal observer: {terminal_obs}")
    print(f"  - Aggregate observer: {aggregate_obs}")
    print("Handler would report progress to all attached observers")


def example_6_aggregate_tracking():
    """Example 6: Track aggregate progress across multiple operations."""
    print("\n=== Example 6: Aggregate Progress Tracking ===")

    from PyMaSC.core.progress_migration import enable_aggregate_tracking

    # Enable aggregate tracking
    observer = enable_aggregate_tracking(["chr1", "chr2", "chr3"])

    # Simulate processing multiple chromosomes
    from PyMaSC.utils.progress import ProgressBar

    for chrom in ["chr1", "chr2", "chr3"]:
        progress = ProgressBar()
        progress.set(chrom, 100)
        for i in range(0, 101, 50):
            progress.update(i)
        progress.clean()

    # Get summary
    summary = observer.get_summary()
    print(f"\nAggregate Summary:")
    print(f"  Total sources: {summary['total_sources']}")
    print(f"  Completed: {summary['completed']}")
    print(f"  Completion rate: {summary['completion_rate']:.1%}")

    enable_progress_migration(False)


def main():
    """Run all examples."""
    print("PyMaSC Progress Observer Examples")
    print("=" * 50)

    examples = [
        example_1_global_migration,
        example_2_context_manager,
        example_3_file_logging,
        example_4_decorator,
        example_5_handler_with_observers,
        example_6_aggregate_tracking
    ]

    for example in examples:
        try:
            example()
        except Exception as e:
            print(f"Error in {example.__name__}: {e}")
        print()

    print("\nAll examples completed!")
    print("\nKey Benefits of Observer Pattern:")
    print("- Decouple progress generation from display")
    print("- Support multiple simultaneous observers")
    print("- Enable progress logging and analysis")
    print("- Maintain full backward compatibility")
    print("- Gradual migration path for existing code")


if __name__ == "__main__":
    main()