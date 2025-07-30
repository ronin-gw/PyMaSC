"""Test progress migration utilities."""

from unittest.mock import Mock, patch
import tempfile
import os

from PyMaSC.core.progress_migration import (
    enable_progress_migration, with_progress_observer,
    ProgressMigrationContext, create_handler_with_observers,
    enable_file_logging, enable_aggregate_tracking,
    auto_attach_observer
)
from PyMaSC.core.observer import (
    ProgressObserver, ProgressEventType, FileProgressObserver
)
from PyMaSC.core.progress_adapter import ProgressBarAdapter, ReadCountProgressBarAdapter


class TestProgressMigration:
    """Test progress migration functionality."""

    def setup_method(self):
        """Reset migration state before each test."""
        enable_progress_migration(False)

    def teardown_method(self):
        """Ensure migration is disabled after tests."""
        enable_progress_migration(False)

    def test_migration_toggle(self):
        """Test enabling and disabling migration."""
        import PyMaSC.utils.progress as progress_module

        # Store original class
        original_class = progress_module.ProgressBar

        # Enable migration
        enable_progress_migration(True)

        # Create progress bar - should be adapter
        pb = progress_module.ProgressBar()
        assert isinstance(pb, ProgressBarAdapter)

        # Disable migration
        enable_progress_migration(False)

        # Create progress bar - should be original
        pb2 = progress_module.ProgressBar()
        assert pb2.__class__.__name__ == 'ProgressBar'
        assert not isinstance(pb2, ProgressBarAdapter)

    def test_auto_attach_observer(self):
        """Test automatic observer attachment."""
        # Create test observer
        events = []

        class TestObserver(ProgressObserver):
            def update(self, event):
                events.append(event)

        observer = TestObserver()

        # Enable migration and auto-attach
        enable_progress_migration(True)
        auto_attach_observer(observer)

        # Create progress bar
        import PyMaSC.utils.progress as progress_module
        pb = progress_module.ProgressBar()

        # Use progress bar
        pb.set("test", 100)

        # Check observer was notified
        assert len(events) == 1
        assert events[0].event_type == ProgressEventType.STARTED

    def test_progress_migration_context(self):
        """Test migration context manager."""
        import PyMaSC.utils.progress as progress_module

        # Outside context - should be original
        pb1 = progress_module.ProgressBar()
        assert not isinstance(pb1, ProgressBarAdapter)

        # Inside context - should be adapter
        with ProgressMigrationContext() as ctx:
            pb2 = progress_module.ProgressBar()
            assert isinstance(pb2, ProgressBarAdapter)

        # After context - should be original again
        pb3 = progress_module.ProgressBar()
        assert not isinstance(pb3, ProgressBarAdapter)

    def test_context_with_observer(self):
        """Test context manager with observer attachment."""
        events = []

        class TestObserver(ProgressObserver):
            def update(self, event):
                events.append(event)

        observer = TestObserver()

        # Create context and immediately attach observer
        ctx = ProgressMigrationContext()
        ctx.attach_observer(observer)

        with ctx:
            import PyMaSC.utils.progress as progress_module
            pb = progress_module.ProgressBar()
            pb.set("test", 100)
            pb.clean()

        # Observer should have received events
        assert len(events) >= 1  # At least start event
        event_types = [e.event_type for e in events]
        assert ProgressEventType.STARTED in event_types

    def test_with_progress_observer_decorator(self):
        """Test decorator for automatic observer attachment."""
        # Track events
        events = []

        class TestObserver(ProgressObserver):
            def update(self, event):
                events.append(event)

        @with_progress_observer(TestObserver)
        def process_with_progress():
            import PyMaSC.utils.progress as progress_module
            pb = progress_module.ProgressBar()
            pb.set("decorated", 50)
            pb.update(25)
            pb.clean()
            return "done"

        # Call decorated function
        result = process_with_progress()
        assert result == "done"

        # Should have recorded events
        assert len(events) >= 2  # At least start and complete

    def test_create_handler_with_observers(self):
        """Test creating handler with observers."""
        # Mock handler class
        class MockHandler:
            def __init__(self, path):
                self.path = path
                self._observers = []

            def attach_progress_observer(self, observer):
                self._observers.append(observer)

        # Create observers
        obs1 = Mock(spec=ProgressObserver)
        obs2 = Mock(spec=ProgressObserver)

        # Create handler with observers
        handler = create_handler_with_observers(
            MockHandler,
            [obs1, obs2],
            path="test.bam"
        )

        assert handler.path == "test.bam"
        assert len(handler._observers) == 2
        assert obs1 in handler._observers
        assert obs2 in handler._observers

    def test_enable_file_logging(self):
        """Test enabling file logging convenience function."""
        # Test that function works without errors
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            temp_path = f.name

        try:
            # Enable file logging
            observer = enable_file_logging(temp_path, 'json')

            # Should be a FileProgressObserver
            assert isinstance(observer, FileProgressObserver)

            # Migration should be enabled
            import PyMaSC.core.progress_migration as pm
            assert pm._migration_enabled is True

        finally:
            try:
                os.unlink(temp_path)
            except:
                pass

    def test_enable_aggregate_tracking(self):
        """Test enabling aggregate tracking."""
        # Enable aggregate tracking
        observer = enable_aggregate_tracking(["chr1", "chr2"])

        # Create progress bars
        import PyMaSC.utils.progress as progress_module

        pb1 = progress_module.ProgressBar()
        pb1.set("chr1", 100)
        pb1.update(50)
        pb1.clean()

        pb2 = progress_module.ProgressBar()
        pb2.set("chr2", 200)
        pb2.update(100)
        pb2.clean()

        # Check aggregate summary
        summary = observer.get_summary()
        assert summary['completed'] == 2

    def test_read_count_progress_adapter(self):
        """Test ReadCountProgressBar adapter."""
        enable_progress_migration(True)

        events = []

        class TestObserver(ProgressObserver):
            def update(self, event):
                events.append(event)

        observer = TestObserver()
        auto_attach_observer(observer)

        # Create ReadCountProgressBar
        import PyMaSC.utils.progress as progress_module
        rcpb = progress_module.ReadCountProgressBar()

        # Should be adapter
        assert isinstance(rcpb, ReadCountProgressBarAdapter)

        # Attach observer directly
        rcpb.attach_observer(observer)

        # Use it
        rcpb.set_genome(1000000)
        rcpb.set_chrom(100000, "chr1")
        rcpb.update(50000)
        rcpb.finish()

        # Check events
        event_types = [e.event_type for e in events]
        assert ProgressEventType.STARTED in event_types
        assert ProgressEventType.UPDATED in event_types
        assert ProgressEventType.COMPLETED in event_types

        # Check sources
        sources = [e.source for e in events]
        assert "genome" in sources
        assert "chromosome:chr1" in sources


class TestMigrationIntegration:
    """Test integration with existing PyMaSC components."""

    def setup_method(self):
        """Reset state."""
        enable_progress_migration(False)

    def teardown_method(self):
        """Clean up."""
        enable_progress_migration(False)

    def test_readlen_integration(self):
        """Test integration with read length estimation."""
        from PyMaSC.core.readlen_enhanced import ReadLengthEstimator

        # Test that estimator can be created and observers attached
        estimator = ReadLengthEstimator(use_observer=True)

        # Track events
        events = []

        class TestObserver(ProgressObserver):
            def update(self, event):
                events.append(event)

        # Should be able to attach observer
        estimator.attach_observer(TestObserver())
        assert len(estimator._observers) == 1

        # Test observer can be detached
        estimator.detach_observer(estimator._observers[0])
        assert len(estimator._observers) == 0

    def test_unified_handler_integration(self):
        """Test integration with CalcHandler."""
        from PyMaSC.handler.calc import CalcHandler
        from PyMaSC.core.models import (
            CalculationConfig, CalculationTarget, ImplementationAlgorithm,
            ExecutionConfig, ExecutionMode
        )

        # Mock file
        with patch('pysam.AlignmentFile') as mock_af:
            mock_bam = Mock()
            mock_bam.references = ['chr1']
            mock_bam.lengths = [1000]
            mock_bam.has_index.return_value = False
            mock_af.return_value = mock_bam

            # Create handler
            config = CalculationConfig(
                target=CalculationTarget.NCC,
                implementation=ImplementationAlgorithm.SUCCESSIVE,
                max_shift=200,
                mapq_criteria=20
            )

            exec_config = ExecutionConfig(
                mode=ExecutionMode.SINGLE_PROCESS,
                worker_count=1
            )

            handler = CalcHandler("test.bam", config, exec_config)

            # Check observer support
            assert hasattr(handler, 'attach_progress_observer')
            assert hasattr(handler, 'use_observer_progress')
            assert handler.use_observer_progress is True

            # Attach observer
            observer = Mock(spec=ProgressObserver)
            handler.attach_progress_observer(observer)
            assert observer in handler._progress_observers