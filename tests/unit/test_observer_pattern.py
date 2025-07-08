"""Test observer pattern implementation for progress reporting."""

import pytest
from unittest.mock import Mock, patch
import json
import tempfile
import os

from PyMaSC.core.observer import (
    ProgressEvent, ProgressEventType, ProgressObserver,
    ProgressSubject, TerminalProgressObserver,
    MultiprocessProgressObserver, FileProgressObserver,
    AggregateProgressObserver
)
from PyMaSC.core.progress_adapter import (
    ProgressBarAdapter,
    LegacyProgressObserver, ProgressManager,
    get_progress_manager, reset_progress_manager
)


class TestProgressEvent:
    """Test ProgressEvent data structure."""

    def test_event_creation_basic(self):
        """Test basic event creation."""
        event = ProgressEvent(
            event_type=ProgressEventType.STARTED,
            source="chr1",
            current=0,
            total=1000
        )

        assert event.event_type == ProgressEventType.STARTED
        assert event.source == "chr1"
        assert event.current == 0
        assert event.total == 1000
        assert event.percentage == 0.0
        assert event.timestamp is not None

    def test_event_auto_percentage(self):
        """Test automatic percentage calculation."""
        event = ProgressEvent(
            event_type=ProgressEventType.UPDATED,
            source="chr2",
            current=500,
            total=1000
        )

        assert event.percentage == 50.0

    def test_event_with_metadata(self):
        """Test event with metadata."""
        metadata = {"worker_id": 1, "chromosome_length": 100000}
        event = ProgressEvent(
            event_type=ProgressEventType.UPDATED,
            source="chr3",
            metadata=metadata
        )

        assert event.metadata == metadata


class TestProgressObserver:
    """Test ProgressObserver abstract base class."""

    def test_observer_interface(self):
        """Test that observer defines required interface."""
        # Should not be able to instantiate abstract class
        with pytest.raises(TypeError):
            ProgressObserver()

    def test_concrete_observer(self):
        """Test concrete observer implementation."""
        class TestObserver(ProgressObserver):
            def __init__(self):
                self.events = []

            def update(self, event):
                self.events.append(event)

        observer = TestObserver()
        event = ProgressEvent(ProgressEventType.STARTED, "test")
        observer.update(event)

        assert len(observer.events) == 1
        assert observer.events[0] == event


class TestProgressSubject:
    """Test ProgressSubject implementation."""

    def test_subject_initialization(self):
        """Test subject initialization."""
        subject = ProgressSubject()
        assert subject._observers == []
        assert subject._event_history is None

        subject_with_history = ProgressSubject(track_history=True)
        assert subject_with_history._event_history == []

    def test_attach_detach_observer(self):
        """Test attaching and detaching observers."""
        subject = ProgressSubject()
        observer = Mock(spec=ProgressObserver)

        # Attach observer
        subject.attach(observer)
        assert observer in subject._observers

        # Attach same observer again (should not duplicate)
        subject.attach(observer)
        assert subject._observers.count(observer) == 1

        # Detach observer
        subject.detach(observer)
        assert observer not in subject._observers

    def test_notify_observers(self):
        """Test notifying observers."""
        subject = ProgressSubject()
        observer1 = Mock(spec=ProgressObserver)
        observer1.should_handle.return_value = True
        observer2 = Mock(spec=ProgressObserver)
        observer2.should_handle.return_value = True

        subject.attach(observer1)
        subject.attach(observer2)

        event = ProgressEvent(ProgressEventType.STARTED, "test")
        subject.notify(event)

        observer1.update.assert_called_once_with(event)
        observer2.update.assert_called_once_with(event)

    def test_selective_notification(self):
        """Test selective notification based on should_handle."""
        subject = ProgressSubject()
        observer1 = Mock(spec=ProgressObserver)
        observer1.should_handle.return_value = True
        observer2 = Mock(spec=ProgressObserver)
        observer2.should_handle.return_value = False

        subject.attach(observer1)
        subject.attach(observer2)

        event = ProgressEvent(ProgressEventType.UPDATED, "test")
        subject.notify(event)

        observer1.update.assert_called_once_with(event)
        observer2.update.assert_not_called()

    def test_convenience_methods(self):
        """Test convenience notification methods."""
        subject = ProgressSubject()
        observer = Mock(spec=ProgressObserver)
        observer.should_handle.return_value = True
        subject.attach(observer)

        # Test notify_start
        subject.notify_start("chr1", total=1000, message="Starting")
        assert observer.update.call_count == 1
        event = observer.update.call_args[0][0]
        assert event.event_type == ProgressEventType.STARTED
        assert event.source == "chr1"
        assert event.total == 1000
        assert event.message == "Starting"

        # Test notify_progress
        subject.notify_progress("chr1", current=500, total=1000)
        assert observer.update.call_count == 2
        event = observer.update.call_args[0][0]
        assert event.event_type == ProgressEventType.UPDATED
        assert event.current == 500

        # Test notify_complete
        subject.notify_complete("chr1", message="Done")
        assert observer.update.call_count == 3
        event = observer.update.call_args[0][0]
        assert event.event_type == ProgressEventType.COMPLETED

        # Test notify_failure
        subject.notify_failure("chr1", message="Error")
        assert observer.update.call_count == 4
        event = observer.update.call_args[0][0]
        assert event.event_type == ProgressEventType.FAILED

    def test_event_history(self):
        """Test event history tracking."""
        subject = ProgressSubject(track_history=True)

        subject.notify_start("chr1", total=100)
        subject.notify_progress("chr1", 50, 100)
        subject.notify_complete("chr1")

        history = subject.get_event_history()
        assert len(history) == 3
        assert history[0].event_type == ProgressEventType.STARTED
        assert history[1].event_type == ProgressEventType.UPDATED
        assert history[2].event_type == ProgressEventType.COMPLETED


class TestTerminalProgressObserver:
    """Test TerminalProgressObserver."""

    @patch('PyMaSC.core.observer.logger')
    def test_terminal_observer_messages(self, mock_logger):
        """Test terminal observer handles messages."""
        observer = TerminalProgressObserver(use_progress_bar=False)

        # Start event
        event = ProgressEvent(
            ProgressEventType.STARTED,
            "chr1",
            message="Starting chromosome 1"
        )
        observer.update(event)
        mock_logger.info.assert_called_with("chr1: Starting chromosome 1")

        # Complete event
        event = ProgressEvent(
            ProgressEventType.COMPLETED,
            "chr1",
            message="Completed chromosome 1"
        )
        observer.update(event)
        assert mock_logger.info.call_count == 2

        # Failure event
        event = ProgressEvent(
            ProgressEventType.FAILED,
            "chr2",
            message="Failed to process"
        )
        observer.update(event)
        mock_logger.error.assert_called_with("chr2: Failed to process")


class TestMultiprocessProgressObserver:
    """Test MultiprocessProgressObserver."""

    def test_multiprocess_observer_tracking(self):
        """Test multiprocess observer tracks multiple sources."""
        observer = MultiprocessProgressObserver()

        # Start multiple sources
        observer.update(ProgressEvent(ProgressEventType.STARTED, "chr1", total=100))
        observer.update(ProgressEvent(ProgressEventType.STARTED, "chr2", total=200))

        assert "chr1" in observer._active_sources
        assert "chr2" in observer._active_sources

        # Update progress
        observer.update(ProgressEvent(ProgressEventType.UPDATED, "chr1", current=50))
        assert observer._active_sources["chr1"]["current"] == 50

        # Complete one source
        observer.update(ProgressEvent(ProgressEventType.COMPLETED, "chr1"))
        assert "chr1" not in observer._active_sources
        assert "chr2" in observer._active_sources


class TestFileProgressObserver:
    """Test FileProgressObserver."""

    def test_file_observer_json_format(self):
        """Test file observer with JSON format."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            temp_path = f.name

        try:
            with FileProgressObserver(temp_path, format_type='json') as observer:
                event = ProgressEvent(
                    ProgressEventType.STARTED,
                    "chr1",
                    current=0,
                    total=100,
                    message="Starting"
                )
                observer.update(event)

            # Read and verify JSON
            with open(temp_path, 'r') as f:
                data = json.loads(f.readline())
                assert data['event_type'] == 'STARTED'
                assert data['source'] == 'chr1'
                assert data['total'] == 100
                assert data['message'] == 'Starting'
        finally:
            os.unlink(temp_path)

    def test_file_observer_text_format(self):
        """Test file observer with text format."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            temp_path = f.name

        try:
            with FileProgressObserver(temp_path, format_type='text') as observer:
                event = ProgressEvent(
                    ProgressEventType.UPDATED,
                    "chr2",
                    current=50,
                    total=100
                )
                observer.update(event)

            # Read and verify text
            with open(temp_path, 'r') as f:
                line = f.readline()
                assert 'UPDATED' in line
                assert 'chr2' in line
                assert '50.0%' in line
        finally:
            os.unlink(temp_path)


class TestAggregateProgressObserver:
    """Test AggregateProgressObserver."""

    def test_aggregate_observer(self):
        """Test aggregate observer combines progress."""
        observer = AggregateProgressObserver(expected_sources=["chr1", "chr2", "chr3"])

        # Start sources
        observer.update(ProgressEvent(ProgressEventType.STARTED, "chr1", total=100))
        observer.update(ProgressEvent(ProgressEventType.STARTED, "chr2", total=200))

        # Update progress
        observer.update(ProgressEvent(ProgressEventType.UPDATED, "chr1", current=50))
        observer.update(ProgressEvent(ProgressEventType.UPDATED, "chr2", current=100))

        # Complete one, fail another
        observer.update(ProgressEvent(ProgressEventType.COMPLETED, "chr1"))
        observer.update(ProgressEvent(ProgressEventType.FAILED, "chr2"))

        summary = observer.get_summary()
        assert summary['total_sources'] == 2
        assert summary['completed'] == 1
        assert summary['failed'] == 1
        assert summary['in_progress'] == 0


class TestProgressAdapter:
    """Test progress adapter functionality."""

    def test_progress_bar_adapter(self):
        """Test ProgressBarAdapter observer notifications."""
        # Create adapter (don't worry about actual progress bar display)
        adapter = ProgressBarAdapter()

        # Track events
        events = []

        class TestObserver(ProgressObserver):
            def update(self, event):
                events.append(event)

        # Attach observer
        observer = TestObserver()
        adapter._subject.attach(observer)

        # Set should trigger start event
        adapter.set("chr1", 1000)
        assert len(events) == 1
        assert events[0].event_type == ProgressEventType.STARTED
        assert events[0].source == "chr1"
        assert events[0].total == 1000

        # Update should trigger progress event
        adapter.update(500)
        assert len(events) == 2
        assert events[1].event_type == ProgressEventType.UPDATED
        assert events[1].current == 500

        # Clean should trigger complete event
        adapter.clean()
        assert len(events) == 3
        assert events[2].event_type == ProgressEventType.COMPLETED

    def test_legacy_progress_observer(self):
        """Test LegacyProgressObserver."""
        mock_progress_bar = Mock()
        observer = LegacyProgressObserver(mock_progress_bar)

        # Start event
        event = ProgressEvent(ProgressEventType.STARTED, "chr1", total=100)
        observer.update(event)
        mock_progress_bar.set.assert_called_once_with("chr1", 100)

        # Update event
        event = ProgressEvent(ProgressEventType.UPDATED, "chr1", current=50)
        observer.update(event)
        mock_progress_bar.update.assert_called_once_with(50)

        # Complete event
        event = ProgressEvent(ProgressEventType.COMPLETED, "chr1")
        observer.update(event)
        mock_progress_bar.clean.assert_called_once()


class TestProgressManager:
    """Test ProgressManager functionality."""

    def test_progress_manager_creation(self):
        """Test progress manager creation."""
        manager = ProgressManager()
        assert manager.subject is not None
        assert manager._observers == []
        assert manager._legacy_components == []

    def test_create_progress_bar(self):
        """Test creating progress bars."""
        manager = ProgressManager()

        # Observer-aware progress bar
        pb_aware = manager.create_progress_bar(observer_aware=True)
        assert isinstance(pb_aware, ProgressBarAdapter)

        # Legacy progress bar
        pb_legacy = manager.create_progress_bar(observer_aware=False)
        from PyMaSC.utils.progress import ProgressBar
        assert isinstance(pb_legacy, ProgressBar)

    def test_add_observers(self):
        """Test adding different observers."""
        manager = ProgressManager()

        # Add terminal observer
        terminal_obs = manager.add_terminal_observer()
        assert terminal_obs in manager._observers

        # Add multiprocess observer
        mp_obs = manager.add_multiprocess_observer()
        assert mp_obs in manager._observers
        assert len(manager._observers) == 2

    def test_manager_notifications(self):
        """Test manager notification methods."""
        manager = ProgressManager()
        observer = Mock(spec=ProgressObserver)
        observer.should_handle.return_value = True
        manager.subject.attach(observer)

        # Test notifications
        manager.notify_start("test", total=100)
        manager.notify_progress("test", 50, 100)
        manager.notify_complete("test")

        assert observer.update.call_count == 3

    def test_global_progress_manager(self):
        """Test global progress manager."""
        reset_progress_manager()

        manager1 = get_progress_manager()
        manager2 = get_progress_manager()

        assert manager1 is manager2  # Should be same instance

        reset_progress_manager()
        manager3 = get_progress_manager()

        assert manager3 is not manager1  # Should be new instance