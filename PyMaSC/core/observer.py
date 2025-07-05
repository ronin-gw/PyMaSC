"""Observer pattern implementation for event-driven progress reporting.

This module provides a flexible and extensible progress reporting system
using the Observer pattern. It decouples progress generation from display,
allowing multiple observers to track progress in different ways.

Key components:
- ProgressEvent: Different types of progress events
- ProgressObserver: Abstract base class for observers
- ProgressSubject: Subject that manages observers and notifications
- Concrete observers for terminal, file, and multiprocessing
"""
from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum, auto
from typing import Optional, List, Dict, Any, Union
import logging
from datetime import datetime
import threading
from collections import defaultdict

logger = logging.getLogger(__name__)


class ProgressEventType(Enum):
    """Types of progress events."""
    STARTED = auto()
    UPDATED = auto()
    COMPLETED = auto()
    FAILED = auto()
    STAGE_CHANGED = auto()
    MESSAGE = auto()


@dataclass
class ProgressEvent:
    """Progress event data structure.

    Attributes:
        event_type: Type of progress event
        source: Source identifier (e.g., chromosome name)
        current: Current progress value
        total: Total expected value
        percentage: Progress percentage (0-100)
        message: Optional message for the event
        timestamp: Event timestamp
        metadata: Additional event metadata
    """
    event_type: ProgressEventType
    source: str
    current: Optional[int] = None
    total: Optional[int] = None
    percentage: Optional[float] = None
    message: Optional[str] = None
    timestamp: datetime = None
    metadata: Dict[str, Any] = None

    def __post_init__(self):
        """Initialize timestamp and calculate percentage if needed."""
        if self.timestamp is None:
            self.timestamp = datetime.now()

        if self.percentage is None and self.current is not None and self.total:
            self.percentage = (self.current / self.total) * 100


class ProgressObserver(ABC):
    """Abstract base class for progress observers.

    Defines the interface that all progress observers must implement.
    Observers receive notifications about progress events and handle
    them according to their specific implementation.
    """

    @abstractmethod
    def update(self, event: ProgressEvent) -> None:
        """Handle a progress event.

        Args:
            event: Progress event to handle
        """
        pass

    def should_handle(self, event: ProgressEvent) -> bool:
        """Check if this observer should handle the event.

        Args:
            event: Progress event to check

        Returns:
            True if the observer should handle this event
        """
        return True


class ProgressSubject:
    """Subject that manages progress observers.

    This class implements the subject side of the Observer pattern.
    It manages a list of observers and notifies them about progress
    events. Thread-safe for concurrent operations.

    Attributes:
        observers: List of registered observers
        _lock: Thread lock for observer list modifications
        _event_history: Optional event history tracking
    """

    def __init__(self, track_history: bool = False):
        """Initialize progress subject.

        Args:
            track_history: Whether to keep event history
        """
        self._observers: List[ProgressObserver] = []
        self._lock = threading.Lock()
        self._track_history = track_history
        self._event_history: List[ProgressEvent] = [] if track_history else None

    def attach(self, observer: ProgressObserver) -> None:
        """Attach an observer to receive notifications.

        Args:
            observer: Observer to attach
        """
        with self._lock:
            if observer not in self._observers:
                self._observers.append(observer)
                logger.debug(f"Attached observer: {observer.__class__.__name__}")

    def detach(self, observer: ProgressObserver) -> None:
        """Detach an observer from receiving notifications.

        Args:
            observer: Observer to detach
        """
        with self._lock:
            if observer in self._observers:
                self._observers.remove(observer)
                logger.debug(f"Detached observer: {observer.__class__.__name__}")

    def notify(self, event: ProgressEvent) -> None:
        """Notify all observers about a progress event.

        Args:
            event: Progress event to broadcast
        """
        if self._track_history and self._event_history is not None:
            self._event_history.append(event)

        with self._lock:
            observers = self._observers.copy()

        for observer in observers:
            try:
                if observer.should_handle(event):
                    observer.update(event)
            except Exception as e:
                logger.error(f"Observer {observer.__class__.__name__} "
                           f"failed to handle event: {e}")

    def notify_start(self, source: str, total: Optional[int] = None,
                     message: Optional[str] = None, **metadata) -> None:
        """Convenience method to notify start of operation.

        Args:
            source: Source identifier
            total: Total expected progress value
            message: Optional start message
            **metadata: Additional metadata
        """
        event = ProgressEvent(
            event_type=ProgressEventType.STARTED,
            source=source,
            current=0,
            total=total,
            message=message,
            metadata=metadata
        )
        self.notify(event)

    def notify_progress(self, source: str, current: int, total: Optional[int] = None,
                       message: Optional[str] = None, **metadata) -> None:
        """Convenience method to notify progress update.

        Args:
            source: Source identifier
            current: Current progress value
            total: Total expected value
            message: Optional progress message
            **metadata: Additional metadata
        """
        event = ProgressEvent(
            event_type=ProgressEventType.UPDATED,
            source=source,
            current=current,
            total=total,
            message=message,
            metadata=metadata
        )
        self.notify(event)

    def notify_complete(self, source: str, message: Optional[str] = None,
                       **metadata) -> None:
        """Convenience method to notify completion.

        Args:
            source: Source identifier
            message: Optional completion message
            **metadata: Additional metadata
        """
        event = ProgressEvent(
            event_type=ProgressEventType.COMPLETED,
            source=source,
            message=message,
            metadata=metadata
        )
        self.notify(event)

    def notify_failure(self, source: str, message: Optional[str] = None,
                      **metadata) -> None:
        """Convenience method to notify failure.

        Args:
            source: Source identifier
            message: Optional failure message
            **metadata: Additional metadata
        """
        event = ProgressEvent(
            event_type=ProgressEventType.FAILED,
            source=source,
            message=message,
            metadata=metadata
        )
        self.notify(event)

    def get_event_history(self) -> Optional[List[ProgressEvent]]:
        """Get event history if tracking is enabled.

        Returns:
            List of events or None if not tracking
        """
        if self._track_history:
            with self._lock:
                return self._event_history.copy()
        return None


class TerminalProgressObserver(ProgressObserver):
    """Terminal-based progress display observer.

    Displays progress information in the terminal using the
    existing ProgressBar functionality. Handles single-source
    progress tracking with visual progress bar.
    """

    def __init__(self, use_progress_bar: bool = True):
        """Initialize terminal progress observer.

        Args:
            use_progress_bar: Whether to use visual progress bar
        """
        self.use_progress_bar = use_progress_bar
        self._active_sources: Dict[str, Any] = {}

        if use_progress_bar:
            from PyMaSC.utils.progress import ProgressBar
            self._progress_bar = ProgressBar()
        else:
            self._progress_bar = None

    def update(self, event: ProgressEvent) -> None:
        """Handle progress event with terminal display.

        Args:
            event: Progress event to display
        """
        if event.event_type == ProgressEventType.STARTED:
            self._handle_start(event)
        elif event.event_type == ProgressEventType.UPDATED:
            self._handle_update(event)
        elif event.event_type == ProgressEventType.COMPLETED:
            self._handle_complete(event)
        elif event.event_type == ProgressEventType.FAILED:
            self._handle_failure(event)
        elif event.event_type == ProgressEventType.MESSAGE:
            self._handle_message(event)

    def _handle_start(self, event: ProgressEvent) -> None:
        """Handle start event."""
        self._active_sources[event.source] = {
            'total': event.total,
            'start_time': event.timestamp
        }

        if self._progress_bar and event.total:
            self._progress_bar.set(event.source, event.total)

        if event.message:
            logger.info(f"{event.source}: {event.message}")

    def _handle_update(self, event: ProgressEvent) -> None:
        """Handle progress update."""
        if self._progress_bar and event.current is not None:
            self._progress_bar.update(event.current)

    def _handle_complete(self, event: ProgressEvent) -> None:
        """Handle completion event."""
        if event.source in self._active_sources:
            del self._active_sources[event.source]

        if self._progress_bar:
            self._progress_bar.clean()

        if event.message:
            logger.info(f"{event.source}: {event.message}")

    def _handle_failure(self, event: ProgressEvent) -> None:
        """Handle failure event."""
        if event.source in self._active_sources:
            del self._active_sources[event.source]

        if self._progress_bar:
            self._progress_bar.clean()

        logger.error(f"{event.source}: {event.message or 'Failed'}")

    def _handle_message(self, event: ProgressEvent) -> None:
        """Handle message event."""
        if event.message:
            logger.info(f"{event.source}: {event.message}")


class MultiprocessProgressObserver(ProgressObserver):
    """Multiprocess-aware progress observer.

    Handles progress events from multiple sources (processes/threads)
    simultaneously. Uses the MultiLineProgressManager for concurrent
    progress display.
    """

    def __init__(self, queue=None):
        """Initialize multiprocess progress observer.

        Args:
            queue: Optional queue for inter-process communication
        """
        from PyMaSC.utils.progress import MultiLineProgressManager
        self._manager = MultiLineProgressManager()
        self._queue = queue
        self._active_sources: Dict[str, Dict[str, Any]] = defaultdict(dict)

    def update(self, event: ProgressEvent) -> None:
        """Handle progress event from multiple sources.

        Args:
            event: Progress event to handle
        """
        if event.event_type == ProgressEventType.STARTED:
            self._handle_start(event)
        elif event.event_type == ProgressEventType.UPDATED:
            self._handle_update(event)
        elif event.event_type == ProgressEventType.COMPLETED:
            self._handle_complete(event)
        elif event.event_type == ProgressEventType.FAILED:
            self._handle_failure(event)

    def _handle_start(self, event: ProgressEvent) -> None:
        """Handle start event for a source."""
        self._active_sources[event.source] = {
            'total': event.total,
            'current': 0,
            'start_time': event.timestamp
        }
        self._update_display(event.source)

    def _handle_update(self, event: ProgressEvent) -> None:
        """Handle progress update for a source."""
        if event.source in self._active_sources:
            self._active_sources[event.source]['current'] = event.current
            self._update_display(event.source)

    def _handle_complete(self, event: ProgressEvent) -> None:
        """Handle completion for a source."""
        if event.source in self._active_sources:
            self._manager.erase(event.source)
            del self._active_sources[event.source]

    def _handle_failure(self, event: ProgressEvent) -> None:
        """Handle failure for a source."""
        if event.source in self._active_sources:
            self._manager.erase(event.source)
            del self._active_sources[event.source]

    def _update_display(self, source: str) -> None:
        """Update display for a specific source."""
        info = self._active_sources[source]
        if info.get('total'):
            percentage = (info['current'] / info['total']) * 100
            bar_length = 50
            filled = int(bar_length * percentage / 100)
            bar = '=' * filled + '-' * (bar_length - filled)
            display = f"[{bar}] {percentage:.1f}%"
        else:
            display = f"Processing... {info['current']} items"

        self._manager.update(source, display)


class FileProgressObserver(ProgressObserver):
    """File-based progress logging observer.

    Logs progress events to a file for persistent tracking
    and post-analysis. Useful for debugging and monitoring
    long-running operations.
    """

    def __init__(self, log_path: str, format_type: str = 'json'):
        """Initialize file progress observer.

        Args:
            log_path: Path to log file
            format_type: Log format ('json' or 'text')
        """
        self.log_path = log_path
        self.format_type = format_type
        self._file = None
        self._lock = threading.Lock()

    def __enter__(self):
        """Context manager entry."""
        self._file = open(self.log_path, 'a')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        if self._file:
            self._file.close()

    def update(self, event: ProgressEvent) -> None:
        """Log progress event to file.

        Args:
            event: Progress event to log
        """
        with self._lock:
            if self.format_type == 'json':
                self._write_json(event)
            else:
                self._write_text(event)

    def _write_json(self, event: ProgressEvent) -> None:
        """Write event in JSON format."""
        import json

        data = {
            'timestamp': event.timestamp.isoformat(),
            'event_type': event.event_type.name,
            'source': event.source,
            'current': event.current,
            'total': event.total,
            'percentage': event.percentage,
            'message': event.message,
            'metadata': event.metadata
        }

        if self._file:
            json.dump(data, self._file)
            self._file.write('\n')
            self._file.flush()

    def _write_text(self, event: ProgressEvent) -> None:
        """Write event in text format."""
        parts = [
            event.timestamp.strftime('%Y-%m-%d %H:%M:%S'),
            event.event_type.name,
            event.source
        ]

        if event.percentage is not None:
            parts.append(f"{event.percentage:.1f}%")
        elif event.current is not None:
            parts.append(f"{event.current}/{event.total or '?'}")

        if event.message:
            parts.append(event.message)

        if self._file:
            self._file.write(' | '.join(parts) + '\n')
            self._file.flush()


class AggregateProgressObserver(ProgressObserver):
    """Aggregate progress from multiple sources.

    Combines progress from multiple sources into a single
    overall progress indicator. Useful for tracking overall
    job progress across multiple chromosomes or files.
    """

    def __init__(self, expected_sources: Optional[List[str]] = None):
        """Initialize aggregate progress observer.

        Args:
            expected_sources: List of expected source identifiers
        """
        self.expected_sources = set(expected_sources) if expected_sources else None
        self._source_progress: Dict[str, Dict[str, Any]] = {}
        self._completed_sources: set = set()
        self._failed_sources: set = set()
        self._lock = threading.Lock()

    def update(self, event: ProgressEvent) -> None:
        """Update aggregate progress based on event.

        Args:
            event: Progress event to aggregate
        """
        with self._lock:
            if event.event_type == ProgressEventType.STARTED:
                self._source_progress[event.source] = {
                    'current': 0,
                    'total': event.total
                }
            elif event.event_type == ProgressEventType.UPDATED:
                if event.source in self._source_progress:
                    self._source_progress[event.source]['current'] = event.current
            elif event.event_type == ProgressEventType.COMPLETED:
                self._completed_sources.add(event.source)
                if event.source in self._source_progress:
                    self._source_progress[event.source]['current'] = \
                        self._source_progress[event.source].get('total', 0)
            elif event.event_type == ProgressEventType.FAILED:
                self._failed_sources.add(event.source)

        # Log aggregate progress
        self._log_aggregate_progress()

    def _log_aggregate_progress(self) -> None:
        """Calculate and log aggregate progress."""
        total_current = sum(p['current'] for p in self._source_progress.values())
        total_expected = sum(p['total'] or 0 for p in self._source_progress.values())

        if total_expected > 0:
            percentage = (total_current / total_expected) * 100
            logger.info(f"Overall progress: {percentage:.1f}% "
                       f"({len(self._completed_sources)} completed, "
                       f"{len(self._failed_sources)} failed)")

    def get_summary(self) -> Dict[str, Any]:
        """Get summary of aggregate progress.

        Returns:
            Dictionary with progress summary
        """
        with self._lock:
            total_sources = len(self._source_progress)
            return {
                'total_sources': total_sources,
                'completed': len(self._completed_sources),
                'failed': len(self._failed_sources),
                'in_progress': total_sources - len(self._completed_sources) - len(self._failed_sources),
                'completion_rate': len(self._completed_sources) / total_sources if total_sources > 0 else 0
            }