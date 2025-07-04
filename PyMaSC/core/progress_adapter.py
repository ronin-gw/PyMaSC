"""Adapter to bridge new observer-based progress with legacy progress system.

This module provides adapters and utilities to integrate the new observer-based
progress reporting system with the existing progress bar implementation. It
ensures backward compatibility while enabling gradual migration.
"""
import logging
from typing import Optional, Dict, Any, Union
from collections import defaultdict

from .observer import (
    ProgressObserver, ProgressSubject, ProgressEvent,
    ProgressEventType, TerminalProgressObserver,
    MultiprocessProgressObserver
)
from PyMaSC.utils.progress import ProgressBar, ProgressHook, MultiLineProgressManager

logger = logging.getLogger(__name__)


class ProgressBarAdapter(ProgressBar):
    """Adapter that makes ProgressBar work as a ProgressSubject.
    
    This adapter wraps the existing ProgressBar class and adds
    observer pattern capabilities, allowing it to notify observers
    about progress events while maintaining the original API.
    """
    
    def __init__(self, *args, **kwargs):
        """Initialize progress bar adapter."""
        # Store subject before calling super
        self._subject = ProgressSubject()
        self._current_source = None
        self._current_total = None
        self._current_value = 0
        
        # Initialize parent
        super().__init__(*args, **kwargs)
        
        # Store references to actual implementation methods
        self._do_update = self._update if hasattr(self, '_update') else lambda x: None
        self._do_clean = self._clean if hasattr(self, '_clean') else lambda: None
        self._do_format = self._format if hasattr(self, '_format') else lambda x: None
        
        # Force restore our methods that may have been overridden
        # Save the class method references
        self.update = self.__class__.update.__get__(self, self.__class__)
        self.clean = self.__class__.clean.__get__(self, self.__class__)
        self.set = self.__class__.set.__get__(self, self.__class__)
    
    def attach_observer(self, observer: ProgressObserver) -> None:
        """Attach an observer to this progress bar.
        
        Args:
            observer: Observer to attach
        """
        self._subject.attach(observer)
    
    def detach_observer(self, observer: ProgressObserver) -> None:
        """Detach an observer from this progress bar.
        
        Args:
            observer: Observer to detach
        """
        self._subject.detach(observer)
    
    def set(self, name: str, maxval: int) -> None:
        """Set up progress bar and notify start event.
        
        Args:
            name: Progress source name
            maxval: Maximum value for progress
        """
        super().set(name, maxval)
        self._current_source = name
        self._current_total = maxval
        
        # Notify observers about start
        self._subject.notify_start(
            source=name,
            total=maxval,
            message=f"Starting {name}"
        )
    
    def update(self, val: int) -> None:
        """Update progress and notify observers.
        
        Args:
            val: Current progress value
        """
        # Store current value
        self._current_value = val
        
        # Update the actual progress bar display if enabled
        if self.global_switch and self._do_update:
            try:
                self._do_update(val)
            except:
                pass  # Ignore display errors
        
        # Always notify observers about progress
        if self._current_source:
            self._subject.notify_progress(
                source=self._current_source,
                current=val,
                total=self._current_total
            )
    
    def clean(self) -> None:
        """Clean up progress bar and notify completion."""
        # Clean the actual progress bar display if enabled
        if self.global_switch and self._do_clean:
            try:
                self._do_clean()
            except:
                pass  # Ignore display errors
        
        # Always notify observers about completion
        if self._current_source:
            self._subject.notify_complete(
                source=self._current_source,
                message=f"Completed {self._current_source}"
            )
            self._current_source = None
            self._current_total = None
            self._current_value = 0


class ProgressHookAdapter(ProgressHook):
    """Adapter for ProgressHook with observer pattern support.
    
    Extends ProgressHook to support the observer pattern while
    maintaining compatibility with multiprocessing queue communication.
    """
    
    def __init__(self, queue, *args, **kwargs):
        """Initialize progress hook adapter."""
        super().__init__(queue, *args, **kwargs)
        self._subject = ProgressSubject()
    
    def attach_observer(self, observer: ProgressObserver) -> None:
        """Attach observer to this progress hook."""
        self._subject.attach(observer)
    
    def format(self, s: str) -> None:
        """Format and send progress through queue and to observers."""
        super().format(s)
        
        # Also notify local observers if any
        if self.name and hasattr(self, '_unit'):
            current = len(s) * self._unit if hasattr(self, '_unit') else 0
            self._subject.notify_progress(
                source=self.name,
                current=int(current),
                total=getattr(self, '_maxval', None)
            )


class LegacyProgressObserver(ProgressObserver):
    """Observer that delegates to legacy progress components.
    
    This observer receives events from the new observer system
    and delegates them to legacy progress bar components.
    """
    
    def __init__(self, legacy_progress: Union[ProgressBar, MultiLineProgressManager]):
        """Initialize legacy progress observer.
        
        Args:
            legacy_progress: Legacy progress component to delegate to
        """
        self.legacy_progress = legacy_progress
        self._source_info: Dict[str, Dict[str, Any]] = {}
    
    def update(self, event: ProgressEvent) -> None:
        """Handle progress event by delegating to legacy component.
        
        Args:
            event: Progress event to handle
        """
        # Check if it's a ProgressBar-like object (has set/update/clean methods)
        if hasattr(self.legacy_progress, 'set') and hasattr(self.legacy_progress, 'update'):
            self._handle_progress_bar_event(event)
        # Check if it's a MultiLineProgressManager-like object (has update/erase methods)
        elif hasattr(self.legacy_progress, 'update') and hasattr(self.legacy_progress, 'erase'):
            self._handle_multiline_event(event)
    
    def _handle_progress_bar_event(self, event: ProgressEvent) -> None:
        """Handle event for ProgressBar."""
        if event.event_type == ProgressEventType.STARTED and event.total:
            self.legacy_progress.set(event.source, event.total)
            self._source_info[event.source] = {'total': event.total}
        
        elif event.event_type == ProgressEventType.UPDATED and event.current is not None:
            self.legacy_progress.update(event.current)
        
        elif event.event_type in (ProgressEventType.COMPLETED, ProgressEventType.FAILED):
            self.legacy_progress.clean()
            if event.source in self._source_info:
                del self._source_info[event.source]
    
    def _handle_multiline_event(self, event: ProgressEvent) -> None:
        """Handle event for MultiLineProgressManager."""
        if event.event_type == ProgressEventType.STARTED:
            self._source_info[event.source] = {
                'total': event.total,
                'current': 0
            }
            self._update_multiline_display(event.source)
        
        elif event.event_type == ProgressEventType.UPDATED:
            if event.source in self._source_info:
                self._source_info[event.source]['current'] = event.current
                self._update_multiline_display(event.source)
        
        elif event.event_type in (ProgressEventType.COMPLETED, ProgressEventType.FAILED):
            self.legacy_progress.erase(event.source)
            if event.source in self._source_info:
                del self._source_info[event.source]
    
    def _update_multiline_display(self, source: str) -> None:
        """Update multiline display for a source."""
        info = self._source_info.get(source, {})
        current = info.get('current', 0)
        total = info.get('total', 0)
        
        if total > 0:
            percentage = (current / total) * 100
            display = f"[{'=' * int(percentage / 2):50s}] {percentage:.1f}%"
        else:
            display = f"Processing... {current} items"
        
        self.legacy_progress.update(source, display)


class ProgressManager:
    """Central progress management with observer pattern.
    
    This class provides a unified interface for progress management,
    supporting both the new observer pattern and legacy progress
    components through adapters.
    """
    
    def __init__(self, use_legacy: bool = True):
        """Initialize progress manager.
        
        Args:
            use_legacy: Whether to use legacy progress components
        """
        self.subject = ProgressSubject(track_history=True)
        self.use_legacy = use_legacy
        self._observers = []
        self._legacy_components = []
    
    def create_progress_bar(self, observer_aware: bool = True) -> Union[ProgressBar, ProgressBarAdapter]:
        """Create a progress bar with optional observer support.
        
        Args:
            observer_aware: Whether to create observer-aware adapter
            
        Returns:
            ProgressBar or ProgressBarAdapter instance
        """
        if observer_aware:
            progress_bar = ProgressBarAdapter()
            # Connect to central subject
            progress_bar.attach_observer(self.subject)
            return progress_bar
        else:
            return ProgressBar()
    
    def add_terminal_observer(self) -> TerminalProgressObserver:
        """Add a terminal progress observer.
        
        Returns:
            The created observer
        """
        observer = TerminalProgressObserver()
        self.subject.attach(observer)
        self._observers.append(observer)
        return observer
    
    def add_multiprocess_observer(self, queue=None) -> MultiprocessProgressObserver:
        """Add a multiprocess progress observer.
        
        Args:
            queue: Optional queue for IPC
            
        Returns:
            The created observer
        """
        observer = MultiprocessProgressObserver(queue)
        self.subject.attach(observer)
        self._observers.append(observer)
        return observer
    
    def bridge_legacy_component(self, legacy_component: Union[ProgressBar, MultiLineProgressManager]) -> None:
        """Bridge a legacy progress component with observer pattern.
        
        Args:
            legacy_component: Legacy component to bridge
        """
        observer = LegacyProgressObserver(legacy_component)
        self.subject.attach(observer)
        self._legacy_components.append(legacy_component)
    
    def notify_start(self, source: str, total: Optional[int] = None, **kwargs) -> None:
        """Notify start of operation.
        
        Args:
            source: Source identifier
            total: Total expected value
            **kwargs: Additional metadata
        """
        self.subject.notify_start(source, total, **kwargs)
    
    def notify_progress(self, source: str, current: int, total: Optional[int] = None, **kwargs) -> None:
        """Notify progress update.
        
        Args:
            source: Source identifier
            current: Current value
            total: Total expected value
            **kwargs: Additional metadata
        """
        self.subject.notify_progress(source, current, total, **kwargs)
    
    def notify_complete(self, source: str, **kwargs) -> None:
        """Notify completion.
        
        Args:
            source: Source identifier
            **kwargs: Additional metadata
        """
        self.subject.notify_complete(source, **kwargs)
    
    def notify_failure(self, source: str, **kwargs) -> None:
        """Notify failure.
        
        Args:
            source: Source identifier
            **kwargs: Additional metadata
        """
        self.subject.notify_failure(source, **kwargs)
    
    def get_summary(self) -> Dict[str, Any]:
        """Get progress summary.
        
        Returns:
            Summary of all progress events
        """
        history = self.subject.get_event_history() or []
        
        sources = defaultdict(lambda: {'events': 0, 'last_event': None})
        for event in history:
            sources[event.source]['events'] += 1
            sources[event.source]['last_event'] = event.event_type.name
        
        return {
            'total_events': len(history),
            'sources': dict(sources),
            'observers': len(self._observers),
            'legacy_components': len(self._legacy_components)
        }


class ReadCountProgressBarAdapter:
    """Adapter for ReadCountProgressBar with observer support.
    
    This adapter wraps ReadCountProgressBar to add observer pattern
    capabilities while maintaining the original API.
    """
    
    def __init__(self, *args, **kwargs):
        """Initialize adapter."""
        from PyMaSC.utils.progress import ReadCountProgressBar
        self._bar = ReadCountProgressBar(*args, **kwargs)
        self._subject = ProgressSubject()
        self._current_chrom = None
        self._current_genome_total = None
        self._current_chrom_total = None
        self._genome_progress = 0
    
    def attach_observer(self, observer: ProgressObserver) -> None:
        """Attach an observer."""
        self._subject.attach(observer)
    
    def detach_observer(self, observer: ProgressObserver) -> None:
        """Detach an observer."""
        self._subject.detach(observer)
    
    def set_genome(self, maxval: int) -> None:
        """Set genome size and notify observers."""
        self._bar.set_genome(maxval)
        self._current_genome_total = maxval
        
        # Notify about genome-level start
        self._subject.notify_start(
            source="genome",
            total=maxval,
            message="Starting genome-wide read counting"
        )
    
    def set_chrom(self, maxval: int, name: str) -> None:
        """Set chromosome and notify observers."""
        # Complete previous chromosome if any
        if self._current_chrom:
            self._subject.notify_complete(
                source=f"chromosome:{self._current_chrom}",
                message=f"Completed {self._current_chrom}"
            )
        
        self._bar.set_chrom(maxval, name)
        self._current_chrom = name
        self._current_chrom_total = maxval
        
        # Notify about chromosome start
        self._subject.notify_start(
            source=f"chromosome:{name}",
            total=maxval,
            message=f"Processing {name}"
        )
    
    def update(self, val: int) -> None:
        """Update progress and notify observers."""
        self._bar.update(val)
        
        # Update chromosome progress
        if self._current_chrom:
            self._subject.notify_progress(
                source=f"chromosome:{self._current_chrom}",
                current=val,
                total=self._current_chrom_total
            )
        
        # Update genome progress
        self._genome_progress += 1
        if self._current_genome_total:
            self._subject.notify_progress(
                source="genome",
                current=self._genome_progress,
                total=self._current_genome_total
            )
    
    def finish(self) -> None:
        """Finish and notify completion."""
        self._bar.finish()
        
        # Complete current chromosome
        if self._current_chrom:
            self._subject.notify_complete(
                source=f"chromosome:{self._current_chrom}",
                message=f"Completed {self._current_chrom}"
            )
        
        # Complete genome
        self._subject.notify_complete(
            source="genome",
            message="Completed genome-wide read counting"
        )
    
    def disable_bar(self) -> None:
        """Disable visual progress bar but keep observer notifications."""
        self._bar.disable_bar()
    
    def __getattr__(self, name):
        """Delegate unknown attributes to wrapped bar."""
        return getattr(self._bar, name)


# Global progress manager instance
_global_progress_manager = None


def get_progress_manager() -> ProgressManager:
    """Get or create global progress manager.
    
    Returns:
        Global ProgressManager instance
    """
    global _global_progress_manager
    if _global_progress_manager is None:
        _global_progress_manager = ProgressManager()
    return _global_progress_manager


def reset_progress_manager() -> None:
    """Reset global progress manager."""
    global _global_progress_manager
    _global_progress_manager = None