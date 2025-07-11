"""Progress coordination for single and multi-process execution.

This module provides a unified interface for progress reporting
across different execution modes, abstracting the complexity of
progress management from the main calculation logic.
"""
from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum
from typing import Optional, Dict, Any, List, TYPE_CHECKING

if TYPE_CHECKING:
    from PyMaSC.core.observer import ProgressObserver
    from multiprocessing import Queue, Lock

from PyMaSC.utils.progress import ProgressBar, ProgressHook, MultiLineProgressManager
from PyMaSC.core.progress_adapter import ProgressBarAdapter, get_progress_manager

logger = logging.getLogger(__name__)


class ProgressMode(Enum):
    """Progress reporting mode."""
    SINGLE_LINE = "single_line"
    MULTI_LINE = "multi_line"
    OBSERVER = "observer"
    DISABLED = "disabled"


@dataclass
class ProgressConfig:
    """Configuration for progress reporting.
    
    Attributes:
        mode: Progress reporting mode
        enable_global: Whether to enable global progress switches
        observers: List of progress observers (for observer mode)
    """
    mode: ProgressMode = ProgressMode.SINGLE_LINE
    enable_global: bool = True
    observers: Optional[List['ProgressObserver']] = None


class ProgressReporter(ABC):
    """Abstract base for progress reporters."""
    
    @abstractmethod
    def start_chromosome(self, chrom: str, total: int) -> None:
        """Start progress for a chromosome."""
        pass
    
    @abstractmethod
    def update(self, chrom: str, position: int) -> None:
        """Update progress position."""
        pass
    
    @abstractmethod
    def finish_chromosome(self, chrom: str) -> None:
        """Finish progress for a chromosome."""
        pass
    
    @abstractmethod
    def cleanup(self) -> None:
        """Clean up progress display."""
        pass


class SingleLineProgressReporter(ProgressReporter):
    """Single-line progress reporter for single process mode."""
    
    def __init__(self):
        self._progress = ProgressBar()
        self._current_chrom: Optional[str] = None
    
    def start_chromosome(self, chrom: str, total: int) -> None:
        if self._current_chrom:
            self._progress.clean()
        self._current_chrom = chrom
        self._progress.set(chrom, total)
    
    def update(self, chrom: str, position: int) -> None:
        if chrom == self._current_chrom:
            self._progress.update(position)
    
    def finish_chromosome(self, chrom: str) -> None:
        if chrom == self._current_chrom:
            self._progress.clean()
            self._current_chrom = None
    
    def cleanup(self) -> None:
        self._progress.clean()


class MultiLineProgressReporter(ProgressReporter):
    """Multi-line progress reporter for multi-process mode."""
    
    def __init__(self, logger_lock: Optional['Lock'] = None):
        self._progress = MultiLineProgressManager()
        self._logger_lock = logger_lock
    
    def start_chromosome(self, chrom: str, total: int) -> None:
        # MultiLine manager handles this automatically
        pass
    
    def update(self, chrom: str, position: int) -> None:
        if self._logger_lock:
            with self._logger_lock:
                self._progress.update(chrom, position)
        else:
            self._progress.update(chrom, position)
    
    def finish_chromosome(self, chrom: str) -> None:
        if self._logger_lock:
            with self._logger_lock:
                self._progress.erase(chrom)
        else:
            self._progress.erase(chrom)
    
    def cleanup(self) -> None:
        self._progress.clean()


class ObserverProgressReporter(ProgressReporter):
    """Observer-based progress reporter."""
    
    def __init__(self, observers: Optional[List['ProgressObserver']] = None):
        self._adapter = ProgressBarAdapter()
        self._observers = observers or []
        self._current_chrom: Optional[str] = None
        
        # Attach observers
        for observer in self._observers:
            self._adapter._subject.attach(observer)
    
    def start_chromosome(self, chrom: str, total: int) -> None:
        if self._current_chrom:
            self._adapter.clean()
        self._current_chrom = chrom
        self._adapter.set(chrom, total)
    
    def update(self, chrom: str, position: int) -> None:
        if chrom == self._current_chrom:
            self._adapter.update(position)
    
    def finish_chromosome(self, chrom: str) -> None:
        if chrom == self._current_chrom:
            self._adapter.clean()
            self._current_chrom = None
    
    def cleanup(self) -> None:
        self._adapter.clean()


class NullProgressReporter(ProgressReporter):
    """No-op progress reporter."""
    
    def start_chromosome(self, chrom: str, total: int) -> None:
        pass
    
    def update(self, chrom: str, position: int) -> None:
        pass
    
    def finish_chromosome(self, chrom: str) -> None:
        pass
    
    def cleanup(self) -> None:
        pass


class ProgressCoordinator:
    """Coordinates progress reporting across execution modes.
    
    This class provides a unified interface for progress reporting,
    abstracting the differences between single-process and multi-process
    execution modes.
    """
    
    def __init__(self, config: ProgressConfig):
        """Initialize progress coordinator.
        
        Args:
            config: Progress configuration
        """
        self.config = config
        self._reporter: Optional[ProgressReporter] = None
        self._multiprocess_queue: Optional['Queue'] = None
        self._logger_lock: Optional['Lock'] = None
    
    def setup_single_process(self) -> ProgressReporter:
        """Set up progress for single process mode.
        
        Returns:
            Progress reporter for single process
        """
        if self.config.mode == ProgressMode.OBSERVER:
            self._reporter = ObserverProgressReporter(self.config.observers)
        elif self.config.mode == ProgressMode.DISABLED:
            self._reporter = NullProgressReporter()
        else:
            self._reporter = SingleLineProgressReporter()
        
        return self._reporter
    
    def setup_multiprocess(self, 
                          report_queue: 'Queue',
                          logger_lock: 'Lock') -> None:
        """Set up progress for multiprocess mode.
        
        Args:
            report_queue: Queue for progress reports
            logger_lock: Lock for thread-safe logging
        """
        self._multiprocess_queue = report_queue
        self._logger_lock = logger_lock
        
        # Configure global switches for multiprocess
        if self.config.enable_global:
            ProgressBar.global_switch = False
            ProgressHook.global_switch = True
        
        # Create multi-line reporter
        self._reporter = MultiLineProgressReporter(logger_lock)
    
    def handle_multiprocess_report(self, chrom: Optional[str], data: Any) -> bool:
        """Handle progress report from worker process.
        
        Args:
            chrom: Chromosome name (None for progress updates)
            data: Report data
            
        Returns:
            True if this was a progress update, False if result
        """
        if chrom is None and self._reporter:
            # Progress update
            actual_chrom, position = data
            self._reporter.update(actual_chrom, position)
            return True
        return False
    
    def get_reporter(self) -> Optional[ProgressReporter]:
        """Get current progress reporter."""
        return self._reporter
    
    def cleanup(self) -> None:
        """Clean up progress reporting."""
        if self._reporter:
            self._reporter.cleanup()
    
    @classmethod
    def create_for_single_process(cls, 
                                 use_observer: bool = False,
                                 observers: Optional[List['ProgressObserver']] = None,
                                 disabled: bool = False) -> 'ProgressCoordinator':
        """Create coordinator for single process execution.
        
        Args:
            use_observer: Whether to use observer pattern
            observers: List of progress observers
            disabled: Whether to disable progress
            
        Returns:
            Configured progress coordinator
        """
        if disabled:
            mode = ProgressMode.DISABLED
        elif use_observer:
            mode = ProgressMode.OBSERVER
        else:
            mode = ProgressMode.SINGLE_LINE
        
        config = ProgressConfig(mode=mode, observers=observers)
        return cls(config)
    
    @classmethod
    def create_for_multiprocess(cls) -> 'ProgressCoordinator':
        """Create coordinator for multiprocess execution.
        
        Returns:
            Configured progress coordinator
        """
        config = ProgressConfig(mode=ProgressMode.MULTI_LINE)
        return cls(config)