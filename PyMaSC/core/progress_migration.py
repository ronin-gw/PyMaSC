"""Migration utilities for transitioning to observer-based progress.

This module provides utilities and monkey-patches to help migrate
existing code to use the new observer-based progress system while
maintaining backward compatibility.
"""
import logging
from typing import Optional, Type, Callable, Any
from functools import wraps

from PyMaSC.utils.progress import ProgressBar, ProgressHook, MultiLineProgressManager
from .progress_adapter import (
    ProgressBarAdapter, ProgressHookAdapter, ReadCountProgressBarAdapter,
    get_progress_manager, ProgressManager
)
from .observer import ProgressObserver, ProgressSubject

logger = logging.getLogger(__name__)


# Global flag to control migration
_migration_enabled = False
_original_classes = {}


def enable_progress_migration(enabled: bool = True) -> None:
    """Enable or disable global progress migration.
    
    When enabled, all progress bar instantiations will automatically
    use the observer-aware adapters instead of the original classes.
    
    Args:
        enabled: Whether to enable migration
    """
    global _migration_enabled
    _migration_enabled = enabled
    
    if enabled:
        _apply_migration_patches()
        logger.info("Progress migration enabled - using observer adapters")
    else:
        _restore_original_classes()
        logger.info("Progress migration disabled - using original classes")


def _apply_migration_patches():
    """Apply monkey patches to replace progress classes with adapters."""
    import PyMaSC.utils.progress as progress_module
    
    # Store originals if not already stored
    if not _original_classes:
        _original_classes['ProgressBar'] = progress_module.ProgressBar
        _original_classes['ProgressHook'] = progress_module.ProgressHook
        _original_classes['ReadCountProgressBar'] = progress_module.ReadCountProgressBar
    
    # Replace with adapter factories
    progress_module.ProgressBar = _create_progress_bar_factory()
    progress_module.ProgressHook = _create_progress_hook_factory()
    progress_module.ReadCountProgressBar = _create_read_count_factory()


def _restore_original_classes():
    """Restore original progress classes."""
    if _original_classes:
        import PyMaSC.utils.progress as progress_module
        progress_module.ProgressBar = _original_classes['ProgressBar']
        progress_module.ProgressHook = _original_classes['ProgressHook']
        progress_module.ReadCountProgressBar = _original_classes['ReadCountProgressBar']


def _create_progress_bar_factory():
    """Create a factory that returns ProgressBar or adapter based on config."""
    def factory(*args, **kwargs):
        if _migration_enabled:
            adapter = ProgressBarAdapter(*args, **kwargs)
            # Auto-attach to global manager if available
            if hasattr(factory, '_auto_attach_observers'):
                manager = get_progress_manager()
                for observer in factory._auto_attach_observers:
                    adapter._subject.attach(observer)
            return adapter
        else:
            return _original_classes['ProgressBar'](*args, **kwargs)
    
    factory._auto_attach_observers = []
    return factory


def _create_progress_hook_factory():
    """Create a factory that returns ProgressHook or adapter based on config."""
    def factory(*args, **kwargs):
        if _migration_enabled:
            return ProgressHookAdapter(*args, **kwargs)
        else:
            return _original_classes['ProgressHook'](*args, **kwargs)
    return factory


def _create_read_count_factory():
    """Create a factory that returns ReadCountProgressBar or adapter."""
    def factory(*args, **kwargs):
        if _migration_enabled and 'ReadCountProgressBar' in _original_classes:
            # Use original class directly in adapter to avoid recursion
            original_class = _original_classes['ReadCountProgressBar']
            
            # Create adapter manually to avoid recursion
            adapter = object.__new__(ReadCountProgressBarAdapter)
            adapter._bar = original_class(*args, **kwargs)
            adapter._subject = ProgressSubject()
            adapter._current_chrom = None
            adapter._current_genome_total = None
            adapter._current_chrom_total = None
            adapter._genome_progress = 0
            
            # Auto-attach observers
            if hasattr(factory, '_auto_attach_observers'):
                for observer in factory._auto_attach_observers:
                    adapter.attach_observer(observer)
            return adapter
        else:
            return _original_classes.get('ReadCountProgressBar', 
                                        ReadCountProgressBar)(*args, **kwargs)
    
    factory._auto_attach_observers = []
    return factory


def auto_attach_observer(observer: ProgressObserver, 
                        progress_class: Optional[Type] = None) -> None:
    """Automatically attach an observer to all new progress instances.
    
    Args:
        observer: Observer to auto-attach
        progress_class: Specific progress class to target (None for all)
    """
    import PyMaSC.utils.progress as progress_module
    
    if progress_class is None or progress_class == ProgressBar:
        if hasattr(progress_module.ProgressBar, '_auto_attach_observers'):
            progress_module.ProgressBar._auto_attach_observers.append(observer)
    
    if progress_class is None or progress_class.__name__ == 'ReadCountProgressBar':
        if hasattr(progress_module.ReadCountProgressBar, '_auto_attach_observers'):
            progress_module.ReadCountProgressBar._auto_attach_observers.append(observer)


def with_progress_observer(observer_class: Type[ProgressObserver], 
                          **observer_kwargs):
    """Decorator to automatically attach observers to progress bars in a function.
    
    Args:
        observer_class: Observer class to instantiate
        **observer_kwargs: Arguments for observer initialization
        
    Example:
        >>> from PyMaSC.core.observer import FileProgressObserver
        >>> @with_progress_observer(FileProgressObserver, log_path="progress.log")
        ... def process_data():
        ...     progress = ProgressBar()
        ...     # Observer automatically attached
    """
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Create observer
            observer = observer_class(**observer_kwargs)
            
            # Temporarily enable migration and auto-attach
            original_state = _migration_enabled
            enable_progress_migration(True)
            
            try:
                auto_attach_observer(observer)
                result = func(*args, **kwargs)
                return result
            finally:
                # Restore original state
                if not original_state:
                    enable_progress_migration(False)
        
        return wrapper
    return decorator


class ProgressMigrationContext:
    """Context manager for temporary progress migration.
    
    Example:
        >>> with ProgressMigrationContext() as ctx:
        ...     ctx.attach_observer(FileProgressObserver("log.txt"))
        ...     # All progress bars created here will have observer
        ...     progress = ProgressBar()
    """
    
    def __init__(self, enable: bool = True):
        """Initialize migration context.
        
        Args:
            enable: Whether to enable migration in this context
        """
        self.enable = enable
        self._original_state = None
        self._observers = []
    
    def attach_observer(self, observer: ProgressObserver) -> None:
        """Attach an observer for this context.
        
        Args:
            observer: Observer to attach
        """
        self._observers.append(observer)
    
    def __enter__(self):
        """Enter context and enable migration."""
        self._original_state = _migration_enabled
        self._original_auto_attach = []
        
        if self.enable:
            enable_progress_migration(True)
            
            # Store and update auto-attach observers
            import PyMaSC.utils.progress as progress_module
            
            # For each observer, add to auto-attach list
            for observer in self._observers:
                if hasattr(progress_module.ProgressBar, '_auto_attach_observers'):
                    progress_module.ProgressBar._auto_attach_observers.append(observer)
                if hasattr(progress_module.ReadCountProgressBar, '_auto_attach_observers'):
                    progress_module.ReadCountProgressBar._auto_attach_observers.append(observer)
        
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Exit context and restore state."""
        # Clean up auto-attach observers
        if self.enable and self._observers:
            import PyMaSC.utils.progress as progress_module
            
            for observer in self._observers:
                if hasattr(progress_module.ProgressBar, '_auto_attach_observers'):
                    try:
                        progress_module.ProgressBar._auto_attach_observers.remove(observer)
                    except ValueError:
                        pass
                if hasattr(progress_module.ReadCountProgressBar, '_auto_attach_observers'):
                    try:
                        progress_module.ReadCountProgressBar._auto_attach_observers.remove(observer)
                    except ValueError:
                        pass
        
        # Restore original state
        if self._original_state is not None:
            enable_progress_migration(self._original_state)


def create_handler_with_observers(handler_class: Type,
                                 observers: list,
                                 *args, **kwargs) -> Any:
    """Create a handler instance with progress observers attached.
    
    Args:
        handler_class: Handler class to instantiate
        observers: List of observers to attach
        *args, **kwargs: Arguments for handler initialization
        
    Returns:
        Handler instance with observers attached
        
    Example:
        >>> from PyMaSC.handler.unified import UnifiedCalcHandler
        >>> from PyMaSC.core.observer import TerminalProgressObserver
        >>> handler = create_handler_with_observers(
        ...     UnifiedCalcHandler,
        ...     [TerminalProgressObserver()],
        ...     path="input.bam",
        ...     config=calc_config
        ... )
    """
    handler = handler_class(*args, **kwargs)
    
    # Check if handler supports observer attachment
    if hasattr(handler, 'attach_progress_observer'):
        for observer in observers:
            handler.attach_progress_observer(observer)
    else:
        logger.warning(f"{handler_class.__name__} does not support progress observers")
    
    return handler


# Convenience functions for common migration scenarios
def enable_file_logging(log_path: str, format_type: str = 'json') -> ProgressObserver:
    """Enable file logging for all progress events.
    
    Args:
        log_path: Path to log file
        format_type: Log format ('json' or 'text')
        
    Returns:
        The created file observer
    """
    from .observer import FileProgressObserver
    
    observer = FileProgressObserver(log_path, format_type)
    enable_progress_migration(True)
    auto_attach_observer(observer)
    
    return observer


def enable_aggregate_tracking(expected_sources: Optional[list] = None) -> ProgressObserver:
    """Enable aggregate progress tracking.
    
    Args:
        expected_sources: Optional list of expected sources
        
    Returns:
        The created aggregate observer
    """
    from .observer import AggregateProgressObserver
    
    observer = AggregateProgressObserver(expected_sources)
    enable_progress_migration(True)
    auto_attach_observer(observer)
    
    return observer