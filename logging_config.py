"""
Enhanced logging configuration for BioOps application.
Provides structured logging, performance tracking, and analysis monitoring.
"""

import os
import sys
import json
import time
import logging
import logging.config
from pathlib import Path
from typing import Dict, Any, Optional, Union
from datetime import datetime
from functools import wraps
from contextlib import contextmanager
import traceback

import yaml
import structlog
from pythonjsonlogger import jsonlogger
from prometheus_client import Counter, Histogram, Gauge
import sentry_sdk
from sentry_sdk.integrations.logging import LoggingIntegration

# Metrics
analysis_counter = Counter('bioops_analyses_total', 'Total number of analyses', ['protein_1', 'protein_2', 'status'])
analysis_duration = Histogram('bioops_analysis_duration_seconds', 'Analysis duration in seconds', ['step'])
active_analyses = Gauge('bioops_active_analyses', 'Number of active analyses')
error_counter = Counter('bioops_errors_total', 'Total number of errors', ['error_type', 'function'])

class PerformanceFilter(logging.Filter):
    """Add performance metrics to log records"""
    
    def filter(self, record):
        # Add memory usage
        try:
            import psutil
            process = psutil.Process(os.getpid())
            record.memory_mb = process.memory_info().rss / 1024 / 1024
            record.cpu_percent = process.cpu_percent(interval=0.1)
        except:
            record.memory_mb = 0
            record.cpu_percent = 0
        
        # Add timestamp
        record.timestamp = datetime.now().isoformat()
        
        return True

class AnalysisContextFilter(logging.Filter):
    """Add analysis context to log records"""
    
    def __init__(self):
        super().__init__()
        self.context = {}
    
    def set_context(self, **kwargs):
        """Set context variables"""
        self.context.update(kwargs)
    
    def clear_context(self):
        """Clear context variables"""
        self.context.clear()
    
    def filter(self, record):
        # Add context to record
        for key, value in self.context.items():
            setattr(record, key, value)
        return True

class StructuredFormatter(jsonlogger.JsonFormatter):
    """Enhanced JSON formatter with additional fields"""
    
    def add_fields(self, log_record, record, message_dict):
        super().add_fields(log_record, record, message_dict)
        
        # Add custom fields without overwriting existing ones
        if 'timestamp' not in log_record:
            log_record['timestamp'] = datetime.now().isoformat()
        if 'level' not in log_record:
            log_record['level'] = record.levelname
        if 'logger' not in log_record:
            log_record['logger'] = record.name
        
        # Add source location
        if 'source' not in log_record:
            log_record['source'] = {
                'file': record.pathname,
                'line': record.lineno,
                'function': record.funcName
            }
        
        # Add any extra fields from the record
        for key, value in record.__dict__.items():
            # Skip standard LogRecord attributes and 'message' to avoid conflicts
            if key not in ['name', 'msg', 'args', 'created', 'filename', 'funcName', 
                          'levelname', 'levelno', 'lineno', 'module', 'msecs', 
                          'message', 'pathname', 'process', 'processName', 
                          'relativeCreated', 'thread', 'threadName', 'exc_info',
                          'exc_text', 'stack_info', 'extra_fields'] and not key.startswith('_'):
                if key not in log_record:
                    log_record[key] = value

class BioOpsLogger:
    """Central logging configuration for BioOps"""
    
    def __init__(self, config_path: Optional[str] = None):
        self.config_path = config_path or "config/logging.yaml"
        self.context_filter = AnalysisContextFilter()
        self.performance_filter = PerformanceFilter()
        self._setup_structlog()
        self._setup_logging()
        self._setup_sentry()
        
    
    def _setup_logging(self):
        """Setup logging from configuration file"""
        try:
            # Load configuration
            config_file = Path(self.config_path)
            if config_file.exists():
                with open(config_file, 'r') as f:
                    config = yaml.safe_load(f)
                
                # Ensure log directory exists
                log_dir = Path("logs")
                log_dir.mkdir(exist_ok=True)
                
                # Apply configuration
                logging.config.dictConfig(config)
            else:
                # Fallback to basic configuration
                self._setup_basic_logging()
            
            # Add custom filters to all handlers
            for handler in logging.root.handlers:
                handler.addFilter(self.context_filter)
                handler.addFilter(self.performance_filter)
                
        except Exception as e:
            print(f"Error setting up logging: {e}")
            self._setup_basic_logging()
    
    def _setup_basic_logging(self):
        """Setup basic logging configuration as fallback"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.StreamHandler(sys.stdout),
                logging.FileHandler('logs/bioops.log', mode='a')
            ]
        )
    
    def _setup_sentry(self):
        """Setup Sentry error tracking if configured"""
        sentry_dsn = os.getenv('SENTRY_DSN')
        if sentry_dsn:
            sentry_logging = LoggingIntegration(
                level=logging.INFO,
                event_level=logging.ERROR
            )
            
            sentry_sdk.init(
                dsn=sentry_dsn,
                integrations=[sentry_logging],
                traces_sample_rate=float(os.getenv('SENTRY_TRACES_SAMPLE_RATE', 0.1)),
                environment=os.getenv('SENTRY_ENVIRONMENT', 'production')
            )
    
    def _setup_structlog(self):
        """Setup structured logging"""
        structlog.configure(
            processors=[
                structlog.stdlib.filter_by_level,
                structlog.stdlib.add_logger_name,
                structlog.stdlib.add_log_level,
                structlog.stdlib.PositionalArgumentsFormatter(),
                structlog.processors.TimeStamper(fmt="iso"),
                structlog.processors.StackInfoRenderer(),
                structlog.processors.format_exc_info,
                structlog.processors.UnicodeDecoder(),
                structlog.processors.CallsiteParameterAdder(
                    parameters=[
                        structlog.processors.CallsiteParameter.FILENAME,
                        structlog.processors.CallsiteParameter.FUNC_NAME,
                        structlog.processors.CallsiteParameter.LINENO,
                    ]
                ),
                structlog.processors.dict_tracebacks,
                structlog.dev.ConsoleRenderer() if os.getenv('APP_ENV') == 'development' else structlog.processors.JSONRenderer()
            ],
            context_class=dict,
            logger_factory=structlog.stdlib.LoggerFactory(),
            cache_logger_on_first_use=True,
        )
    
    def get_logger(self, name: str) -> logging.Logger:
        """Get a logger instance with custom configuration"""
        logger = logging.getLogger(name)
        return logger
    
    def set_context(self, **kwargs):
        """Set analysis context for all subsequent logs"""
        self.context_filter.set_context(**kwargs)
    
    def clear_context(self):
        """Clear analysis context"""
        self.context_filter.clear_context()

# Decorators for logging

def log_performance(logger_name: str = "bioops.performance"):
    """Decorator to log function performance"""
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            logger = logging.getLogger(logger_name)
            start_time = time.time()
            
            # Log function entry
            logger.debug(f"Entering {func.__name__}", extra={
                'function': func.__name__,
                'args': str(args)[:100],  # Truncate long args
                'kwargs': str(kwargs)[:100]
            })
            
            try:
                # Execute function
                result = func(*args, **kwargs)
                
                # Log success
                duration = time.time() - start_time
                logger.info(f"Completed {func.__name__}", extra={
                    'function': func.__name__,
                    'duration_seconds': duration,
                    'status': 'success'
                })
                
                # Update metrics
                analysis_duration.labels(step=func.__name__).observe(duration)
                
                return result
                
            except Exception as e:
                # Log error
                duration = time.time() - start_time
                logger.error(f"Error in {func.__name__}: {str(e)}", extra={
                    'function': func.__name__,
                    'duration_seconds': duration,
                    'status': 'error',
                    'error_type': type(e).__name__,
                    'error_message': str(e),
                    'traceback': traceback.format_exc()
                })
                
                # Update error metrics
                error_counter.labels(
                    error_type=type(e).__name__,
                    function=func.__name__
                ).inc()
                
                raise
        
        return wrapper
    return decorator

def log_analysis_step(step_name: str, logger_name: str = "bioops.analysis"):
    """Decorator to log analysis steps"""
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            logger = logging.getLogger(logger_name)
            
            # Log step start
            logger.info(f"Starting analysis step: {step_name}", extra={
                'step': step_name,
                'function': func.__name__
            })
            
            try:
                with analysis_duration.labels(step=step_name).time():
                    result = func(*args, **kwargs)
                
                # Log step completion
                logger.info(f"Completed analysis step: {step_name}", extra={
                    'step': step_name,
                    'function': func.__name__,
                    'status': 'success'
                })
                
                return result
                
            except Exception as e:
                # Log step failure
                logger.error(f"Failed analysis step: {step_name}", extra={
                    'step': step_name,
                    'function': func.__name__,
                    'status': 'error',
                    'error': str(e),
                    'traceback': traceback.format_exc()
                })
                raise
        
        return wrapper
    return decorator

@contextmanager
def log_context(**kwargs):
    """Context manager for temporary logging context"""
    logger_config = BioOpsLogger()
    logger_config.set_context(**kwargs)
    try:
        yield
    finally:
        logger_config.clear_context()

# Utility functions

def log_analysis_start(pdb_1_id: str, pdb_2_id: str, input_method: str):
    """Log the start of a new analysis"""
    logger = logging.getLogger("bioops.analysis")
    
    analysis_id = f"{pdb_1_id}_{pdb_2_id}_{int(time.time())}"
    
    logger.info("Starting new analysis", extra={
        'analysis_id': analysis_id,
        'pdb_1_id': pdb_1_id,
        'pdb_2_id': pdb_2_id,
        'input_method': input_method,
        'event': 'analysis_start'
    })
    
    # Update metrics
    active_analyses.inc()
    
    return analysis_id

def log_analysis_complete(analysis_id: str, pdb_1_id: str, pdb_2_id: str, success: bool):
    """Log the completion of an analysis"""
    logger = logging.getLogger("bioops.analysis")
    
    status = 'success' if success else 'failed'
    
    logger.info("Analysis complete", extra={
        'analysis_id': analysis_id,
        'pdb_1_id': pdb_1_id,
        'pdb_2_id': pdb_2_id,
        'status': status,
        'event': 'analysis_complete'
    })
    
    # Update metrics
    analysis_counter.labels(
        protein_1=pdb_1_id,
        protein_2=pdb_2_id,
        status=status
    ).inc()
    active_analyses.dec()

def log_tool_execution(tool_name: str, command: str, return_code: int, 
                      stdout: str, stderr: str, duration: float):
    """Log external tool execution"""
    logger = logging.getLogger("bioops.tools")
    
    log_data = {
        'tool': tool_name,
        'command': command[:200],  # Truncate long commands
        'return_code': return_code,
        'duration_seconds': duration,
        'success': return_code == 0 or (tool_name == "FATCAT" and return_code == 1)
    }
    
    if return_code != 0 and not (tool_name == "FATCAT" and return_code == 1):
        log_data['stderr'] = stderr[:500]  # Include stderr for failures
        logger.error(f"{tool_name} execution failed", extra=log_data)
    else:
        logger.info(f"{tool_name} execution completed", extra=log_data)

# Initialize global logger configuration
bioops_logger = BioOpsLogger()

# Export commonly used functions
__all__ = [
    'BioOpsLogger',
    'log_performance',
    'log_analysis_step',
    'log_context',
    'log_analysis_start',
    'log_analysis_complete',
    'log_tool_execution',
    'bioops_logger'
]
