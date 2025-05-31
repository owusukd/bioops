"""
Tests for logging_config.py - Enhanced logging configuration
"""

import pytest
from unittest.mock import Mock, patch, MagicMock
import logging
import time
from pathlib import Path
import json

# Import the module to test
import sys
sys.path.insert(0, '..')
import logging_config


class TestLoggingConfiguration:
    """Test logging configuration setup"""
    
    @patch('pathlib.Path.exists')
    @patch('builtins.open')
    @patch('yaml.safe_load')
    def test_bioops_logger_init_with_config(self, mock_yaml, mock_open, mock_exists):
        """Test BioOpsLogger initialization with config file"""
        mock_exists.return_value = True
        mock_yaml.return_value = {
            'version': 1,
            'handlers': {},
            'loggers': {}
        }
        
        with patch('logging.config.dictConfig') as mock_dict_config:
            logger = logging_config.BioOpsLogger()
            
            mock_dict_config.assert_called_once()
            assert logger.context_filter is not None
            assert logger.performance_filter is not None
    
    @patch('pathlib.Path.exists')
    def test_bioops_logger_init_without_config(self, mock_exists):
        """Test BioOpsLogger initialization without config file"""
        mock_exists.return_value = False
        
        with patch('logging.basicConfig') as mock_basic_config:
            logger = logging_config.BioOpsLogger()
            
            mock_basic_config.assert_called_once()
    
    def test_get_logger(self):
        """Test logger retrieval"""
        bioops_logger = logging_config.BioOpsLogger()
        logger = bioops_logger.get_logger("test_logger")
        
        assert isinstance(logger, logging.Logger)
        assert logger.name == "test_logger"


class TestFilters:
    """Test custom logging filters"""
    
    @patch('psutil.Process')
    def test_performance_filter(self, mock_process):
        """Test PerformanceFilter adds metrics"""
        # Mock process metrics
        mock_proc_instance = Mock()
        mock_proc_instance.memory_info.return_value.rss = 1024 * 1024 * 100  # 100MB
        mock_proc_instance.cpu_percent.return_value = 50.0
        mock_process.return_value = mock_proc_instance
        
        filter_obj = logging_config.PerformanceFilter()
        record = logging.LogRecord(
            name="test",
            level=logging.INFO,
            pathname="test.py",
            lineno=1,
            msg="Test message",
            args=(),
            exc_info=None
        )
        
        result = filter_obj.filter(record)
        
        assert result == True
        assert hasattr(record, 'memory_mb')
        assert hasattr(record, 'cpu_percent')
        assert hasattr(record, 'timestamp')
        assert record.memory_mb == 100.0
        assert record.cpu_percent == 50.0
    
    def test_analysis_context_filter(self):
        """Test AnalysisContextFilter"""
        filter_obj = logging_config.AnalysisContextFilter()
        
        # Set context
        filter_obj.set_context(analysis_id="test123", pdb_1="5CXV")
        
        record = logging.LogRecord(
            name="test",
            level=logging.INFO,
            pathname="test.py",
            lineno=1,
            msg="Test message",
            args=(),
            exc_info=None
        )
        
        result = filter_obj.filter(record)
        
        assert result == True
        assert hasattr(record, 'analysis_id')
        assert hasattr(record, 'pdb_1')
        assert record.analysis_id == "test123"
        assert record.pdb_1 == "5CXV"
        
        # Clear context
        filter_obj.clear_context()
        
        new_record = logging.LogRecord(
            name="test",
            level=logging.INFO,
            pathname="test.py",
            lineno=1,
            msg="Test message",
            args=(),
            exc_info=None
        )
        
        filter_obj.filter(new_record)
        assert not hasattr(new_record, 'analysis_id')


class TestDecorators:
    """Test logging decorators"""
    
    @patch('logging.getLogger')
    def test_log_performance_success(self, mock_get_logger):
        """Test performance logging decorator on success"""
        mock_logger = Mock()
        mock_get_logger.return_value = mock_logger
        
        @logging_config.log_performance()
        def test_function(x, y):
            return x + y
        
        result = test_function(1, 2)
        
        assert result == 3
        assert mock_logger.debug.call_count == 1  # Entry log
        assert mock_logger.info.call_count == 1   # Success log
        
        # Check logged data
        info_call = mock_logger.info.call_args
        assert "Completed test_function" in info_call[0][0]
        assert info_call[1]['extra']['status'] == 'success'
        assert 'duration_seconds' in info_call[1]['extra']
    
    @patch('logging.getLogger')
    def test_log_performance_failure(self, mock_get_logger):
        """Test performance logging decorator on failure"""
        mock_logger = Mock()
        mock_get_logger.return_value = mock_logger
        
        @logging_config.log_performance()
        def test_function():
            raise ValueError("Test error")
        
        with pytest.raises(ValueError):
            test_function()
        
        assert mock_logger.debug.call_count == 1  # Entry log
        assert mock_logger.error.call_count == 1  # Error log
        
        # Check error logging
        error_call = mock_logger.error.call_args
        assert "Error in test_function" in error_call[0][0]
        assert error_call[1]['extra']['status'] == 'error'
        assert error_call[1]['extra']['error_type'] == 'ValueError'
    
    @patch('logging.getLogger')
    def test_log_analysis_step_success(self, mock_get_logger):
        """Test analysis step logging decorator"""
        mock_logger = Mock()
        mock_get_logger.return_value = mock_logger
        
        @logging_config.log_analysis_step("Test Step")
        def test_step():
            return "success"
        
        with patch('logging_config.analysis_duration.labels') as mock_labels:
            mock_timer = Mock()
            mock_labels.return_value.time.return_value.__enter__ = Mock()
            mock_labels.return_value.time.return_value.__exit__ = Mock()
            
            result = test_step()
        
        assert result == "success"
        assert mock_logger.info.call_count == 2  # Start and complete
    
    def test_log_context_context_manager(self):
        """Test log_context context manager"""
        logger_config = Mock()
        
        with patch('logging_config.BioOpsLogger', return_value=logger_config):
            with logging_config.log_context(test_id="123"):
                logger_config.set_context.assert_called_once_with(test_id="123")
            
            logger_config.clear_context.assert_called_once()


class TestUtilityFunctions:
    """Test utility functions"""
    
    @patch('logging.getLogger')
    @patch('time.time')
    def test_log_analysis_start(self, mock_time, mock_get_logger):
        """Test analysis start logging"""
        mock_time.return_value = 1234567890
        mock_logger = Mock()
        mock_get_logger.return_value = mock_logger
        
        with patch('logging_config.active_analyses.inc') as mock_inc:
            analysis_id = logging_config.log_analysis_start("5CXV", "5DSG", "upload")
        
        assert analysis_id == "5CXV_5DSG_1234567890"
        mock_logger.info.assert_called_once()
        mock_inc.assert_called_once()
        
        # Check logged data
        info_call = mock_logger.info.call_args
        assert info_call[1]['extra']['analysis_id'] == analysis_id
        assert info_call[1]['extra']['pdb_1_id'] == "5CXV"
        assert info_call[1]['extra']['pdb_2_id'] == "5DSG"
    
    @patch('logging.getLogger')
    def test_log_analysis_complete_success(self, mock_get_logger):
        """Test analysis completion logging - success"""
        mock_logger = Mock()
        mock_get_logger.return_value = mock_logger
        
        with patch('logging_config.analysis_counter.labels') as mock_counter:
            with patch('logging_config.active_analyses.dec') as mock_dec:
                logging_config.log_analysis_complete("test_id", "5CXV", "5DSG", True)
        
        mock_logger.info.assert_called_once()
        mock_counter.assert_called_once_with(protein_1="5CXV", protein_2="5DSG", status="success")
        mock_dec.assert_called_once()
    
    @patch('logging.getLogger')
    def test_log_tool_execution_success(self, mock_get_logger):
        """Test tool execution logging - success"""
        mock_logger = Mock()
        mock_get_logger.return_value = mock_logger
        
        logging_config.log_tool_execution(
            "FATCAT",
            "FATCAT -p1 test1.pdb -p2 test2.pdb",
            0,
            "Success output",
            "",
            5.23
        )
        
        mock_logger.info.assert_called_once()
        info_call = mock_logger.info.call_args
        assert "FATCAT execution completed" in info_call[0][0]
        assert info_call[1]['extra']['tool'] == "FATCAT"
        assert info_call[1]['extra']['success'] == True
        assert info_call[1]['extra']['duration_seconds'] == 5.23
    
    @patch('logging.getLogger')
    def test_log_tool_execution_failure(self, mock_get_logger):
        """Test tool execution logging - failure"""
        mock_logger = Mock()
        mock_get_logger.return_value = mock_logger
        
        logging_config.log_tool_execution(
            "P2RANK",
            "prank predict -f test.pdb",
            1,
            "",
            "Error: File not found",
            2.5
        )
        
        mock_logger.error.assert_called_once()
        error_call = mock_logger.error.call_args
        assert "P2RANK execution failed" in error_call[0][0]
        assert error_call[1]['extra']['stderr'] == "Error: File not found"


class TestFormatters:
    """Test custom formatters"""
    
    def test_structured_formatter(self):
        """Test StructuredFormatter"""
        formatter = logging_config.StructuredFormatter()
        
        record = logging.LogRecord(
            name="test.logger",
            level=logging.INFO,
            pathname="/app/test.py",
            lineno=42,
            msg="Test message",
            args=(),
            exc_info=None
        )
        
        # Add extra fields
        record.extra_fields = {"custom_field": "custom_value"}
        
        # Format the record
        formatted = formatter.format(record)
        
        # Parse JSON output
        parsed = json.loads(formatted)
        
        assert parsed['level'] == 'INFO'
        assert parsed['logger'] == 'test.logger'
        assert parsed['message'] == 'Test message'
        assert parsed['source']['file'] == '/app/test.py'
        assert parsed['source']['line'] == 42
        assert parsed['source']['function'] == '<module>'
        assert 'timestamp' in parsed


class TestMetrics:
    """Test Prometheus metrics"""
    
    def test_metrics_defined(self):
        """Test that all metrics are properly defined"""
        assert logging_config.analysis_counter is not None
        assert logging_config.analysis_duration is not None
        assert logging_config.active_analyses is not None
        assert logging_config.error_counter is not None
        
        # Check metric types
        assert hasattr(logging_config.analysis_counter, 'inc')
        assert hasattr(logging_config.analysis_duration, 'observe')
        assert hasattr(logging_config.active_analyses, 'inc')
        assert hasattr(logging_config.active_analyses, 'dec')
        assert hasattr(logging_config.error_counter, 'inc')


class TestSentryIntegration:
    """Test Sentry integration"""
    
    @patch.dict('os.environ', {'SENTRY_DSN': 'https://test@sentry.io/123'})
    @patch('sentry_sdk.init')
    def test_sentry_setup_with_dsn(self, mock_sentry_init):
        """Test Sentry setup when DSN is provided"""
        logging_config.BioOpsLogger()
        
        mock_sentry_init.assert_called_once()
        call_args = mock_sentry_init.call_args[1]
        assert call_args['dsn'] == 'https://test@sentry.io/123'
        assert 'integrations' in call_args
    
    @patch.dict('os.environ', {}, clear=True)
    @patch('sentry_sdk.init')
    def test_sentry_setup_without_dsn(self, mock_sentry_init):
        """Test Sentry setup when DSN is not provided"""
        logging_config.BioOpsLogger()
        
        mock_sentry_init.assert_not_called()


# Fixtures
@pytest.fixture(autouse=True)
def reset_logging():
    """Reset logging configuration between tests"""
    # Remove all handlers
    logger = logging.getLogger()
    handlers = logger.handlers[:]
    for handler in handlers:
        logger.removeHandler(handler)
    
    yield
    
    # Clean up after tests
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
        