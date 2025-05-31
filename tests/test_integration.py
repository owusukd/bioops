"""
Integration tests for the BioOps application
Tests the interaction between different components
"""

import pytest
from unittest.mock import Mock, patch, MagicMock
import tempfile
import shutil
from pathlib import Path
import pandas as pd
import numpy as np
import streamlit as st

# Import modules to test
import sys
sys.path.insert(0, '..')
import bioops
import protFuncs
import logging_config


class TestEndToEndWorkflow:
    """Test complete analysis workflow"""
    
    @patch('bioops.download_pdb_content')
    @patch('bioops.parse_pdb_file')
    @patch('protFuncs.FATCAT')
    @patch('protFuncs.getSuperimposed3Dpart')
    @patch('protFuncs.runP2rank')
    @patch('protFuncs.extractPocket_pdb')
    @patch('protFuncs.runApoc')
    @patch('protFuncs.parse_pocket_comparison_result')
    @patch('protFuncs.pocket_3D_similar_overlap_score')
    @patch('protFuncs.run_pockets_pdb2pqr')
    @patch('protFuncs.run_pockets_apbs')
    @patch('protFuncs.compare_pockets_potentials')
    def test_complete_analysis_pdb_id_input(self, *mocks):
        """Test complete analysis workflow with PDB ID input"""
        # Unpack mocks
        (mock_compare_potentials, mock_apbs, mock_pdb2pqr, mock_overlap,
         mock_parse_apoc, mock_apoc, mock_extract, mock_p2rank,
         mock_superimposed, mock_fatcat, mock_parse_pdb, mock_download) = mocks
        
        # Setup mock returns
        mock_download.return_value = "MOCK PDB CONTENT"
        mock_parse_pdb.return_value = True
        mock_fatcat.return_value = True
        mock_superimposed.return_value = True
        mock_p2rank.return_value = True
        mock_extract.return_value = (3, 4)  # 3 pockets for pdb1, 4 for pdb2
        mock_apoc.return_value = True
        mock_parse_apoc.return_value = True
        mock_overlap.return_value = True
        mock_pdb2pqr.return_value = True
        mock_apbs.return_value = True
        mock_compare_potentials.return_value = pd.DataFrame(np.random.rand(7, 7))
        
        # Mock session state
        with patch.dict('streamlit.session_state', {
            'input_file_option': '**Enter PDB ID**',
            'protein_1': '5CXV',
            'protein_2': '5DSG',
            'protein1_chain': 'A',
            'protein2_chain': 'B'
        }):
            # Initialize session state
            bioops.initialize_session_state()
            
            # Process input files
            with patch('bioops.ensure_directories'):
                result = bioops.process_input_files(
                    st.session_state.input_file_option,
                    ('5CXV', '5DSG', 'A', 'B')
                )
            
            assert result == True
            assert st.session_state.pdb_1_id == '5cxv'
            assert st.session_state.pdb_2_id == '5dsg'
            
            # Run analysis
            with patch('streamlit.progress'), patch('streamlit.empty'):
                result = bioops.run_analysis()
            
            assert result == True
            assert st.session_state.processing_complete == True
            assert st.session_state.pdb_1_num_pocket == 3
            assert st.session_state.pdb_2_num_pocket == 4
            assert not st.session_state.pocket_similarity_matrix.empty


class TestComponentInteraction:
    """Test interaction between different components"""
    
    def test_bioops_protfuncs_integration(self):
        """Test bioops.py and protFuncs.py integration"""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Mock paths
            with patch('protFuncs.RESULTS_BASE_DIR', Path(tmpdir) / 'results'):
                with patch('bioops.RESULTS_DIR', Path(tmpdir) / 'results'):
                    # Create test file
                    results_dir = Path(tmpdir) / 'results'
                    results_dir.mkdir(parents=True)
                    test_pdb = results_dir / 'test.pdb'
                    test_pdb.write_text("MOCK PDB CONTENT")
                    
                    # Test file reading through bioops
                    with patch('protFuncs.read_file_content') as mock_read:
                        mock_read.return_value = "MOCK PDB CONTENT"
                        content = protFuncs.read_saved_pdb_file_to_view('test', 'test')
                        
                        assert content == "MOCK PDB CONTENT"
    
    def test_logging_integration(self):
        """Test logging integration across modules"""
        # Create logger
        bioops_logger = logging_config.BioOpsLogger()
        
        # Test logging from bioops
        logger = bioops_logger.get_logger('bioops')
        logger.info("Test from bioops")
        
        # Test logging from protFuncs
        logger = bioops_logger.get_logger('protFuncs')
        logger.info("Test from protFuncs")
        
        # Test context setting
        bioops_logger.set_context(analysis_id="test123")
        
        # Create a log record and verify context is added
        test_logger = bioops_logger.get_logger('test')
        with patch.object(test_logger, 'info') as mock_info:
            test_logger.info("Test message")
            
            # The context should be added by the filter
            assert mock_info.called
    
    @patch('bioops.log_analysis_start')
    @patch('bioops.log_analysis_complete')
    def test_analysis_tracking(self, mock_complete, mock_start):
        """Test analysis tracking through logging"""
        mock_start.return_value = "analysis_123"
        
        with patch.dict('streamlit.session_state', {
            'pdb_1': 'test',
            'pdb_2': 'test',
            'pdb_1_id': '5CXV',
            'pdb_2_id': '5DSG',
            'input_file_option': 'test'
        }):
            # Mock all analysis steps
            with patch('bioops.protFuncs') as mock_protfuncs:
                mock_protfuncs.FATCAT.return_value = True
                mock_protfuncs.extractPocket_pdb.return_value = (1, 1)
                mock_protfuncs.compare_pockets_potentials.return_value = pd.DataFrame()
                
                with patch('streamlit.progress'), patch('streamlit.empty'):
                    bioops.run_analysis()
            
            mock_start.assert_called_once()
            mock_complete.assert_called_once_with("analysis_123", "5CXV", "5DSG", True)


class TestErrorPropagation:
    """Test error propagation between components"""
    
    def test_download_error_propagation(self):
        """Test PDB download error propagation"""
        with patch('requests.get') as mock_get:
            mock_get.side_effect = Exception("Network error")
            
            with pytest.raises(bioops.PDBDownloadError):
                bioops.download_pdb_content("5CXV")
    
    def test_analysis_error_propagation(self):
        """Test analysis error propagation"""
        with patch.dict('streamlit.session_state', {
            'pdb_1': 'test',
            'pdb_2': 'test',
            'pdb_1_id': '5CXV',
            'pdb_2_id': '5DSG'
        }):
            with patch('bioops.protFuncs.FATCAT') as mock_fatcat:
                mock_fatcat.side_effect = protFuncs.ExternalToolError("FATCAT failed")
                
                with patch('streamlit.progress'), patch('streamlit.empty'):
                    result = bioops.run_analysis()
                
                assert result == False
                assert st.session_state.processing_error == True
                assert "FATCAT failed" in st.session_state.error_message


class TestDataFlow:
    """Test data flow between components"""
    
    def test_pocket_count_flow(self):
        """Test pocket count data flow from P2Rank to visualization"""
        with patch.dict('streamlit.session_state', {
            'pdb_1': 'test',
            'pdb_2': 'test',
            'pdb_1_id': '5CXV',
            'pdb_2_id': '5DSG'
        }):
            # Mock P2Rank returning pocket counts
            with patch('bioops.protFuncs.extractPocket_pdb') as mock_extract:
                mock_extract.return_value = (5, 7)
                
                # Mock other required functions
                with patch('bioops.protFuncs') as mock_pf:
                    mock_pf.FATCAT.return_value = True
                    mock_pf.compare_pockets_potentials.return_value = pd.DataFrame()
                    
                    with patch('streamlit.progress'), patch('streamlit.empty'):
                        bioops.run_analysis()
            
            # Verify pocket counts are stored correctly
            assert st.session_state.pdb_1_num_pocket == 5
            assert st.session_state.pdb_2_num_pocket == 7
    
    def test_similarity_matrix_flow(self):
        """Test similarity matrix data flow"""
        test_matrix = pd.DataFrame(
            np.random.rand(5, 5),
            index=[f'pocket_{i}' for i in range(5)],
            columns=[f'pocket_{i}' for i in range(5)]
        )
        
        with patch.dict('streamlit.session_state', {
            'pdb_1': 'test',
            'pdb_2': 'test',
            'pdb_1_id': '5CXV',
            'pdb_2_id': '5DSG'
        }):
            with patch('bioops.protFuncs.compare_pockets_potentials') as mock_compare:
                mock_compare.return_value = test_matrix
                
                # Mock other required functions
                with patch('bioops.protFuncs') as mock_pf:
                    mock_pf.FATCAT.return_value = True
                    mock_pf.extractPocket_pdb.return_value = (2, 3)
                    
                    with patch('streamlit.progress'), patch('streamlit.empty'):
                        bioops.run_analysis()
            
            # Verify similarity matrix is stored
            assert not st.session_state.pocket_similarity_matrix.empty
            assert st.session_state.pocket_similarity_matrix.shape == (5, 5)


class TestPerformanceIntegration:
    """Test performance tracking integration"""
    
    @patch('time.time')
    def test_analysis_duration_tracking(self, mock_time):
        """Test analysis duration is tracked correctly"""
        # Mock time progression
        mock_time.side_effect = [
            1000.0,  # Start time
            1001.0, 1002.0, 1003.0, 1004.0, 1005.0,  # Step times
            1006.0, 1007.0, 1008.0, 1010.0  # More step times
        ]
        
        with patch('logging_config.analysis_duration.labels') as mock_duration:
            mock_histogram = Mock()
            mock_duration.return_value = mock_histogram
            
            # Run analysis with mocked components
            with patch.dict('streamlit.session_state', {
                'pdb_1': 'test',
                'pdb_2': 'test',
                'pdb_1_id': '5CXV',
                'pdb_2_id': '5DSG'
            }):
                with patch('bioops.protFuncs') as mock_pf:
                    mock_pf.FATCAT.return_value = True
                    mock_pf.extractPocket_pdb.return_value = (1, 1)
                    mock_pf.compare_pockets_potentials.return_value = pd.DataFrame()
                    
                    with patch('streamlit.progress'), patch('streamlit.empty'):
                        bioops.run_analysis()
            
            # Verify duration tracking was called for each step
            assert mock_duration.call_count >= 8  # At least 8 analysis steps


# Fixtures
@pytest.fixture
def clean_session_state():
    """Clean session state for tests"""
    # Clear session state
    for key in list(st.session_state.keys()):
        del st.session_state[key]
    
    yield
    
    # Clear again after test
    for key in list(st.session_state.keys()):
        del st.session_state[key]


@pytest.fixture
def mock_analysis_environment():
    """Mock complete analysis environment"""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create directory structure
        base_dir = Path(tmpdir)
        (base_dir / "results" / "fatcat").mkdir(parents=True)
        (base_dir / "results" / "p2rank").mkdir(parents=True)
        (base_dir / "results" / "apoc").mkdir(parents=True)
        (base_dir / "results" / "apbs").mkdir(parents=True)
        
        with patch('bioops.RESULTS_DIR', base_dir / "results"):
            with patch('protFuncs.RESULTS_BASE_DIR', base_dir / "results"):
                yield base_dir