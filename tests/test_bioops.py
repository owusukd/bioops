"""
Tests for bioops.py - Main Streamlit application
"""

import pytest
import streamlit as st
from unittest.mock import Mock, patch, MagicMock
import tempfile
import shutil
from pathlib import Path
import pandas as pd

# Import the module to test
import sys
sys.path.insert(0, '..')
import bioops


class TestBioOpsConfiguration:
    """Test configuration and setup functions"""
    
    def test_get_app_config(self):
        """Test application configuration retrieval"""
        config = bioops.get_app_config()
        
        assert "logo_main" in config
        assert "logo_icon" in config
        assert "temp_dir" in config
        assert "results_dir" in config
        assert isinstance(config["temp_dir"], Path)
        assert isinstance(config["results_dir"], Path)
    
    def test_ensure_directories(self):
        """Test directory creation"""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Mock the config to use temp directory
            with patch('bioops.get_app_config') as mock_config:
                mock_config.return_value = {
                    "temp_dir": Path(tmpdir) / "temp",
                    "results_dir": Path(tmpdir) / "results"
                }
                
                bioops.ensure_directories()
                
                assert (Path(tmpdir) / "temp").exists()
                assert (Path(tmpdir) / "results").exists()


class TestInputValidation:
    """Test input validation functions"""
    
    def test_validate_pdb_id_valid(self):
        """Test valid PDB ID validation"""
        assert bioops.validate_pdb_id("5CXV") == True
        assert bioops.validate_pdb_id("1ABC") == True
        assert bioops.validate_pdb_id("7xyz") == True
    
    def test_validate_pdb_id_invalid(self):
        """Test invalid PDB ID validation"""
        assert bioops.validate_pdb_id("") == False
        assert bioops.validate_pdb_id("ABC") == False
        assert bioops.validate_pdb_id("12345") == False
        assert bioops.validate_pdb_id("AB CD") == False
        assert bioops.validate_pdb_id("ABCDE") == False
    
    def test_validate_chain_id_valid(self):
        """Test valid chain ID validation"""
        assert bioops.validate_chain_id("A") == True
        assert bioops.validate_chain_id("B") == True
        assert bioops.validate_chain_id("Z") == True
        assert bioops.validate_chain_id("") == True  # Empty is valid (defaults to A)
    
    def test_validate_chain_id_invalid(self):
        """Test invalid chain ID validation"""
        assert bioops.validate_chain_id("AB") == False
        assert bioops.validate_chain_id("1") == False
        assert bioops.validate_chain_id(" ") == False
    
    def test_validate_uploaded_file_valid(self):
        """Test valid uploaded file validation"""
        mock_file = Mock()
        mock_file.name = "protein.pdb"
        mock_file.size = 1024 * 1024  # 1MB
        
        assert bioops.validate_uploaded_file(mock_file) == True
    
    def test_validate_uploaded_file_invalid(self):
        """Test invalid uploaded file validation"""
        # None file
        assert bioops.validate_uploaded_file(None) == False
        
        # Wrong extension
        mock_file = Mock()
        mock_file.name = "protein.txt"
        mock_file.size = 1024
        assert bioops.validate_uploaded_file(mock_file) == False
        
        # Too large
        mock_file = Mock()
        mock_file.name = "protein.pdb"
        mock_file.size = 100 * 1024 * 1024  # 100MB
        assert bioops.validate_uploaded_file(mock_file) == False


class TestErrorHandling:
    """Test error handling classes and functions"""
    
    def test_bioops_error_hierarchy(self):
        """Test custom exception hierarchy"""
        assert issubclass(bioops.FileProcessingError, bioops.BioOpsError)
        assert issubclass(bioops.PDBDownloadError, bioops.BioOpsError)
        assert issubclass(bioops.AnalysisError, bioops.BioOpsError)
    
    @patch('streamlit.error')
    def test_show_error(self, mock_error):
        """Test error display function"""
        bioops.show_error("Test Error", "This is a test error")
        mock_error.assert_called_once()
        
    @patch('streamlit.warning')
    def test_show_warning(self, mock_warning):
        """Test warning display function"""
        bioops.show_error("Test Warning", "This is a test warning", "warning")
        mock_warning.assert_called_once()


class TestPDBFileHandling:
    """Test PDB file handling functions"""
    
    def test_chain_select(self):
        """Test ChainSelect class"""
        selector = bioops.ChainSelect("A")
        
        # Mock chain
        mock_chain = Mock()
        mock_chain.get_id.return_value = "A"
        assert selector.accept_chain(mock_chain) == True
        
        mock_chain.get_id.return_value = "B"
        assert selector.accept_chain(mock_chain) == False
    
    @patch('requests.get')
    def test_download_pdb_content_success(self, mock_get):
        """Test successful PDB download"""
        mock_response = Mock()
        mock_response.text = "MOCK PDB CONTENT"
        mock_response.raise_for_status = Mock()
        mock_get.return_value = mock_response
        
        content = bioops.download_pdb_content("5CXV")
        
        assert content == "MOCK PDB CONTENT"
        mock_get.assert_called_once_with(
            "https://files.rcsb.org/download/5cxv.pdb",
            timeout=30
        )
    
    @patch('requests.get')
    def test_download_pdb_content_failure(self, mock_get):
        """Test failed PDB download"""
        mock_get.side_effect = Exception("Network error")
        
        with pytest.raises(bioops.PDBDownloadError):
            bioops.download_pdb_content("5CXV")
    
    @patch('bioops.ensure_directories')
    @patch('bioops.RESULTS_DIR', Path('/tmp/results'))
    def test_parse_pdb_file_success(self, mock_ensure):
        """Test successful PDB file parsing"""
        pdb_content = """ATOM      1  N   MET A   1      20.154  29.699   5.276  1.00 49.05           N
ATOM      2  CA  MET A   1      21.260  30.420   5.898  1.00 49.05           C
END"""
        
        with patch('bioops.PDBIO') as mock_pdbio:
            mock_io_instance = Mock()
            mock_pdbio.return_value = mock_io_instance
            
            result = bioops.parse_pdb_file(pdb_content, "A", "test")
            
            assert result == True
            mock_io_instance.save.assert_called_once()


class TestSessionState:
    """Test session state management"""
    
    @patch.dict('streamlit.session_state', {})
    def test_initialize_session_state(self):
        """Test session state initialization"""
        bioops.initialize_session_state()
        
        assert 'pdb_1' in st.session_state
        assert 'pdb_2' in st.session_state
        assert 'processing_complete' in st.session_state
        assert 'processing_error' in st.session_state
        assert st.session_state['processing_complete'] == False
    
    @patch.dict('streamlit.session_state', {
        'processing_complete': True,
        'pdb_1': 'test',
        'pocket_similarity_matrix': pd.DataFrame({'a': [1, 2]})
    })
    def test_reset_session_state(self):
        """Test session state reset"""
        bioops.reset_session_state()
        
        assert st.session_state['processing_complete'] == False
        assert st.session_state['pdb_1'] is None
        assert st.session_state['pocket_similarity_matrix'].empty


class TestProcessingLogic:
    """Test main processing functions"""
    
    def test_validate_inputs_pdb_id_valid(self):
        """Test input validation for PDB ID mode"""
        with patch('bioops.validate_pdb_id', return_value=True):
            with patch('bioops.validate_chain_id', return_value=True):
                result = bioops.validate_inputs(
                    "**Enter PDB ID**", 
                    ("5CXV", "5DSG", "A", "B")
                )
                assert result == True
    
    def test_validate_inputs_pdb_id_invalid(self):
        """Test input validation for invalid PDB ID"""
        with patch('bioops.validate_pdb_id', side_effect=[False, True]):
            with patch('bioops.show_error') as mock_error:
                result = bioops.validate_inputs(
                    "**Enter PDB ID**",
                    ("BAD", "5DSG", "A", "B")
                )
                assert result == False
                mock_error.assert_called_once()
    
    def test_validate_inputs_upload_valid(self):
        """Test input validation for file upload mode"""
        mock_file1 = Mock()
        mock_file1.name = "protein1.pdb"
        mock_file1.size = 1024
        
        mock_file2 = Mock()
        mock_file2.name = "protein2.pdb"
        mock_file2.size = 1024
        
        with patch('bioops.validate_uploaded_file', return_value=True):
            with patch('bioops.validate_chain_id', return_value=True):
                result = bioops.validate_inputs(
                    "**Upload Files**",
                    (mock_file1, mock_file2, "A", "B")
                )
                assert result == True


class TestVisualization:
    """Test visualization functions"""
    
    @patch('py3Dmol.view')
    @patch('stmol.showmol')
    def test_view_pdb_success(self, mock_showmol, mock_view):
        """Test successful PDB visualization"""
        mock_view_instance = Mock()
        mock_view.return_value = mock_view_instance
        
        with patch('bioops.protFuncs.read_saved_pdb_file_to_view', return_value="PDB DATA"):
            bioops.view_pdb(["5CXV"], "5CXV", 0.5, False)
            
            mock_view_instance.addModel.assert_called_once_with("PDB DATA", "pdb")
            mock_view_instance.setStyle.assert_called_once()
            mock_view_instance.zoomTo.assert_called_once()
            mock_showmol.assert_called_once()
    
    @patch('streamlit.info')
    def test_view_pdb_no_selection(self, mock_info):
        """Test PDB visualization with no selection"""
        bioops.view_pdb([], "5CXV", 0.5, False)
        mock_info.assert_called_once()


class TestAnalysisExecution:
    """Test analysis execution functions"""
    
    @patch('bioops.protFuncs')
    @patch('streamlit.progress')
    @patch('streamlit.empty')
    @patch.dict('streamlit.session_state', {
        'pdb_1': 'test1',
        'pdb_2': 'test2',
        'pdb_1_id': '5CXV',
        'pdb_2_id': '5DSG',
        'input_file_option': '**Enter PDB ID**'
    })
    def test_run_analysis_success(self, mock_empty, mock_progress, mock_protfuncs):
        """Test successful analysis run"""
        # Mock all the analysis functions
        mock_protfuncs.FATCAT.return_value = True
        mock_protfuncs.getSuperimposed3Dpart.return_value = True
        mock_protfuncs.runP2rank.return_value = True
        mock_protfuncs.extractPocket_pdb.return_value = (3, 4)
        mock_protfuncs.runApoc.return_value = True
        mock_protfuncs.parse_pocket_comparison_result.return_value = True
        mock_protfuncs.compare_pockets_potentials.return_value = pd.DataFrame()
        
        # Mock progress bar
        mock_progress_bar = Mock()
        mock_progress.return_value = mock_progress_bar
        
        # Mock status text
        mock_status = Mock()
        mock_empty.return_value = mock_status
        
        result = bioops.run_analysis()
        
        assert result == True
        assert st.session_state.pdb_1_num_pocket == 3
        assert st.session_state.pdb_2_num_pocket == 4
        assert st.session_state.processing_complete == True
    
    @patch('bioops.protFuncs')
    @patch.dict('streamlit.session_state', {
        'pdb_1': None,
        'pdb_2': None
    })
    def test_run_analysis_missing_data(self, mock_protfuncs):
        """Test analysis with missing data"""
        result = bioops.run_analysis()
        
        assert result == False
        assert st.session_state.processing_error == True


class TestIntegration:
    """Integration tests"""
    
    @patch('streamlit.sidebar')
    @patch('streamlit.columns')
    @patch('streamlit.button')
    def test_render_input_ui(self, mock_button, mock_columns, mock_sidebar):
        """Test input UI rendering"""
        # Mock sidebar methods
        mock_sidebar.radio.return_value = "**Enter PDB ID**"
        mock_sidebar.columns.return_value = [Mock(), Mock()]
        
        bioops.render_input_ui()
        
        mock_sidebar.radio.assert_called_once()
        assert mock_sidebar.columns.call_count >= 2
        mock_sidebar.button.assert_called_once()
    
    @patch('streamlit.title')
    @patch('streamlit.tabs')
    def test_render_results_ui(self, mock_tabs, mock_title):
        """Test results UI rendering"""
        mock_tabs.return_value = [Mock(), Mock(), Mock()]
        
        with patch('bioops.render_structure_analysis_tab'):
            with patch('bioops.render_electrostatic_analysis_tab'):
                with patch('bioops.render_summary_tab'):
                    bioops.render_results_ui()
        
        mock_title.assert_called_once()
        mock_tabs.assert_called_once_with(["Structure Analysis", "Electrostatic Analysis", "Summary"])


# Fixtures for testing
@pytest.fixture
def mock_streamlit():
    """Mock streamlit for testing"""
    with patch.multiple('streamlit',
                       error=Mock(),
                       warning=Mock(),
                       info=Mock(),
                       success=Mock(),
                       button=Mock(return_value=False),
                       text_input=Mock(),
                       file_uploader=Mock(),
                       radio=Mock(),
                       sidebar=Mock(),
                       columns=Mock(return_value=[Mock(), Mock()]),
                       progress=Mock(),
                       empty=Mock()):
        yield


@pytest.fixture
def temp_results_dir():
    """Create temporary results directory for testing"""
    with tempfile.TemporaryDirectory() as tmpdir:
        results_dir = Path(tmpdir) / "results"
        results_dir.mkdir()
        yield results_dir