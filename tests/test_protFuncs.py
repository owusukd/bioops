"""
Tests for protFuncs.py - Core analysis functions
"""

import pytest
from unittest.mock import Mock, patch, MagicMock, mock_open
import tempfile
import shutil
from pathlib import Path
import pandas as pd
import numpy as np
import subprocess

# Import the module to test
import sys
sys.path.insert(0, '..')
import protFuncs


class TestUtilityFunctions:
    """Test utility functions"""
    
    def test_ensure_directory(self):
        """Test directory creation"""
        with tempfile.TemporaryDirectory() as tmpdir:
            test_path = Path(tmpdir) / "test_dir"
            result = protFuncs.ensure_directory(test_path)
            
            assert test_path.exists()
            assert test_path.is_dir()
            assert result == test_path
    
    @patch('subprocess.run')
    def test_run_command_success(self, mock_run):
        """Test successful command execution"""
        mock_result = Mock()
        mock_result.returncode = 0
        mock_result.stdout = "Success"
        mock_result.stderr = ""
        mock_run.return_value = mock_result
        
        returncode, stdout, stderr = protFuncs.run_command("echo test")
        
        assert returncode == 0
        assert stdout == "Success"
        assert stderr == ""
        mock_run.assert_called_once()
    
    @patch('subprocess.run')
    def test_run_command_timeout(self, mock_run):
        """Test command timeout"""
        mock_run.side_effect = subprocess.TimeoutExpired("cmd", 30)
        
        with pytest.raises(protFuncs.ExternalToolError):
            protFuncs.run_command("long_running_command", timeout=30)
    
    def test_validate_file_exists_success(self):
        """Test file validation for existing file"""
        with tempfile.NamedTemporaryFile() as tmp:
            tmp.write(b"test content")
            tmp.flush()
            
            result = protFuncs.validate_file_exists(tmp.name)
            assert result.exists()
    
    def test_validate_file_exists_not_found(self):
        """Test file validation for non-existent file"""
        with pytest.raises(protFuncs.FileNotFoundError):
            protFuncs.validate_file_exists("/nonexistent/file.pdb")
    
    def test_validate_file_exists_empty(self):
        """Test file validation for empty file"""
        with tempfile.NamedTemporaryFile() as tmp:
            # Empty file
            with pytest.raises(protFuncs.FileNotFoundError):
                protFuncs.validate_file_exists(tmp.name)


class TestCoreAnalysisFunctions:
    """Test core analysis functions"""
    
    @patch('protFuncs.validate_file_exists')
    @patch('protFuncs.run_command')
    @patch('protFuncs.ensure_directory')
    def test_FATCAT_success(self, mock_ensure, mock_run, mock_validate):
        """Test successful FATCAT execution"""
        mock_ensure.return_value = Path("/tmp/fatcat")
        mock_validate.return_value = Path("/tmp/test.pdb")
        mock_run.return_value = (1, "Success", "")  # FATCAT returns 1 on success
        
        result = protFuncs.FATCAT("5CXV", "5DSG")
        
        assert result == True
        assert mock_run.call_count == 1
        assert mock_validate.call_count == 2
    
    def test_get_residue_id(self):
        """Test residue ID extraction from PDB"""
        pdb_string = """
ATOM      1  N   MET A   1      20.154  29.699   5.276  1.00 49.05           N
ATOM      2  CA  MET A   1      21.260  30.420   5.898  1.00 49.05           C
ATOM      3  C   MET A   1      21.167  30.779   7.650  1.00 49.05           C
ATOM      4  N   ALA A   2      22.154  31.699   6.276  1.00 49.05           N
ATOM      5  CA  ALA A   2      23.260  32.420   6.898  1.00 49.05           C
"""
        residues = protFuncs.get_residue_id(pdb_string)
        
        assert residues == [1, 2]
        assert len(residues) == 2
    
    @patch('protFuncs.validate_file_exists')
    @patch('builtins.open', new_callable=mock_open)
    def test_getSuperimposed3Dpart(self, mock_file, mock_validate):
        """Test superimposed 3D part extraction"""
        # Mock alignment file content
        aln_content = """
Chain 1:    1  MKVL
               ||||
Chain 2:    1  MKVL
"""
        # Mock PDB content
        pdb_content = """
ATOM      1  N   MET A   1      20.154  29.699   5.276  1.00 49.05           N
ATOM      2  CA  MET A   1      21.260  30.420   5.898  1.00 49.05           C
"""
        
        mock_validate.return_value = Path("/tmp/test")
        mock_file.return_value.readlines.return_value = aln_content.split('\n')
        mock_file.return_value.read.return_value = pdb_content
        
        with patch('protFuncs.get_residue_id', return_value=[1, 2]):
            with patch('protFuncs.write_flexible_pdb'):
                result = protFuncs.getSuperimposed3Dpart("5CXV", "5DSG")
                assert result == True
    
    @patch('protFuncs.validate_file_exists')
    @patch('protFuncs.run_command')
    @patch('protFuncs.ensure_directory')
    def test_runP2rank_success(self, mock_ensure, mock_run, mock_validate):
        """Test successful P2Rank execution"""
        mock_ensure.return_value = Path("/tmp/p2rank")
        mock_validate.return_value = Path("/tmp/test.pdb")
        mock_run.return_value = (0, "Success", "")
        
        # Mock CSV file check
        with patch('pathlib.Path.exists', return_value=True):
            result = protFuncs.runP2rank("5CXV", "5DSG")
        
        assert result == True
        assert mock_run.call_count == 2  # Called for both proteins
    
    def test_write_pocket_pdb(self):
        """Test pocket PDB writing"""
        pocket_predictions = [
            "1,pocket1,10.5,5.5,3.2,100,0.95,0.85,10,A_1 A_2 A_3"
        ]
        
        pdb_lines = [
            "ATOM      1  N   MET A   1      20.154  29.699   5.276  1.00 49.05           N\n",
            "ATOM      2  CA  MET A   1      21.260  30.420   5.898  1.00 49.05           C\n",
            "ATOM      3  N   ALA A   2      22.154  31.699   6.276  1.00 49.05           N\n",
            "END\n"
        ]
        
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch('protFuncs.validate_file_exists'):
                with patch('protFuncs.ensure_directory', return_value=Path(tmpdir)):
                    with patch('builtins.open', mock_open()) as mock_file:
                        mock_file.return_value.readlines.return_value = pdb_lines
                        
                        result = protFuncs.write_pocket_pdb("test", pocket_predictions, 1)
                        
                        assert result == True
                        assert mock_file.call_count >= 2  # Read and write


class TestElectrostaticAnalysis:
    """Test electrostatic analysis functions"""
    
    @patch('protFuncs.validate_file_exists')
    @patch('protFuncs.run_command')
    @patch('protFuncs.ensure_directory')
    def test_run_pockets_pdb2pqr(self, mock_ensure, mock_run, mock_validate):
        """Test PDB2PQR execution"""
        mock_ensure.return_value = Path("/tmp/apbs/test")
        mock_validate.return_value = Path("/tmp/test.pdb")
        mock_run.return_value = (0, "Success", "")
        
        with patch('pathlib.Path.exists', return_value=True):
            result = protFuncs.run_pockets_pdb2pqr("test", 3)
        
        assert result == True
        assert mock_run.call_count == 3  # Once for each pocket
    
    @patch('protFuncs.validate_file_exists')
    @patch('protFuncs.run_command')
    @patch('protFuncs.ensure_directory')
    def test_run_pockets_apbs(self, mock_ensure, mock_run, mock_validate):
        """Test APBS execution"""
        mock_ensure.return_value = Path("/tmp/apbs/test")
        mock_run.return_value = (0, "Success", "")
        
        with patch('pathlib.Path.exists', return_value=True):
            with patch('builtins.open', mock_open()):
                result = protFuncs.run_pockets_apbs("test", 2)
        
        assert result == True
        assert mock_run.call_count == 2  # Once for each pocket
    
    @patch('protFuncs.parse_dx_file')
    @patch('pathlib.Path.exists')
    def test_compare_pockets_potentials(self, mock_exists, mock_parse):
        """Test pocket potential comparison"""
        mock_exists.return_value = True
        mock_parse.return_value = np.random.rand(10, 10, 10)
        
        result = protFuncs.compare_pockets_potentials("pdb1", "pdb2", 2, 2)
        
        assert isinstance(result, pd.DataFrame)
        assert result.shape == (4, 4)  # 2+2 pockets total


class TestVisualizationFunctions:
    """Test visualization functions"""
    
    def test_heatmap_pocket_corr(self):
        """Test heatmap generation"""
        # Create test similarity matrix
        data = np.random.rand(5, 5)
        similarity_matrix = pd.DataFrame(data, 
                                       index=[f"pocket_{i}" for i in range(5)],
                                       columns=[f"pocket_{i}" for i in range(5)])
        
        with patch('matplotlib.pyplot.savefig'):
            fig = protFuncs.heatmap_pocket_corr(similarity_matrix, "pdb1", "pdb2")
            
            assert fig is not None
            assert hasattr(fig, 'axes')
    
    def test_cluster_pockets(self):
        """Test dendrogram generation"""
        # Create test similarity matrix
        data = np.array([[1.0, 0.8, 0.3],
                        [0.8, 1.0, 0.5],
                        [0.3, 0.5, 1.0]])
        similarity_matrix = pd.DataFrame(data,
                                       index=["p1", "p2", "p3"],
                                       columns=["p1", "p2", "p3"])
        
        with patch('matplotlib.pyplot.savefig'):
            fig = protFuncs.cluster_pockets(similarity_matrix, "pdb1", "pdb2")
            
            assert fig is not None


class TestDataExportFunctions:
    """Test data export functions"""
    
    def test_export_analysis_summary(self):
        """Test analysis summary export"""
        summary = protFuncs.export_analysis_summary("pdb1", "pdb2", 3, 4)
        
        assert "proteins" in summary
        assert summary["proteins"]["protein_1"]["id"] == "pdb1"
        assert summary["proteins"]["protein_1"]["num_pockets"] == 3
        assert summary["proteins"]["protein_2"]["id"] == "pdb2"
        assert summary["proteins"]["protein_2"]["num_pockets"] == 4
        assert "analysis_folders" in summary
        assert "timestamp" in summary
    
    @patch('protFuncs.RESULTS_BASE_DIR', Path("/tmp/results"))
    @patch('pathlib.Path.exists')
    @patch('zipfile.ZipFile')
    def test_create_results_package(self, mock_zip, mock_exists, mock_streamlit):
        """Test results package creation"""
        mock_exists.return_value = True
        
        # Mock the ZipFile context manager
        mock_zip_instance = MagicMock()
        mock_zip.return_value.__enter__.return_value = mock_zip_instance
        
        with patch('protFuncs.add_folder_to_zip', return_value=True):
            with patch('protFuncs.export_analysis_summary', return_value={"test": "data"}):
                result = protFuncs.create_results_package("pdb1", "pdb2", 3, 4)
        
        assert result is not None
        mock_zip_instance.writestr.assert_called_once()  # For summary JSON


class TestAPocFunctions:
    """Test APoc-related functions"""
    
    @patch('protFuncs.validate_file_exists')
    @patch('protFuncs.run_command')
    def test_runApoc_success(self, mock_run, mock_validate):
        """Test successful APoc execution"""
        mock_validate.return_value = Path("/tmp/test.pdb")
        mock_run.return_value = (0, "Success", "")
        
        with patch('pathlib.Path.exists', return_value=True):
            with patch('pathlib.Path.stat') as mock_stat:
                mock_stat.return_value.st_size = 100
                result = protFuncs.runApoc("pdb1", "pdb2")
        
        assert result == True
    
    @patch('protFuncs.validate_file_exists')
    @patch('builtins.open', mock_open())
    def test_parse_pocket_comparison_result(self, mock_validate):
        """Test APoc result parsing"""
        mock_content = """
>>>>>>>>>>>>>>>>>>>>>>>>> Pocket 1 vs Pocket 2
Pocket: pdb1_pkt_1
Pocket: pdb2_pkt_1
PS-score = 0.85  P-value = 1.2e-5
Number of aligned residues = 25
RMSD = 1.23  Seq identity = 0.68
"""
        
        mock_validate.return_value = Path("/tmp/test.txt")
        
        with patch('builtins.open', mock_open(read_data=mock_content)):
            result = protFuncs.parse_pocket_comparison_result("pdb1", "pdb2")
        
        assert result == True


# Test fixtures
@pytest.fixture
def mock_streamlit():
    """Mock streamlit for testing"""
    with patch.multiple('streamlit',
                       spinner=Mock(side_effect=lambda x: MagicMock()),
                       success=Mock(),
                       error=Mock()):
        yield


@pytest.fixture
def temp_analysis_dir():
    """Create temporary analysis directory structure"""
    with tempfile.TemporaryDirectory() as tmpdir:
        base = Path(tmpdir)
        (base / "results" / "fatcat").mkdir(parents=True)
        (base / "results" / "p2rank").mkdir(parents=True)
        (base / "results" / "apoc").mkdir(parents=True)
        (base / "results" / "apbs").mkdir(parents=True)
        (base / "results" / "overlap_score").mkdir(parents=True)
        yield base