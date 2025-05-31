"""
Pytest configuration and shared fixtures for BioOps tests
"""

import pytest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock, patch
import pandas as pd
import numpy as np
import streamlit as st
import logging
import sys

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))


# Fixtures for temporary directories
@pytest.fixture
def temp_dir():
    """Provide a temporary directory that's automatically cleaned up"""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def bioops_test_env(temp_dir):
    """Create a complete BioOps test environment"""
    # Create directory structure
    (temp_dir / "results" / "fatcat").mkdir(parents=True)
    (temp_dir / "results" / "p2rank").mkdir(parents=True)
    (temp_dir / "results" / "apoc").mkdir(parents=True)
    (temp_dir / "results" / "apbs").mkdir(parents=True)
    (temp_dir / "results" / "overlap_score").mkdir(parents=True)
    (temp_dir / "logs").mkdir(parents=True)
    
    # Create mock config files
    (temp_dir / "config").mkdir()
    (temp_dir / "config" / "logging.yaml").write_text("""
version: 1
handlers:
  console:
    class: logging.StreamHandler
    level: DEBUG
loggers:
  bioops:
    level: DEBUG
    handlers: [console]
""")
    
    yield temp_dir


# Mock data fixtures
@pytest.fixture
def mock_pdb_content():
    """Provide mock PDB file content"""
    return """HEADER    TRANSFERASE                             20-JAN-99   1ABC
ATOM      1  N   MET A   1      20.154  29.699   5.276  1.00 49.05           N
ATOM      2  CA  MET A   1      21.260  30.420   5.898  1.00 49.05           C
ATOM      3  C   MET A   1      21.167  30.779   7.650  1.00 49.05           C
ATOM      4  O   MET A   1      20.145  30.639   8.160  1.00 49.05           O
ATOM      5  CB  MET A   1      21.362  31.652   5.088  1.00 49.05           C
ATOM      6  CG  MET A   1      22.247  32.693   5.691  1.00 49.05           C
ATOM      7  SD  MET A   1      23.824  33.191   4.980  1.00 49.05           S
ATOM      8  CE  MET A   1      24.359  34.235   6.336  1.00 49.05           C
ATOM      9  N   ALA A   2      22.154  31.699   6.276  1.00 49.05           N
ATOM     10  CA  ALA A   2      23.260  32.420   6.898  1.00 49.05           C
ATOM     11  C   ALA A   2      23.167  32.779   8.650  1.00 49.05           C
ATOM     12  O   ALA A   2      22.145  32.639   9.160  1.00 49.05           O
ATOM     13  CB  ALA A   2      23.362  33.652   6.088  1.00 49.05           C
TER
END"""


@pytest.fixture
def mock_alignment_content():
    """Provide mock FATCAT alignment content"""
    return """FATCAT result
Chain 1:    1  MKVLWAALLVTFLAGCQAKVEQAVE
                ||||||||||||||||||||||||
Chain 2:    1  MKVLWAALLVTFLAGCQAKVEQAVE

Chain 1:   26  TLVQSPATLSVSPGERVTLSCRAS
                |||||||||||||||||||||||||
Chain 2:   26  TLVQSPATLSVSPGERVTLSCRAS
"""


@pytest.fixture
def mock_p2rank_csv():
    """Provide mock P2Rank prediction CSV content"""
    return """rank,name,score,probability,sas_points,surf_atoms,center_x,center_y,center_z,residue_ids
1,pocket1,15.234,0.892,120,85,12.5,23.4,34.5,A_1 A_2 A_3 A_4 A_5
2,pocket2,12.567,0.756,98,72,45.2,12.3,23.4,A_10 A_11 A_12 A_13
3,pocket3,8.234,0.523,65,45,34.5,45.6,12.3,A_20 A_21 A_22
"""


@pytest.fixture
def mock_apoc_result():
    """Provide mock APoc result content"""
    return """
>>>>>>>>>>>>>>>>>>>>>>>>> Pocket 1: 5CXV_pkt_1 <-> 5DSG_pkt_1
Pocket: 5CXV_pkt_1
Pocket: 5DSG_pkt_1
PS-score = 0.875   P-value = 2.34e-8
Number of aligned residues = 35
RMSD = 1.234   Seq identity = 0.743

>>>>>>>>>>>>>>>>>>>>>>>>> Pocket 2: 5CXV_pkt_1 <-> 5DSG_pkt_2
Pocket: 5CXV_pkt_1
Pocket: 5DSG_pkt_2
PS-score = 0.623   P-value = 5.67e-4
Number of aligned residues = 25
RMSD = 2.345   Seq identity = 0.520
"""


@pytest.fixture
def mock_dx_content():
    """Provide mock DX file content"""
    return """object 1 class gridpositions counts 97 97 97
origin 0.000 0.000 0.000
delta 1.000 0.000 0.000
delta 0.000 1.000 0.000
delta 0.000 0.000 1.000
object 2 class gridconnections counts 97 97 97
object 3 class array type double rank 0 items 912673 data follows
0.123 0.234 -0.345 0.456 -0.567 0.678
0.789 -0.890 0.901 -0.012 0.123 0.234
attribute "dep" string "positions"
object "regular positions regular connections" class field
component "positions" value 1
component "connections" value 2
component "data" value 3
"""


@pytest.fixture
def mock_similarity_matrix():
    """Provide mock similarity matrix"""
    data = np.array([
        [1.0, 0.85, 0.62, 0.45, 0.33],
        [0.85, 1.0, 0.73, 0.52, 0.41],
        [0.62, 0.73, 1.0, 0.68, 0.55],
        [0.45, 0.52, 0.68, 1.0, 0.77],
        [0.33, 0.41, 0.55, 0.77, 1.0]
    ])
    
    index = ['5CXV_pkt_1', '5CXV_pkt_2', '5CXV_pkt_3', '5DSG_pkt_1', '5DSG_pkt_2']
    return pd.DataFrame(data, index=index, columns=index)


# Streamlit mocking
@pytest.fixture
def mock_streamlit(monkeypatch):
    """Mock streamlit components"""
    mock_st = Mock()
    
    # Mock session state
    mock_session_state = {}
    mock_st.session_state = mock_session_state
    
    # Mock UI components
    mock_st.error = Mock()
    mock_st.warning = Mock()
    mock_st.info = Mock()
    mock_st.success = Mock()
    mock_st.spinner = Mock(side_effect=lambda x: Mock())
    mock_st.progress = Mock(return_value=Mock())
    mock_st.empty = Mock(return_value=Mock())
    mock_st.button = Mock(return_value=False)
    mock_st.text_input = Mock()
    mock_st.file_uploader = Mock()
    mock_st.radio = Mock()
    mock_st.columns = Mock(return_value=[Mock(), Mock()])
    
    # Mock sidebar
    mock_sidebar = Mock()
    mock_sidebar.radio = Mock()
    mock_sidebar.button = Mock(return_value=False)
    mock_sidebar.columns = Mock(return_value=[Mock(), Mock()])
    mock_st.sidebar = mock_sidebar
    
    monkeypatch.setattr('streamlit', mock_st)
    return mock_st


# External tool mocking
@pytest.fixture
def mock_external_tools(monkeypatch):
    """Mock external bioinformatics tools"""
    def mock_run_command(cmd, *args, **kwargs):
        """Mock successful execution of external tools"""
        if 'FATCAT' in cmd:
            return (1, "FATCAT success", "")  # FATCAT returns 1 on success
        elif 'prank' in cmd:
            return (0, "P2Rank success", "")
        elif 'apoc' in cmd:
            return (0, "APoc success", "")
        elif 'pdb2pqr' in cmd:
            return (0, "PDB2PQR success", "")
        elif 'apbs' in cmd:
            return (0, "APBS success", "")
        else:
            return (0, "Success", "")
    
    monkeypatch.setattr('protFuncs.run_command', mock_run_command)


# Logging setup
@pytest.fixture(autouse=True)
def setup_logging():
    """Setup logging for tests"""
    logging.basicConfig(level=logging.DEBUG)
    yield
    # Clean up handlers
    logger = logging.getLogger()
    handlers = logger.handlers[:]
    for handler in handlers:
        logger.removeHandler(handler)


# Session state cleanup
@pytest.fixture
def clean_session_state():
    """Clean streamlit session state"""
    # Clear before test
    if hasattr(st, 'session_state'):
        for key in list(st.session_state.keys()):
            del st.session_state[key]
    
    yield
    
    # Clear after test
    if hasattr(st, 'session_state'):
        for key in list(st.session_state.keys()):
            del st.session_state[key]


# Network mocking
@pytest.fixture
def mock_requests(monkeypatch):
    """Mock requests library for PDB downloads"""
    mock_response = Mock()
    mock_response.text = "MOCK PDB CONTENT"
    mock_response.status_code = 200
    mock_response.raise_for_status = Mock()
    
    mock_get = Mock(return_value=mock_response)
    monkeypatch.setattr('requests.get', mock_get)
    
    return mock_get


# File system helpers
@pytest.fixture
def create_test_files(temp_dir):
    """Helper to create test files"""
    def _create_files(files_dict):
        """Create files from a dictionary of path: content"""
        created_files = {}
        for path, content in files_dict.items():
            full_path = temp_dir / path
            full_path.parent.mkdir(parents=True, exist_ok=True)
            full_path.write_text(content)
            created_files[path] = full_path
        return created_files
    
    return _create_files


# Pytest markers
def pytest_configure(config):
    """Register custom markers"""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "integration: marks tests as integration tests"
    )
    config.addinivalue_line(
        "markers", "unit: marks tests as unit tests"
    )
    config.addinivalue_line(
        "markers", "requires_tools: marks tests that require external tools"
    )


# Skip tests if tools not available
@pytest.fixture
def skip_if_tools_missing():
    """Skip test if required external tools are not available"""
    import shutil
    
    required_tools = ['FATCAT', 'prank', 'apoc', 'apbs', 'pdb2pqr']
    missing_tools = []
    
    for tool in required_tools:
        if not shutil.which(tool):
            missing_tools.append(tool)
    
    if missing_tools:
        pytest.skip(f"Required tools not found: {', '.join(missing_tools)}")


# Performance testing helpers
@pytest.fixture
def benchmark_timer():
    """Simple timer for performance benchmarking"""
    import time
    
    class Timer:
        def __init__(self):
            self.times = {}
        
        def start(self, name):
            self.times[name] = {'start': time.time()}
        
        def stop(self, name):
            if name in self.times:
                self.times[name]['end'] = time.time()
                self.times[name]['duration'] = self.times[name]['end'] - self.times[name]['start']
        
        def get_duration(self, name):
            return self.times.get(name, {}).get('duration', 0)
    
    return Timer()


# Test data generators
@pytest.fixture
def generate_test_pdb():
    """Generate test PDB content with specified parameters"""
    def _generate(num_residues=10, chain='A'):
        lines = ["HEADER    TEST PROTEIN                            01-JAN-00   TEST"]
        atom_num = 1
        
        for res_num in range(1, num_residues + 1):
            # Add backbone atoms
            for atom_name, x, y, z in [
                ('N', 20.0 + res_num, 30.0, 5.0),
                ('CA', 21.0 + res_num, 30.5, 5.5),
                ('C', 21.0 + res_num, 30.5, 6.5),
                ('O', 20.5 + res_num, 30.5, 7.0)
            ]:
                lines.append(
                    f"ATOM  {atom_num:5d}  {atom_name:<3s} ALA {chain}{res_num:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 50.00           {atom_name[0]}"
                )
                atom_num += 1
        
        lines.append("TER")
        lines.append("END")
        return "\n".join(lines)
    
    return _generate