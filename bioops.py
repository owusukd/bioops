# Import enhanced logging
from logging_config import (bioops_logger, log_analysis_complete,
                            log_analysis_start, log_analysis_step, log_context,
                            log_performance)

# Get logger instance
logger = bioops_logger.get_logger(__name__)

import os
import re
from datetime import datetime
from io import StringIO
from pathlib import Path
from typing import List

import pandas as pd
import py3Dmol
import requests
import streamlit as st
from Bio.PDB import PDBIO, PDBParser, Select
from stmol import showmol

import protFuncs

# Constants
TEMP_DIR = Path("tempDownloadDir")
RESULTS_DIR = TEMP_DIR / "results"
DEFAULT_CHAIN = "A"
PDB_URL_TEMPLATE = "https://files.rcsb.org/download/{}.pdb"

# Streamlit Configuration
st.set_page_config(
    layout="wide",
    page_title="BioOps - Protein Analysis",
    initial_sidebar_state="expanded",
    page_icon="logo/BioOps_Logo_icon.png",
)

# Log application startup
logger.info(
    "BioOps application started",
    extra={
        "event": "app_startup",
        "streamlit_version": st.__version__,
        "temp_dir": str(TEMP_DIR),
        "results_dir": str(RESULTS_DIR),
    },
)


# Configuration and Setup
@st.cache_data
def get_app_config():
    """Get application configuration"""
    return {
        "logo_main": "logo/BioOps_Main_Logo.png",
        "logo_icon": "logo/BioOps_Logo_icon.png",
        "temp_dir": TEMP_DIR,
        "results_dir": RESULTS_DIR,
    }


@log_performance()
def ensure_directories():
    """Ensure required directories exist"""
    config = get_app_config()
    config["temp_dir"].mkdir(parents=True, exist_ok=True)
    config["results_dir"].mkdir(parents=True, exist_ok=True)
    logger.debug(f"Directories ensured: {config['temp_dir']}, {config['results_dir']}")


# Error Handling
class BioOpsError(Exception):
    """Base exception for BioOps application"""

    pass


class FileProcessingError(BioOpsError):
    """Error in file processing"""

    pass


class PDBDownloadError(BioOpsError):
    """Error in PDB file download"""

    pass


class AnalysisError(BioOpsError):
    """Error in protein analysis"""

    pass


def show_error(title: str, message: str, error_type: str = "error"):
    """Show error message with consistent styling"""
    logger.warning(
        f"Showing {error_type} to user",
        extra={
            "title": title,
            "error_message": message,  # Changed from 'message' to 'error_message'
            "error_type": error_type,
        },
    )

    if error_type == "error":
        st.error(f"**{title}**: {message}", icon="🚨")
    elif error_type == "warning":
        st.warning(f"**{title}**: {message}", icon="⚠️")
    else:
        st.info(f"**{title}**: {message}", icon="ℹ️")


@st.dialog("Input Error")
def input_file_error():
    logger.error("Input file error dialog shown")
    st.error("Please ensure you have:")
    st.markdown(
        """
    - Uploaded two valid PDB files, OR
    - Entered correct PDB IDs and Chain IDs
    - All required fields are filled
    """
    )


@st.dialog("File Processing Error")
def file_processing_error(details: str = ""):
    logger.error(f"File processing error dialog shown: {details}")
    st.error("Error processing PDB file!")
    if details:
        st.code(details)


@st.dialog("Analysis Error")
def analysis_error(details: str = ""):
    logger.error(f"Analysis error dialog shown: {details}")
    st.error("Error during protein analysis!")
    if details:
        st.code(details)


# Input Validation
@log_performance()
def validate_pdb_id(pdb_id: str) -> bool:
    """Validate PDB ID format"""
    if not pdb_id:
        logger.debug("Empty PDB ID provided")
        return False
    # PDB IDs are typically 4 characters (alphanumeric)
    is_valid = bool(re.match(r"^[a-zA-Z0-9]{4}$", pdb_id.strip()))
    logger.debug(f"PDB ID validation: {pdb_id} -> {is_valid}")
    return is_valid


def validate_chain_id(chain_id: str) -> bool:
    """Validate chain ID format"""
    if not chain_id:
        return True
    # Chain IDs are typically single characters
    is_valid = bool(re.match(r"^[a-zA-Z]$", chain_id.strip()))
    logger.debug(f"Chain ID validation: {chain_id} -> {is_valid}")
    return is_valid


def validate_uploaded_file(uploaded_file) -> bool:
    """Validate uploaded PDB file"""
    if uploaded_file is None:
        return False

    # Check file extension
    if not uploaded_file.name.lower().endswith(".pdb"):
        logger.warning(f"Invalid file extension: {uploaded_file.name}")
        return False

    # Check file size (reasonable limit)
    if uploaded_file.size > 50 * 1024 * 1024:  # 50MB limit
        logger.warning(
            f"File too large: {uploaded_file.name} ({uploaded_file.size} bytes)"
        )
        return False

    logger.debug(f"File validation passed: {uploaded_file.name}")
    return True


# UI Setup
def setup_page():
    """Setup page configuration and header"""
    config = get_app_config()

    # Check if logo files exist
    logo_main = config["logo_main"] if os.path.exists(config["logo_main"]) else None
    logo_icon = config["logo_icon"] if os.path.exists(config["logo_icon"]) else None

    if logo_main:
        st.logo(logo_main, icon_image=logo_icon)
    else:
        st.title("🧬 BioOps - Protein Analysis Tool")

    setup_sidebar()


def setup_sidebar():
    """Setup sidebar content"""
    st.sidebar.header("🔬 Protein Analysis")

    st.sidebar.markdown(
        """
    **BioOps offers comprehensive Protein analysis:**
        
    * **3D Structural Comparison** using FATCAT
    * **Binding Pocket Prediction** using P2Rank  
    * **Binding Pocket Comparison** using APoc
    * **3D Structure Similarity & Pocket Overlap Scores**
    * **Electrostatic Properties Comparison**
    """
    )

    st.sidebar.subheader("📁 Input Requirements")
    st.sidebar.markdown(
        """
    **Choose one of the following input methods:**

    * **PDB ID**: Enter 4-character PDB identifier
    * **Upload**: Upload PDB files from your computer
    
    **Note**: Chain ID is optional (defaults to 'A')
    """
    )


# PDB File Handling
class ChainSelect(Select):
    """PDB chain selector for BioPython"""

    def __init__(self, chain_id: str):
        self.chain_id = chain_id.upper()

    def accept_chain(self, chain):
        return chain.get_id() == self.chain_id


@st.cache_data(show_spinner=False)
@log_performance()
def download_pdb_content(pdb_id: str) -> str:
    """Download PDB file content with caching"""
    try:
        url = PDB_URL_TEMPLATE.format(pdb_id.lower())
        logger.info(
            f"Downloading PDB file: {pdb_id}", extra={"pdb_id": pdb_id, "url": url}
        )

        response = requests.get(url, timeout=30)
        response.raise_for_status()

        logger.info(
            f"PDB download successful: {pdb_id}",
            extra={"pdb_id": pdb_id, "content_length": len(response.text)},
        )
        return response.text

    except requests.RequestException as e:
        logger.error(
            f"PDB download failed: {pdb_id}", extra={"pdb_id": pdb_id, "error": str(e)}
        )
        raise PDBDownloadError(f"Failed to download PDB {pdb_id}: {str(e)}")


@log_analysis_step("Parse PDB File")
def parse_pdb_file(pdb_content: str, chain_id: str, pdb_id: str) -> bool:
    """Parse PDB content and save processed file"""
    try:
        # Ensure directories exist
        ensure_directories()

        parser = PDBParser(QUIET=True)
        pdb_io = StringIO(pdb_content)
        structure = parser.get_structure(pdb_id, pdb_io)

        # Save processed PDB file
        io_handler = PDBIO()
        io_handler.set_structure(structure)

        output_path = RESULTS_DIR / f"{pdb_id}.pdb"
        chain_selector = ChainSelect(chain_id or DEFAULT_CHAIN)
        io_handler.save(str(output_path), chain_selector)

        # Verify file was created and has content
        if not output_path.exists() or output_path.stat().st_size == 0:
            raise FileProcessingError(f"Failed to create valid PDB file for {pdb_id}")

        logger.info(
            f"PDB file parsed and saved: {pdb_id}",
            extra={
                "pdb_id": pdb_id,
                "chain_id": chain_id,
                "output_path": str(output_path),
                "file_size": output_path.stat().st_size,
            },
        )

        return True

    except Exception as e:
        logger.error(
            f"Error parsing PDB file {pdb_id}: {str(e)}",
            extra={
                "pdb_id": pdb_id,
                "chain_id": chain_id,
                "error_type": type(e).__name__,
            },
        )
        raise FileProcessingError(f"Unable to parse PDB file: {str(e)}")


@log_analysis_step("Process Uploaded File")
def process_uploaded_file(uploaded_file, chain_id: str, pdb_id: str) -> bool:
    """Process uploaded PDB file"""
    try:
        if not validate_uploaded_file(uploaded_file):
            raise FileProcessingError("Invalid file format or size")

        logger.info(
            f"Processing uploaded file: {uploaded_file.name}",
            extra={
                "filename": uploaded_file.name,
                "file_size": uploaded_file.size,
                "pdb_id": pdb_id,
            },
        )

        # Read file content
        content = uploaded_file.getvalue().decode("utf-8")
        return parse_pdb_file(content, chain_id, pdb_id)

    except UnicodeDecodeError:
        logger.error(f"File encoding error: {uploaded_file.name}")
        raise FileProcessingError(
            "File encoding error - ensure file is plain text PDB format"
        )
    except Exception as e:
        logger.error(
            f"Error processing uploaded file: {str(e)}",
            extra={"filename": uploaded_file.name, "error_type": type(e).__name__},
        )
        raise FileProcessingError(f"Unable to process uploaded file: {str(e)}")


@log_analysis_step("Process PDB ID")
def process_pdb_id(pdb_id: str, chain_id: str) -> bool:
    """Process PDB ID by downloading and parsing"""
    try:
        if not validate_pdb_id(pdb_id):
            raise FileProcessingError(f"Invalid PDB ID format: {pdb_id}")

        # Download PDB content
        content = download_pdb_content(pdb_id)
        return parse_pdb_file(content, chain_id, pdb_id)

    except Exception as e:
        logger.error(f"Error processing PDB ID {pdb_id}: {str(e)}")
        raise


# Session State Management
def initialize_session_state():
    """Initialize session state variables with defaults"""
    defaults = {
        "pdb_1": None,
        "pdb_2": None,
        "pdb_1_id": None,
        "pdb_2_id": None,
        "pdb_1_num_pocket": 0,
        "pdb_2_num_pocket": 0,
        "pocket_similarity_matrix": pd.DataFrame(),
        "checkbox_states_1": {},
        "checkbox_states_2": {},
        "processing_complete": False,
        "processing_error": False,
        "error_message": "",
        "analysis_progress": 0,
        "analysis_id": None,
    }

    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value

    logger.debug("Session state initialized")


def reset_session_state():
    """Reset session state for new analysis"""
    keys_to_reset = [
        "processing_complete",
        "processing_error",
        "error_message",
        "pdb_1",
        "pdb_2",
        "pdb_1_num_pocket",
        "pdb_2_num_pocket",
        "pocket_similarity_matrix",
        "checkbox_states_1",
        "checkbox_states_2",
        "analysis_progress",
    ]

    for key in keys_to_reset:
        if key in st.session_state:
            if isinstance(st.session_state[key], dict):
                st.session_state[key] = {}
            elif isinstance(st.session_state[key], pd.DataFrame):
                st.session_state[key] = pd.DataFrame()
            else:
                st.session_state[key] = 0 if "num" in key or "progress" in key else None

    logger.debug("Session state reset")


# Main Processing Logic
@log_performance()
def validate_inputs(input_option: str, args: tuple) -> bool:
    """Validate input parameters"""
    logger.debug(f"Validating inputs for option: {input_option}")

    if input_option == "**Enter PDB ID**":
        pdb1_id, pdb2_id, pdb1_chain, pdb2_chain = args

        # Validate PDB IDs
        if not validate_pdb_id(pdb1_id) or not validate_pdb_id(pdb2_id):
            show_error("Invalid Input", "Please enter valid 4-character PDB IDs")
            return False

        # Validate chain IDs
        if not validate_chain_id(pdb1_chain) or not validate_chain_id(pdb2_chain):
            show_error("Invalid Input", "Chain IDs must be single characters")
            return False

    elif input_option == "**Upload Files**":
        pdb1_file, pdb2_file, pdb1_chain, pdb2_chain = args

        # Validate uploaded files
        if not validate_uploaded_file(pdb1_file) or not validate_uploaded_file(
            pdb2_file
        ):
            show_error("Invalid Files", "Please upload valid PDB files (max 50MB each)")
            return False

        # Validate chain IDs
        if not validate_chain_id(pdb1_chain) or not validate_chain_id(pdb2_chain):
            show_error("Invalid Input", "Chain IDs must be single characters")
            return False

    return True


@log_analysis_step("Process Input Files")
def process_input_files(input_option: str, args: tuple) -> bool:
    """Process input files with improved error handling"""
    try:
        # Reset error state
        st.session_state.processing_error = False
        st.session_state.error_message = ""

        # Validate inputs first
        if not validate_inputs(input_option, args):
            st.session_state.processing_error = True
            return False

        if input_option == "**Enter PDB ID**":
            pdb1_id, pdb2_id, pdb1_chain, pdb2_chain = args

            if pdb1_id != pdb2_id:
                # Process PDB IDs
                with st.spinner(f"Downloading PDB file for {pdb1_id}..."):
                    if not process_pdb_id(pdb1_id, pdb1_chain):
                        raise FileProcessingError(f"Failed to process PDB {pdb1_id}")

                with st.spinner(f"Downloading PDB file for {pdb2_id}..."):
                    if not process_pdb_id(pdb2_id, pdb2_chain):
                        raise FileProcessingError(f"Failed to process PDB {pdb2_id}")

                st.session_state.pdb_1_id = pdb1_id.lower()
                st.session_state.pdb_2_id = pdb2_id.lower()
            else:
                raise FileProcessingError(
                    f"Enter two different PDB IDs: {pdb1_id} and {pdb2_id}"
                )

        elif input_option == "**Upload Files**":
            pdb1_file, pdb2_file, pdb1_chain, pdb2_chain = args

            # Extract file names without extension
            st.session_state.pdb_1_id = Path(pdb1_file.name).stem.lower()
            st.session_state.pdb_2_id = Path(pdb2_file.name).stem.lower()

            if st.session_state.pdb_1_id != st.session_state.pdb_2_id:
                # Process uploaded files
                with st.spinner(f"Processing uploaded file: {pdb1_file.name}..."):
                    if not process_uploaded_file(
                        pdb1_file, pdb1_chain, st.session_state.pdb_1_id
                    ):
                        raise FileProcessingError(
                            f"Failed to process uploaded file {pdb1_file.name}"
                        )

                with st.spinner(f"Processing uploaded file: {pdb2_file.name}..."):
                    if not process_uploaded_file(
                        pdb2_file, pdb2_chain, st.session_state.pdb_2_id
                    ):
                        raise FileProcessingError(
                            f"Failed to process uploaded file {pdb2_file.name}"
                        )

            else:
                raise FileProcessingError(
                    f"Upload two different PDB file: {st.session_state.pdb_1_id} and {st.session_state.pdb_2_id}"
                )

        # Read processed PDB files and store in session state
        pdb1_path = RESULTS_DIR / f"{st.session_state.pdb_1_id}.pdb"
        pdb2_path = RESULTS_DIR / f"{st.session_state.pdb_2_id}.pdb"

        with open(pdb1_path, "r") as f:
            st.session_state.pdb_1 = f.read()
        with open(pdb2_path, "r") as f:
            st.session_state.pdb_2 = f.read()

        logger.info(
            "Input files processed successfully",
            extra={
                "pdb_1_id": st.session_state.pdb_1_id,
                "pdb_2_id": st.session_state.pdb_2_id,
                "input_method": input_option,
            },
        )

        return True

    except (FileProcessingError, PDBDownloadError) as e:
        st.session_state.processing_error = True
        st.session_state.error_message = str(e)
        show_error("Processing Error", str(e))
        return False
    except Exception as e:
        st.session_state.processing_error = True
        st.session_state.error_message = f"Unexpected error: {str(e)}"
        show_error("Unexpected Error", str(e))
        logger.error(f"Unexpected error in process_input_files: {str(e)}")
        return False


@log_performance()
def run_analysis() -> bool:
    """Run protein analysis with progress tracking"""
    analysis_id = None
    try:
        # Verify PDB data is available
        required_data = [
            st.session_state.pdb_1,
            st.session_state.pdb_2,
            st.session_state.pdb_1_id,
            st.session_state.pdb_2_id,
        ]

        if any(x is None for x in required_data):
            raise AnalysisError("Required PDB data is missing")

        # Start analysis tracking
        analysis_id = log_analysis_start(
            st.session_state.pdb_1_id,
            st.session_state.pdb_2_id,
            st.session_state.input_file_option,
        )
        st.session_state.analysis_id = analysis_id

        # Set logging context for this analysis
        with log_context(
            analysis_id=analysis_id,
            pdb_1_id=st.session_state.pdb_1_id,
            pdb_2_id=st.session_state.pdb_2_id,
        ):
            # Progress tracking
            progress_bar = st.progress(0)
            status_text = st.empty()

            # Analysis steps with progress updates
            analysis_steps = [
                (
                    "Running FATCAT for 3D structural comparison",
                    lambda: protFuncs.FATCAT(
                        st.session_state.pdb_1_id, st.session_state.pdb_2_id
                    ),
                ),
                (
                    "Extracting superimposed 3D parts",
                    lambda: protFuncs.getSuperimposed3Dpart(
                        st.session_state.pdb_1_id, st.session_state.pdb_2_id
                    ),
                ),
                (
                    "Predicting ligand binding pockets with P2rank",
                    lambda: protFuncs.runP2rank(
                        st.session_state.pdb_1_id, st.session_state.pdb_2_id
                    ),
                ),
                (
                    "Extracting pocket PDB files",
                    lambda: protFuncs.extractPocket_pdb(
                        st.session_state.pdb_1_id, st.session_state.pdb_2_id
                    ),
                ),
                (
                    "Comparing binding pockets with APoc",
                    lambda: protFuncs.runApoc(
                        st.session_state.pdb_1_id, st.session_state.pdb_2_id
                    ),
                ),
                (
                    "Parsing pocket comparison results",
                    lambda: protFuncs.parse_pocket_comparison_result(
                        st.session_state.pdb_1_id, st.session_state.pdb_2_id
                    ),
                ),
                ("Calculating overlap scores", lambda: calculate_overlap_scores()),
                (
                    "Running electrostatic analysis",
                    lambda: run_electrostatic_analysis(),
                ),
            ]

            for i, (description, func) in enumerate(analysis_steps):
                status_text.text(description)
                progress = (i + 1) / len(analysis_steps)
                progress_bar.progress(progress)

                logger.info(
                    f"Starting analysis step: {description}",
                    extra={
                        "step_index": i,
                        "step_name": description,
                        "progress": progress,
                    },
                )

                try:
                    if i == 3:  # extractPocket_pdb returns values
                        result = func()
                        if isinstance(result, tuple) and len(result) == 2:
                            (
                                st.session_state.pdb_1_num_pocket,
                                st.session_state.pdb_2_num_pocket,
                            ) = result
                            logger.info(
                                f"Pockets extracted",
                                extra={
                                    "pdb_1_pockets": st.session_state.pdb_1_num_pocket,
                                    "pdb_2_pockets": st.session_state.pdb_2_num_pocket,
                                },
                            )
                    elif i == 7:  # electrostatic analysis returns similarity matrix
                        result = func()
                        if isinstance(result, pd.DataFrame):
                            st.session_state.pocket_similarity_matrix = result
                            logger.info(
                                f"Similarity matrix calculated",
                                extra={"matrix_shape": result.shape},
                            )
                    else:
                        func()
                except Exception as step_error:
                    logger.error(
                        f"Analysis step failed: {description}",
                        extra={
                            "step_name": description,
                            "error": str(step_error),
                            "error_type": type(step_error).__name__,
                        },
                    )
                    raise AnalysisError(
                        f"Failed at step '{description}': {str(step_error)}"
                    )

        # Clear progress indicators
        progress_bar.empty()
        status_text.empty()

        st.session_state.processing_complete = True
        st.success("✅ Analysis completed successfully!")

        # Log analysis completion
        log_analysis_complete(
            analysis_id, st.session_state.pdb_1_id, st.session_state.pdb_2_id, True
        )

        return True

    except AnalysisError as e:
        st.session_state.processing_error = True
        st.session_state.error_message = str(e)
        show_error("Analysis Error", str(e))

        if analysis_id:
            log_analysis_complete(
                analysis_id, st.session_state.pdb_1_id, st.session_state.pdb_2_id, False
            )

        return False
    except Exception as e:
        st.session_state.processing_error = True
        st.session_state.error_message = f"Unexpected analysis error: {str(e)}"
        show_error("Unexpected Error", str(e))
        logger.error(f"Unexpected error in run_analysis: {str(e)}")

        if analysis_id:
            log_analysis_complete(
                analysis_id, st.session_state.pdb_1_id, st.session_state.pdb_2_id, False
            )

        return False


@log_analysis_step("Calculate Overlap Scores")
def calculate_overlap_scores():
    """Calculate overlap scores for both proteins"""
    protFuncs.pocket_3D_similar_overlap_score(
        st.session_state.pdb_1_id,
        st.session_state.pdb_2_id,
        st.session_state.pdb_1_num_pocket,
    )
    protFuncs.pocket_3D_similar_overlap_score(
        st.session_state.pdb_2_id,
        st.session_state.pdb_1_id,
        st.session_state.pdb_2_num_pocket,
    )


@log_analysis_step("Electrostatic Analysis")
def run_electrostatic_analysis():
    """Run electrostatic analysis and return similarity matrix"""
    # PDB2PQR processing
    protFuncs.run_pockets_pdb2pqr(
        st.session_state.pdb_1_id, st.session_state.pdb_1_num_pocket
    )
    protFuncs.run_pockets_pdb2pqr(
        st.session_state.pdb_2_id, st.session_state.pdb_2_num_pocket
    )

    # APBS processing
    protFuncs.run_pockets_apbs(
        st.session_state.pdb_1_id, st.session_state.pdb_1_num_pocket
    )
    protFuncs.run_pockets_apbs(
        st.session_state.pdb_2_id, st.session_state.pdb_2_num_pocket
    )

    # Compare potentials
    similarity_matrix = protFuncs.compare_pockets_potentials(
        st.session_state.pdb_1_id,
        st.session_state.pdb_2_id,
        st.session_state.pdb_1_num_pocket,
        st.session_state.pdb_2_num_pocket,
    )

    return similarity_matrix


@log_performance()
def main():
    """Main processing function"""
    try:
        # Get input values based on selected option
        if st.session_state.input_file_option == "**Enter PDB ID**":
            args = (
                st.session_state.protein_1,
                st.session_state.protein_2,
                st.session_state.protein1_chain,
                st.session_state.protein2_chain,
            )
            args = tuple(arg.strip().lower() for arg in args)
        else:
            args = (
                st.session_state.protein_1_file,
                st.session_state.protein_2_file,
                st.session_state.file1_chain,
                st.session_state.file2_chain,
            )
            args = tuple(arg.strip().lower() for arg in args)
            
        logger.info(
            "Starting main analysis",
            extra={"input_option": st.session_state.input_file_option},
        )

        # Reset previous state
        reset_session_state()

        # Process input files
        if process_input_files(st.session_state.input_file_option, args):
            # Run analysis if file processing was successful
            run_analysis()

    except Exception as e:
        st.session_state.processing_error = True
        st.session_state.error_message = f"Main processing error: {str(e)}"
        show_error("Processing Error", str(e))
        logger.error(f"Error in main: {str(e)}")


# UI Components
def render_input_ui():
    """Render input UI components with improved validation"""
    st.sidebar.radio(
        "Input File Option",
        ["**Enter PDB ID**", "**Upload Files**"],
        key="input_file_option",
        help="Choose how to provide your protein structures",
    )

    # Create input columns
    col1, col2 = st.sidebar.columns(2, vertical_alignment="bottom")
    col3, col4 = st.sidebar.columns(2, vertical_alignment="bottom")

    if st.session_state.input_file_option == "**Enter PDB ID**":
        with col1:
            st.text_input(
                "🧬 Protein 1: PDB ID",
                placeholder="e.g., 5CXV",
                key="protein_1",
                help="Enter 4-character PDB identifier",
                max_chars=4,
            )
        with col2:
            st.text_input(
                "⛓️ Chain ID",
                placeholder="A",
                key="protein1_chain",
                help="Optional: defaults to 'A'",
                max_chars=1,
            )
        with col3:
            st.text_input(
                "🧬 Protein 2: PDB ID",
                placeholder="e.g., 5DSG",
                key="protein_2",
                help="Enter 4-character PDB identifier",
                max_chars=4,
            )
        with col4:
            st.text_input(
                "⛓️ Chain ID",
                placeholder="A",
                key="protein2_chain",
                help="Optional: defaults to 'A'",
                max_chars=1,
            )

    elif st.session_state.input_file_option == "**Upload Files**":
        with col1:
            st.file_uploader(
                "🧬 Protein 1 PDB File",
                key="protein_1_file",
                type=["pdb"],
                help="Upload PDB file (max 50MB)",
            )
        with col2:
            st.text_input(
                "⛓️ Chain ID",
                placeholder="A",
                key="file1_chain",
                help="Optional: defaults to 'A'",
                max_chars=1,
            )
        with col3:
            st.file_uploader(
                "🧬 Protein 2 PDB File",
                key="protein_2_file",
                type=["pdb"],
                help="Upload PDB file (max 50MB)",
            )
        with col4:
            st.text_input(
                "⛓️ Chain ID",
                placeholder="A",
                key="file2_chain",
                help="Optional: defaults to 'A'",
                max_chars=1,
            )

    # Submit button with improved styling
    st.sidebar.markdown("---")
    if st.sidebar.button(
        "Start Analysis",
        key="submit_button",
        on_click=main,
        type="primary",
        use_container_width=True,
    ):
        pass


@log_performance()
def view_pdb(
    selected_pdbs_ids: List[str], main_pdb_id: str, opacity: float, use_surface: bool
):
    """Visualize selected PDB structures with error handling"""
    if not selected_pdbs_ids:
        st.info("👆 Please select at least one structure to display.")
        return

    logger.debug(
        f"Visualizing PDB structures",
        extra={
            "main_pdb_id": main_pdb_id,
            "selected_structures": selected_pdbs_ids,
            "opacity": opacity,
            "use_surface": use_surface,
        },
    )

    try:
        # Initialize py3Dmol view
        view = py3Dmol.view()
        view.setBackgroundColor("white")

        # Color palette for different structures
        colors = [
            "red",
            "blue",
            "green",
            "orange",
            "purple",
            "cyan",
            "magenta",
            "yellow",
        ]

        # Add each selected protein model
        models_added = 0
        for i, pdb_id in enumerate(selected_pdbs_ids):
            try:
                pdb_data = protFuncs.read_saved_pdb_file_to_view(pdb_id, main_pdb_id)
                if pdb_data:
                    view.addModel(pdb_data, "pdb")
                    color = colors[models_added % len(colors)]
                    view.setStyle(
                        {"model": models_added}, {"cartoon": {"color": color}}
                    )

                    if use_surface:
                        view.addSurface(
                            py3Dmol.VDW,
                            {"color": color, "opacity": opacity},
                            {"model": models_added},
                        )

                    models_added += 1
                else:
                    st.warning(f"⚠️ Failed to load PDB data for {pdb_id}")
                    logger.warning(
                        f"Failed to load PDB data for visualization: {pdb_id}"
                    )
            except Exception as e:
                st.warning(f"⚠️ Error loading {pdb_id}: {str(e)}")
                logger.error(
                    f"Error loading PDB for visualization",
                    extra={"pdb_id": pdb_id, "error": str(e)},
                )

        if models_added > 0:
            # Zoom to fit all models
            view.zoomTo()
            # Render the visualization
            showmol(view, height=500, width=800)

            logger.info(f"Successfully visualized {models_added} structures")
        else:
            st.error("❌ No structures could be loaded for visualization")
            logger.error("No structures could be loaded for visualization")

    except Exception as e:
        st.error(f"❌ Visualization error: {str(e)}")
        logger.error(f"Visualization error: {str(e)}", exc_info=True)


def render_protein_view_ui(
    pdb_id: str,
    other_pdb_id: str,
    num_pockets: int,
    checkbox_states_key: str,
    opacity_key: str,
    use_surface_key: str,
):
    """Render protein view UI with improved controls"""
    st.subheader(f"{pdb_id} Structure Analysis")

    # Protein info
    col1, col2 = st.columns(2)
    with col1:
        st.metric("Protein ID", pdb_id)
    with col2:
        st.metric("Predicted Pockets", num_pockets)

    # Visualization controls
    with st.expander("🎛️ Visualization Controls", expanded=True):
        col1, col2 = st.columns(2)
        with col1:
            opacity = st.slider(
                "Surface Opacity",
                min_value=0.0,
                max_value=1.0,
                value=0.5,
                step=0.1,
                key=opacity_key,
                help="Adjust protein surface transparency",
            )
        with col2:
            use_surface = st.checkbox(
                "Show Surface",
                value=True,
                key=use_surface_key,
                help="Display protein surface representation",
            )

    # Structure selection
    st.write("**Select structures to display:**")

    # Initialize checkbox states if not exists
    if checkbox_states_key not in st.session_state:
        st.session_state[checkbox_states_key] = {}

    # Main protein checkbox
    st.session_state[checkbox_states_key][f"{pdb_id}"] = st.checkbox(
        f"🧬 {pdb_id} (Full Structure)",
        value=True,
        key=f"{pdb_id}_view",
        help=f"Show complete structure of {pdb_id}",
    )

    # Pocket selection checkboxes
    if num_pockets > 0:
        st.write("**Select binding pockets:**")
        cols = st.columns(num_pockets)  # Max 4 columns for better layout
        for indx in range(num_pockets):
            col_idx = indx % num_pockets
            with cols[col_idx]:
                pocket_key = f"{pdb_id}_pkt_{indx + 1}"
                st.session_state[checkbox_states_key][pocket_key] = st.checkbox(
                    f"{indx + 1}",
                    value=False,
                    key=pocket_key,
                    help=f"Show binding pocket {indx + 1}",
                )

    # 3D structure comparison checkbox
    comparison_key = f"{pdb_id}_flex_with_{other_pdb_id}"
    st.session_state[checkbox_states_key][comparison_key] = st.checkbox(
        f"Structural Alignment with {other_pdb_id}",
        value=True,
        key=comparison_key,
        help=f"Show structurally similar regions between {pdb_id} and {other_pdb_id}",
    )

    # Get selected PDB IDs from checkbox states
    selected_pdbs_ids = [
        key for key, value in st.session_state[checkbox_states_key].items() if value
    ]

    # Visualize selected structures
    view_pdb(selected_pdbs_ids, pdb_id, opacity, use_surface)


@log_performance()
def view_pocket_dx(
    selected_pkt_dx: str, main_pdb_id: str, pocket_list: dict, pocket_dx_list: dict
):
    """Display electrostatic potential visualization"""
    if not selected_pkt_dx:
        st.info("👆 Please select a pocket to display electrostatic potential.")
        return

    logger.debug(
        f"Viewing pocket electrostatic potential",
        extra={"main_pdb_id": main_pdb_id, "selected_pocket": selected_pkt_dx},
    )

    try:
        # Load PDB and DX files
        selected_pkt_file = pocket_list[selected_pkt_dx]
        selected_pkt_dx_file = pocket_dx_list[selected_pkt_dx]

        pdb_content = protFuncs.parse_dx_file_for_view(selected_pkt_file, main_pdb_id)
        dx_content = protFuncs.parse_dx_file_for_view(selected_pkt_dx_file, main_pdb_id)

        if not pdb_content or not dx_content:
            st.error("❌ Failed to load pocket or electrostatic data")
            logger.error("Failed to load pocket or electrostatic data")
            return

        # Initialize py3Dmol viewer
        view = py3Dmol.view()
        view.setBackgroundColor("white")

        # Add protein structure
        view.addModel(pdb_content, "pqr")
        view.addSurface(
            py3Dmol.VDW,
            {
                "voldata": dx_content,
                "volformat": "dx",
                "volscheme": {"gradient": "rwb", "min": -5, "max": 5},
                "opacity": 0.9,
            },
        )

        # Set visualization style
        view.setStyle({"cartoon": {"color": "spectrum"}})
        view.zoomTo()

        # Render visualization
        showmol(view, height=500, width=800)

        # Add legend
        st.markdown(
            """
        **Electrostatic Potential Legend:**
        - 🔴 **Red**: Negative potential (electron-rich regions)
        - ⚪ **White**: Neutral potential
        - 🔵 **Blue**: Positive potential (electron-poor regions)
        """
        )

        logger.info(
            f"Successfully displayed electrostatic potential for {selected_pkt_dx}"
        )

    except Exception as e:
        st.error(f"❌ Error displaying electrostatic potential: {str(e)}")
        logger.error(
            f"Error displaying electrostatic potential",
            extra={"pocket": selected_pkt_dx, "error": str(e)},
            exc_info=True,
        )


def render_pocket_dx_view_ui(main_pdb_id: str, num_pockets: int, pocket_dx_key: str):
    """Render electrostatic potential visualization UI"""
    st.subheader(f"{main_pdb_id} Electrostatic Potential")

    if num_pockets == 0:
        st.warning("⚠️ No pockets available for electrostatic analysis")
        return

    # Create pocket selection dictionaries
    pocket_dx_list = {}
    pocket_list = {}
    for indx in range(1, num_pockets + 1):
        pocket_name = f"Pocket {indx}"
        pocket_dx_list[pocket_name] = f"{main_pdb_id}_pkt_{indx}_processed_pot.dx"
        pocket_list[pocket_name] = f"{main_pdb_id}_pkt_{indx}_processed.pqr"

    # Pocket selection dropdown
    selected_pkt_dx = st.selectbox(
        "Select pocket for electrostatic analysis:",
        list(pocket_dx_list.keys()),
        index=0,
        key=pocket_dx_key,
        help="Choose a binding pocket to visualize its electrostatic potential",
    )

    # Display visualization
    view_pocket_dx(selected_pkt_dx, main_pdb_id, pocket_list, pocket_dx_list)


def render_results_ui():
    """Render comprehensive results UI"""
    st.title("🧬 Protein Analysis Results", anchor=False)

    # Display analysis ID if available
    if st.session_state.get("analysis_id"):
        st.caption(f"Analysis ID: {st.session_state.analysis_id}")

    # Results tabs
    tab1, tab2, tab3 = st.tabs(
        ["Structure Analysis", "Electrostatic Analysis", "Summary"]
    )

    with tab1:
        render_structure_analysis_tab()

    with tab2:
        render_electrostatic_analysis_tab()

    with tab3:
        render_summary_tab()


def render_structure_analysis_tab():
    """Render structure analysis tab"""
    st.header(
        f"Structural Analysis: {st.session_state.pdb_1_id} vs {st.session_state.pdb_2_id}"
    )

    # Structure visualizations
    col1, col2 = st.columns(2)

    with col1:
        render_protein_view_ui(
            st.session_state.pdb_1_id,
            st.session_state.pdb_2_id,
            st.session_state.pdb_1_num_pocket,
            "checkbox_states_1",
            "prot_1_opacity",
            "prot_1_use_surface",
        )

    with col2:
        render_protein_view_ui(
            st.session_state.pdb_2_id,
            st.session_state.pdb_1_id,
            st.session_state.pdb_2_num_pocket,
            "checkbox_states_2",
            "prot_2_opacity",
            "prot_2_use_surface",
        )

    # Overlap scores section
    st.markdown("---")
    st.subheader("Pocket-Structure Overlap Analysis")

    col1, col2 = st.columns(2)

    with col1:
        st.write(f"**{st.session_state.pdb_1_id} Overlap Scores**")
        st.caption(
            f"Overlap between {st.session_state.pdb_1_id} pockets and structurally similar regions"
        )

        with st.expander("View Detailed Scores", expanded=False):
            try:
                overlap_df = protFuncs.get_pocket_overlap_score_table(
                    st.session_state.pdb_1_id
                )
                st.dataframe(
                    overlap_df,
                    column_config={
                        "3D_Structurally_Similar_Part": st.column_config.TextColumn(
                            "Structural Comparison"
                        ),
                        f"{st.session_state.pdb_1_id}_Pocket_Number": st.column_config.NumberColumn(
                            f"{st.session_state.pdb_1_id} Pocket #"
                        ),
                        "Pocket_and_3D_Similar_Part_Overlap_Score": st.column_config.ProgressColumn(
                            "Overlap Score", min_value=0, max_value=1
                        ),
                    },
                    hide_index=True,
                    use_container_width=True,
                )
            except Exception as e:
                st.error(f"Error loading overlap scores: {str(e)}")
                logger.error(
                    f"Error loading overlap scores for {st.session_state.pdb_1_id}: {str(e)}"
                )

    with col2:
        st.write(f"**{st.session_state.pdb_2_id} Overlap Scores**")
        st.caption(
            f"Overlap between {st.session_state.pdb_2_id} pockets and structurally similar regions"
        )

        with st.expander("View Detailed Scores", expanded=False):
            try:
                overlap_df = protFuncs.get_pocket_overlap_score_table(
                    st.session_state.pdb_2_id
                )
                st.dataframe(
                    overlap_df,
                    column_config={
                        "3D_Structurally_Similar_Part": st.column_config.TextColumn(
                            "Structural Comparison"
                        ),
                        f"{st.session_state.pdb_2_id}_Pocket_Number": st.column_config.NumberColumn(
                            f"{st.session_state.pdb_2_id} Pocket #"
                        ),
                        "Pocket_and_3D_Similar_Part_Overlap_Score": st.column_config.ProgressColumn(
                            "Overlap Score", min_value=0, max_value=1
                        ),
                    },
                    hide_index=True,
                    use_container_width=True,
                )
            except Exception as e:
                st.error(f"Error loading overlap scores: {str(e)}")
                logger.error(
                    f"Error loading overlap scores for {st.session_state.pdb_2_id}: {str(e)}"
                )

    # Pocket comparison section
    st.markdown("---")
    st.subheader("Binding Pocket Comparison (APoc Analysis)")

    with st.expander("Detailed Pocket Comparison Results", expanded=False):
        try:
            comparison_df = protFuncs.get_combine_pockets_comparison_table()

            st.dataframe(
                comparison_df,
                column_config={
                    "pocket_1": st.column_config.TextColumn("Pocket 1"),
                    "pocket_2": st.column_config.TextColumn("Pocket 2"),
                    "PS_score": st.column_config.NumberColumn(
                        "PS Score", format="%.3f"
                    ),
                    "P_value": st.column_config.NumberColumn("P-value", format="%.2e"),
                    "RMSD": st.column_config.NumberColumn("RMSD (Å)", format="%.2f"),
                    "Number_aligned_residues": st.column_config.NumberColumn(
                        "Aligned Residues"
                    ),
                    "Seq_Identity": st.column_config.ProgressColumn(
                        "Sequence Identity", min_value=0, max_value=1
                    ),
                },
                hide_index=True,
                use_container_width=True,
            )

            # Add explanation
            st.info(
                """
            **APoc Analysis Explanation:**
            - **PS Score**: Pocket similarity score (higher = more similar)
            - **P-value**: Statistical significance of the similarity
            - **RMSD**: Root mean square deviation of aligned atoms
            - **Sequence Identity**: Fraction of identical residues in alignment
            """
            )

        except Exception as e:
            st.error(f"Error loading pocket comparison data: {str(e)}")
            logger.error(f"Error loading pocket comparison data: {str(e)}")


def render_electrostatic_analysis_tab():
    """Render electrostatic analysis tab"""
    st.header("Electrostatic Potential Analysis")

    # Heatmap and clustering analysis
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Similarity Heatmap")
        try:
            if not st.session_state.pocket_similarity_matrix.empty:
                heat_map = protFuncs.heatmap_pocket_corr(
                    st.session_state.pocket_similarity_matrix,
                    st.session_state.pdb_1_id,
                    st.session_state.pdb_2_id,
                )
                st.pyplot(heat_map, use_container_width=True)
            else:
                st.warning("⚠️ Similarity matrix not available")
        except Exception as e:
            st.error(f"Error generating heatmap: {str(e)}")
            logger.error(f"Error generating heatmap: {str(e)}")

    with col2:
        st.subheader("Cluster Analysis")
        try:
            if not st.session_state.pocket_similarity_matrix.empty:
                dendrogram = protFuncs.cluster_pockets(
                    st.session_state.pocket_similarity_matrix,
                    st.session_state.pdb_1_id,
                    st.session_state.pdb_2_id,
                )
                st.pyplot(dendrogram, use_container_width=True)
            else:
                st.warning("⚠️ Similarity matrix not available")
        except Exception as e:
            st.error(f"Error generating dendrogram: {str(e)}")
            logger.error(f"Error generating dendrogram: {str(e)}")

    # Individual pocket electrostatic potential visualizations
    st.markdown("---")
    st.subheader("Individual Pocket Electrostatic Potentials")

    col1, col2 = st.columns(2)

    with col1:
        render_pocket_dx_view_ui(
            st.session_state.pdb_1_id, st.session_state.pdb_1_num_pocket, "pocket_dx_1"
        )

    with col2:
        render_pocket_dx_view_ui(
            st.session_state.pdb_2_id, st.session_state.pdb_2_num_pocket, "pocket_dx_2"
        )


@log_performance()
def render_summary_tab():
    """Render analysis summary tab with integrated utilities"""
    st.header("📈 Analysis Summary")

    # Summary metrics
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric(
            "Protein 1",
            st.session_state.pdb_1_id,
            f"{st.session_state.pdb_1_num_pocket} pockets",
        )

    with col2:
        st.metric(
            "Protein 2",
            st.session_state.pdb_2_id,
            f"{st.session_state.pdb_2_num_pocket} pockets",
        )

    with col3:
        total_pockets = (
            st.session_state.pdb_1_num_pocket + st.session_state.pdb_2_num_pocket
        )
        st.metric("Total Pockets", total_pockets)

    with col4:
        total_comparisons = (
            st.session_state.pdb_1_num_pocket * st.session_state.pdb_2_num_pocket
        )
        st.metric("Pocket Comparisons", total_comparisons)

    # Enhanced download section
    st.markdown("---")
    st.subheader("💾 Download Results")

    col1, col2 = st.columns(2)

    with col1:
        # Individual file downloads
        download_options = {
            "Pocket Comparison Results": "tempDownloadDir/results/apoc/Combine_pocket_comp_results.csv",
            "Similarity Matrix": f"tempDownloadDir/results/apbs/{st.session_state.pdb_1_id}_{st.session_state.pdb_2_id}_pocket_corr_similarity_matrix.csv",
        }

        for name, path in download_options.items():
            try:
                if os.path.exists(path):
                    with open(path, "rb") as file:
                        st.download_button(
                            label=f"Download {name}",
                            data=file.read(),
                            file_name=os.path.basename(path),
                            mime="text/csv",
                            use_container_width=True,
                        )
                    logger.debug(f"Download button created for {name}")
            except Exception as e:
                st.warning(f"Could not prepare {name} for download: {str(e)}")
                logger.warning(f"Could not prepare download for {name}: {str(e)}")

    with col2:
        # Complete results package
        st.write("**📦 Complete Results Package**")

        if st.button("Create Results Package", use_container_width=True):
            try:
                package_data = protFuncs.create_results_package(
                    st.session_state.pdb_1_id,
                    st.session_state.pdb_2_id,
                    st.session_state.pdb_1_num_pocket,
                    st.session_state.pdb_2_num_pocket,
                )

                if package_data:
                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    package_filename = f"bioops_results_{st.session_state.pdb_1_id}_{st.session_state.pdb_2_id}_{timestamp}.zip"

                    st.download_button(
                        label="📥 Download Complete Package",
                        data=package_data,
                        file_name=package_filename,
                        mime="application/zip",
                        use_container_width=True,
                    )

                    logger.info(
                        f"Results package created",
                        extra={
                            "filename_name": package_filename,
                            "size_bytes": len(package_data),
                        },
                    )
            except Exception as e:
                st.error(f"Error creating results package: {str(e)}")
                logger.error(f"Error creating results package: {str(e)}", exc_info=True)


def render_welcome_screen():
    """Render welcome screen with instructions"""
    # Show example or logo if available
    config = get_app_config()
    if os.path.exists(config["logo_main"]):
        st.image(config["logo_main"], width=400)

    st.markdown(
        """
    # Welcome to BioOps Protein Analysis Tool
    
    **Comprehensive G-Protein Coupled Receptor (GPCR) Analysis Platform**
    
    ## 🚀 Quick Start Guide
    
    1. **Choose Input Method**: Select either PDB ID entry or file upload
    2. **Provide Structures**: Enter PDB IDs or upload your protein structure files
    3. **Specify Chains**: Optionally specify protein chains (defaults to 'A')
    4. **Start Analysis**: Click the analysis button and wait for results
    
    ## 🔬 Analysis Features
    
    - **Structural Alignment**: Compare 3D protein structures using FATCAT
    - **Pocket Prediction**: Identify binding pockets with P2Rank
    - **Pocket Comparison**: Compare pocket similarities using APoc
    - **Overlap Scoring**: Calculate overlap scores between pockets and structurally similar regions
    - **Electrostatic Analysis**: Calculate and visualize electrostatic potentials
    - **Interactive Visualization**: 3D molecular viewers with customizable displays
    
    ## 📋 Requirements
    
    - **PDB Files**: Standard Protein Data Bank format
    - **File Size**: Maximum 50MB per file
    - **PDB IDs**: 4-character identifiers (e.g., 5CXV, 5DSG)
    - **Chain IDs**: Single character identifiers
    
    ---
    
    **👈 Get started by selecting your input method in the sidebar!**
    """
    )


# APP ENTRY POINT
if __name__ == "__main__":
    try:
        # Initialize logging
        logger.info("=" * 50)
        logger.info("BioOps Application Starting")
        logger.info("=" * 50)

        # Check required tools
        tool_status = protFuncs.check_required_tools()
        missing_tools = [
            tool for tool, status in tool_status.items() if not status["available"]
        ]

        if missing_tools:
            logger.warning(f"Missing required tools: {', '.join(missing_tools)}")
            st.warning(
                f"⚠️ Missing required tools: {', '.join(missing_tools)}. Please run the install_tools script."
            )

        # Initialize session state
        initialize_session_state()

        # Ensure required directories exist
        ensure_directories()

        # Setup page
        setup_page()

        # Render input UI
        render_input_ui()

        # Show results or welcome screen
        if (
            st.session_state.processing_complete
            and not st.session_state.processing_error
        ):
            render_results_ui()
        elif st.session_state.processing_error:
            st.error(f"❌ Analysis failed: {st.session_state.error_message}")
            st.info("👈 Please check your inputs and try again.")
            render_welcome_screen()
        else:
            render_welcome_screen()

        # Footer
        st.markdown("---")
        st.markdown(
            """
        <div style='text-align: center; color: gray; font-size: 1.1em;'>
            <p>BioOps | 
            <a href="https://fatcat.godziklab.org/" target="_blank">FATCAT</a> | 
            <a href="https://doi.org/10.1093/bioinformatics/btt024" target="_blank">APoc</a> | 
            <a href="https://prankweb.cz/" target="_blank">P2rank</a> | 
            <a href"=https://server.poissonboltzmann.org/" target="_blank">PDB2PQR</a> | 
            <a href"=https://server.poissonboltzmann.org/" target="_blank">APBS</a> | 
            <a href="https://streamlit.io/" target="_blank">Streamlit</a> | 
            <a href="https://github.com/owusukd/bioops" target="_blank">GitHub</a> |
            Contact: <a href="mailto:owusukd@yahoo.com" ">owusukd@yahoo.com</a></p>
        </div>
        """,
            unsafe_allow_html=True,
        )

    except Exception as e:
        st.error(f"❌ Application error: {str(e)}")
        logger.error(f"Application error: {str(e)}", exc_info=True)
        st.info("Please refresh the page and try again.")
