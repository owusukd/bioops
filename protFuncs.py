# Import enhanced logging
from logging_config import (bioops_logger, log_analysis_step, log_context,
                            log_performance, log_tool_execution)

# Get logger instance
logger = bioops_logger.get_logger(__name__)

import csv
import io
import json
import os
import re
import shutil
import subprocess
import time
import zipfile
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Tuple, Union

import matplotlib as mpl
import matplotlib.pyplot as plt
# Scientific computing imports
import numpy as np
import pandas as pd
import seaborn as sns
# Streamlit specific imports
import streamlit as st
from gridData import Grid
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D
from scipy.cluster.hierarchy import dendrogram, linkage, optimal_leaf_ordering
from scipy.spatial.distance import squareform

# Constants
RESULTS_BASE_DIR = Path("tempDownloadDir/results")
FATCAT_DIR = RESULTS_BASE_DIR / "fatcat"
P2RANK_DIR = RESULTS_BASE_DIR / "p2rank"
APOC_DIR = RESULTS_BASE_DIR / "apoc"
APBS_DIR = RESULTS_BASE_DIR / "apbs"
OVERLAP_DIR = RESULTS_BASE_DIR / "overlap_score"

# Tool paths configuration
TOOL_PATHS = {
    "FATCAT": [
        "FATCAT",  # In PATH
        "./tools/FATCAT-dist/FATCATMain/FATCAT",
        "/tools/FATCAT-dist/FATCATMain/FATCAT",
        os.path.expanduser("~/tools/FATCAT-dist/FATCATMain/FATCAT"),
        "/usr/local/bin/FATCAT",
    ],
    "prank": [
        "prank",  # In PATH
        "./tools/p2rank_2.5/prank",
        "/tools/p2rank_2.5/prank",
        os.path.expanduser("~/tools/p2rank_2.5/prank"),
        "/usr/local/bin/prank",
    ],
    "apoc": [
        "apoc",  # In PATH
        "./tools/apoc/bin/apoc",
        "/tools/apoc/bin/apoc",
        os.path.expanduser("~/tools/apoc/bin/apoc"),
        "/usr/local/bin/apoc",
    ],
    "apbs": [
        "apbs",  # In PATH
        "./tools/APBS-3.4.1.Linux/bin/apbs",
        "/tools/APBS-3.4.1.Linux/bin/apbs",
        os.path.expanduser("~/tools/APBS-3.4.1.Linux/bin/apbs"),
        "/usr/local/bin/apbs",
    ],
    "pdb2pqr": [
        "pdb2pqr",  # Usually installed system-wide
        "/usr/bin/pdb2pqr",
        "/usr/local/bin/pdb2pqr",
    ],
}


def find_tool_path(tool_name: str) -> Optional[str]:
    """Find the path to a tool"""
    if tool_name not in TOOL_PATHS:
        return None

    for path in TOOL_PATHS[tool_name]:
        if os.path.exists(path) and os.access(path, os.X_OK):
            logger.debug(f"Found {tool_name} at: {path}")
            return path

    # Try which command as fallback
    try:
        result = subprocess.run(["which", tool_name], capture_output=True, text=True)
        if result.returncode == 0 and result.stdout.strip():
            path = result.stdout.strip()
            logger.debug(f"Found {tool_name} via which: {path}")
            return path
    except:
        logger.warning(f"{tool_name} not found in any expected location")
        return None


def check_required_tools() -> dict:
    """Check if all required tools are available"""
    tool_status = {}

    for tool in ["FATCAT", "prank", "apoc", "apbs", "pdb2pqr"]:
        path = find_tool_path(tool)
        tool_status[tool] = {"available": path is not None, "path": path}

        if path:
            logger.info(f"Tool {tool} found at: {path}")
        else:
            logger.warning(f"Tool {tool} not found!")

    return tool_status


class ProteinAnalysisError(Exception):
    """Base exception for protein analysis errors"""

    pass


class FileNotFoundError(ProteinAnalysisError):
    """Custom file not found error"""

    pass


class ExternalToolError(ProteinAnalysisError):
    """Error running external analysis tools"""

    pass


# Utility Functions
@log_performance()
def ensure_directory(path: Union[str, Path]) -> Path:
    """Ensure directory exists and return Path object"""
    path_obj = Path(path)
    path_obj.mkdir(parents=True, exist_ok=True)
    logger.debug(f"Ensured directory exists: {path_obj}")
    return path_obj


def run_command(
    command: str, cwd: Optional[Path] = None, timeout: int = 300
) -> Tuple[int, str, str]:
    """Run shell command with error handling and timeout"""
    start_time = time.time()
    logger.info(
        f"Running command: {command[:100]}...",
        extra={"full_command": command, "cwd": str(cwd)},
    )

    try:
        # Check if command exists and replace with full path
        cmd_parts = command.split()
        cmd_name = cmd_parts[0]

        # Check if this is one of our known tools
        tool_path = find_tool_path(cmd_name)
        if tool_path:
            cmd_parts[0] = tool_path
            command = " ".join(cmd_parts)
            logger.info(f"Using tool path: {tool_path}")
        elif cmd_name in TOOL_PATHS:
            # Tool is required but not found
            raise ExternalToolError(
                f"{cmd_name} not found. Please run install_tools script or check installation."
            )

        result = subprocess.run(
            command,
            shell=True,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=cwd,
        )

        duration = time.time() - start_time

        # Log tool execution
        tool_name = cmd_name.upper()
        log_tool_execution(
            tool_name=tool_name,
            command=command,
            return_code=result.returncode,
            stdout=result.stdout,
            stderr=result.stderr,
            duration=duration,
        )

        return result.returncode, result.stdout, result.stderr

    except subprocess.TimeoutExpired:
        duration = time.time() - start_time
        logger.error(f"Command timed out after {timeout} seconds: {command}")
        raise ExternalToolError(f"Command timed out after {timeout} seconds: {command}")
    except ExternalToolError:
        raise  # Re-raise tool not found errors
    except Exception as e:
        duration = time.time() - start_time
        logger.error(f"Failed to run command: {command}", extra={"error": str(e)})
        raise ExternalToolError(f"Failed to run command '{command}': {str(e)}")


def validate_file_exists(
    file_path: Union[str, Path], description: str = "File"
) -> Path:
    """Validate that a file exists and return Path object"""
    path_obj = Path(file_path)
    if not path_obj.exists():
        logger.error(f"{description} not found: {path_obj}")
        raise FileNotFoundError(f"{description} not found: {path_obj}")
    if path_obj.stat().st_size == 0:
        logger.error(f"{description} is empty: {path_obj}")
        raise FileNotFoundError(f"{description} is empty: {path_obj}")

    logger.debug(f"Validated file exists: {path_obj} ({path_obj.stat().st_size} bytes)")
    return path_obj


@st.cache_data(show_spinner=False)
def read_file_content(file_path: str) -> str:
    """Read file content with caching"""
    try:
        with open(file_path, "r") as file:
            content = file.read()
            logger.debug(f"Read file content: {file_path} ({len(content)} chars)")
            return content
    except Exception as e:
        logger.error(f"Could not read file {file_path}: {str(e)}")
        raise FileNotFoundError(f"Could not read file {file_path}: {str(e)}")


###### Core Analysis Functions


def makeDir(path: Union[str, Path]) -> Path:
    """Create directory if it doesn't exist (legacy compatibility)"""
    return ensure_directory(path)


@log_performance()
def read_saved_pdb_file_to_view(pdb_file: str, main_pdb_id: str) -> Optional[str]:
    """Read PDB file for visualization with improved path handling"""
    try:
        # Determine file path based on file type
        if "flex_with" in pdb_file:
            file_path = FATCAT_DIR / f"{pdb_file}.pdb"
        elif "pkt" in pdb_file:
            file_path = P2RANK_DIR / f"{main_pdb_id}_pockets" / f"{pdb_file}.pdb"
        else:
            file_path = RESULTS_BASE_DIR / f"{pdb_file}.pdb"

        logger.debug(f"Reading PDB file for visualization: {file_path}")

        # Validate and read file
        validate_file_exists(file_path, f"PDB file {pdb_file}")
        return read_file_content(str(file_path))

    except (FileNotFoundError, ProteinAnalysisError) as e:
        st.error(f"Error loading {pdb_file}: {str(e)}")
        logger.error(f"Error in read_saved_pdb_file_to_view: {str(e)}")
        return None
    except Exception as e:
        st.error(f"Unexpected error loading {pdb_file}: {str(e)}")
        logger.error(f"Unexpected error in read_saved_pdb_file_to_view: {str(e)}")
        return None


@log_analysis_step("FATCAT Structural Alignment")
def FATCAT(pdb_1_id: str, pdb_2_id: str) -> bool:
    """Run FATCAT structural alignment with improved error handling"""
    try:
        with log_context(tool="FATCAT", pdb_1=pdb_1_id, pdb_2=pdb_2_id):
            # Ensure output directory exists
            output_dir = ensure_directory(FATCAT_DIR)

            # Validate input files
            pdb_1_path = validate_file_exists(
                RESULTS_BASE_DIR / f"{pdb_1_id}.pdb", f"PDB file for {pdb_1_id}"
            )
            pdb_2_path = validate_file_exists(
                RESULTS_BASE_DIR / f"{pdb_2_id}.pdb", f"PDB file for {pdb_2_id}"
            )

            # Prepare command
            output_path = output_dir / f"{pdb_1_id}_{pdb_2_id}_flex"
            command = f"FATCAT -p1 {pdb_1_path} -p2 {pdb_2_path} -o {output_path} -m"

            # Run FATCAT
            logger.info(f"Running FATCAT alignment: {pdb_1_id} vs {pdb_2_id}")
            returncode, stdout, stderr = run_command(
                command, timeout=600
            )  # 10 minute timeout

            if returncode != 1:  # FATCAT returns 1 on success
                raise ExternalToolError(
                    f"FATCAT failed with return code {returncode}. stderr: {stderr}"
                )

            # Verify output files were created
            expected_files = [
                f"{pdb_1_id}_{pdb_2_id}_flex.aln",
                f"{pdb_1_id}_{pdb_2_id}_flex.pdb",
            ]
            for filename in expected_files:
                output_file = output_dir / filename
                if not output_file.exists():
                    logger.warning(
                        f"Expected FATCAT output file not found: {output_file}"
                    )

            logger.info("FATCAT completed successfully")
            return True

    except (FileNotFoundError, ExternalToolError) as e:
        logger.error(f"FATCAT error: {str(e)}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error in FATCAT: {str(e)}")
        raise ExternalToolError(f"Unexpected FATCAT error: {str(e)}")


def get_residue_id(pdb_string: str) -> List[int]:
    """Extract residue IDs from PDB string with improved parsing"""
    try:
        # More robust regex for different PDB formats
        regex = r"^ATOM\s+\d+\s+\w+\s+\w{3}\s+(\w)\s*(\d+)"
        matches = re.findall(regex, pdb_string, re.MULTILINE)

        # Extract unique residue numbers
        residue_ids = []
        seen = set()

        for chain, res_num in matches:
            res_id = int(res_num)
            if res_id not in seen:
                residue_ids.append(res_id)
                seen.add(res_id)

        logger.debug(f"Extracted {len(residue_ids)} unique residue IDs")
        return sorted(residue_ids)

    except Exception as e:
        logger.error(f"Error parsing residue IDs: {str(e)}")
        raise ProteinAnalysisError(f"Failed to parse residue IDs: {str(e)}")


@log_analysis_step("Extract Superimposed 3D Parts")
def getSuperimposed3Dpart(pdb_1_id: str, pdb_2_id: str) -> bool:
    """Extract superimposed 3D parts with improved error handling"""
    try:
        with log_context(tool="getSuperimposed3Dpart", pdb_1=pdb_1_id, pdb_2=pdb_2_id):
            # Validate input files
            aln_file = validate_file_exists(
                FATCAT_DIR / f"{pdb_1_id}_{pdb_2_id}_flex.aln", "FATCAT alignment file"
            )

            pdb_1_file = validate_file_exists(
                RESULTS_BASE_DIR / f"{pdb_1_id}.pdb", f"PDB file for {pdb_1_id}"
            )

            pdb_2_file = validate_file_exists(
                RESULTS_BASE_DIR / f"{pdb_2_id}.pdb", f"PDB file for {pdb_2_id}"
            )

            logger.info(
                f"Extracting superimposed 3D parts for {pdb_1_id} and {pdb_2_id}"
            )

            # Read files
            with open(aln_file, "r") as f:
                aln_lines = f.readlines()

            with open(pdb_1_file, "r") as f:
                pdb_1_content = f.read()

            with open(pdb_2_file, "r") as f:
                pdb_2_content = f.read()

            # Get residue IDs for both proteins
            position_prot_1_uniq = get_residue_id(pdb_1_content)
            position_prot_2_uniq = get_residue_id(pdb_2_content)

            # Parse alignment and extract aligned regions
            position_chain_1_flex = []
            position_chain_2_flex = []

            # Process alignment file to extract flexible regions
            temp_flex_1 = position_prot_1_uniq.copy()
            temp_flex_2 = position_prot_2_uniq.copy()

            for i, line in enumerate(aln_lines):
                line = line.strip()
                if line.startswith("Chain 1"):
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            start_pos = int(parts[-2])
                            ind_1 = temp_flex_1.index(start_pos)

                            # Handle gaps in alignment
                            alignment_seq = parts[-1].replace("\n", "")
                            gap_indices = [
                                j + ind_1
                                for j, char in enumerate(alignment_seq)
                                if char == "-"
                            ]

                            for gap_idx in gap_indices:
                                if gap_idx < len(temp_flex_1):
                                    temp_flex_1.insert(gap_idx, "-")

                            # Extract aligned positions
                            if i + 1 < len(aln_lines):
                                similarity_line = aln_lines[i + 1].replace("\n", "")
                                aligned_positions = [
                                    j - 14 + ind_1
                                    for j, char in enumerate(similarity_line)
                                    if char != " " and j - 14 + ind_1 < len(temp_flex_1)
                                ]

                                for pos in aligned_positions:
                                    if 0 <= pos < len(temp_flex_1):
                                        position_chain_1_flex.append(temp_flex_1[pos])

                        except (ValueError, IndexError) as e:
                            logger.warning(
                                f"Error parsing alignment line: {line}. Error: {str(e)}"
                            )
                            continue

                elif line.startswith("Chain 2"):
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            start_pos = int(parts[-2])
                            ind_2 = temp_flex_2.index(start_pos)

                            # Handle gaps in alignment
                            alignment_seq = parts[-1].replace("\n", "")
                            gap_indices = [
                                j + ind_2
                                for j, char in enumerate(alignment_seq)
                                if char == "-"
                            ]

                            for gap_idx in gap_indices:
                                if gap_idx < len(temp_flex_2):
                                    temp_flex_2.insert(gap_idx, "-")

                            # Extract aligned positions
                            if i - 1 >= 0:
                                similarity_line = aln_lines[i - 1].replace("\n", "")
                                aligned_positions = [
                                    j - 14 + ind_2
                                    for j, char in enumerate(similarity_line)
                                    if char != " " and j - 14 + ind_2 < len(temp_flex_2)
                                ]

                                for pos in aligned_positions:
                                    if 0 <= pos < len(temp_flex_2):
                                        position_chain_2_flex.append(temp_flex_2[pos])

                        except (ValueError, IndexError) as e:
                            logger.warning(
                                f"Error parsing alignment line: {line}. Error: {str(e)}"
                            )
                            continue

            # Write output PDB files
            output_1 = FATCAT_DIR / f"{pdb_1_id}_flex_with_{pdb_2_id}.pdb"
            output_2 = FATCAT_DIR / f"{pdb_2_id}_flex_with_{pdb_1_id}.pdb"

            # Re-read PDB files as lines for processing
            with open(pdb_1_file, "r") as f:
                pdb_1_lines = f.readlines()
            with open(pdb_2_file, "r") as f:
                pdb_2_lines = f.readlines()

            write_flexible_pdb(output_1, position_chain_1_flex, pdb_1_lines)
            write_flexible_pdb(output_2, position_chain_2_flex, pdb_2_lines)

            logger.info(
                f"Superimposed 3D parts extracted successfully. Flex residues: {len(position_chain_1_flex)}, {len(position_chain_2_flex)}"
            )
            return True

    except (FileNotFoundError, ProteinAnalysisError) as e:
        logger.error(f"Error in getSuperimposed3Dpart: {str(e)}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error in getSuperimposed3Dpart: {str(e)}")
        raise ProteinAnalysisError(f"Failed to extract superimposed parts: {str(e)}")


def write_flexible_pdb(
    output_path: Path, positions: List[Union[int, str]], pdb_lines: List[str]
):
    """Write PDB file containing only flexible/aligned regions"""
    try:
        with open(output_path, "w") as output_file:
            atoms_written = 0
            for position in positions:
                if position != "-":
                    for line in pdb_lines[:-1]:  # Exclude last line (END)
                        if line.startswith("ATOM"):
                            parts = line.split()
                            if len(parts) >= 6:
                                try:
                                    # Handle different PDB formats
                                    if (
                                        len(parts[4]) == 1
                                    ):  # Chain and residue number separate
                                        res_num = int(parts[5])
                                    else:  # Chain and residue number combined
                                        res_num = int(re.sub(r"\D", "", parts[4]))

                                    if res_num == position:
                                        output_file.write(line)
                                        atoms_written += 1
                                except (ValueError, IndexError):
                                    continue

            output_file.write("TER\n")  # End of PDB file

        logger.debug(f"Written flexible PDB: {output_path} ({atoms_written} atoms)")

    except Exception as e:
        logger.error(f"Failed to write flexible PDB file {output_path}: {str(e)}")
        raise ProteinAnalysisError(
            f"Failed to write flexible PDB file {output_path}: {str(e)}"
        )


@log_analysis_step("P2Rank Pocket Prediction")
def runP2rank(pdb_1_id: str, pdb_2_id: str) -> bool:
    """Run P2Rank pocket prediction with improved error handling"""
    try:
        with log_context(tool="P2Rank", pdb_1=pdb_1_id, pdb_2=pdb_2_id):
            # Ensure output directory exists
            output_dir = ensure_directory(P2RANK_DIR)

            # Validate input files
            pdb_1_path = validate_file_exists(
                RESULTS_BASE_DIR / f"{pdb_1_id}.pdb", f"PDB file for {pdb_1_id}"
            )
            pdb_2_path = validate_file_exists(
                RESULTS_BASE_DIR / f"{pdb_2_id}.pdb", f"PDB file for {pdb_2_id}"
            )

            # Run P2Rank for both proteins
            for pdb_id, pdb_path in [(pdb_1_id, pdb_1_path), (pdb_2_id, pdb_2_path)]:
                pocket_output = output_dir / f"{pdb_id}_pockets"
                command = f"prank predict -o {pocket_output} -f {pdb_path} -threads 8 -vis_copy_proteins 0"

                logger.info(f"Running P2Rank for {pdb_id}")
                returncode, stdout, stderr = run_command(
                    command, timeout=900
                )  # 15 minute timeout

                if returncode != 0:
                    raise ExternalToolError(
                        f"P2Rank failed for {pdb_id} with return code {returncode}. stderr: {stderr}"
                    )

                # Verify output was created
                expected_csv = pocket_output / f"{pdb_id}.pdb_predictions.csv"
                if not expected_csv.exists():
                    raise ExternalToolError(
                        f"P2Rank did not generate expected output for {pdb_id}: {expected_csv}"
                    )

                # Log pocket count
                with open(expected_csv, "r") as f:
                    pocket_count = len(f.readlines()) - 1  # Subtract header
                    logger.info(f"P2Rank found {pocket_count} pockets for {pdb_id}")

            logger.info("P2Rank completed successfully for both proteins")
            return True

    except (FileNotFoundError, ExternalToolError) as e:
        logger.error(f"P2Rank error: {str(e)}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error in P2Rank: {str(e)}")
        raise ExternalToolError(f"Unexpected P2Rank error: {str(e)}")


@log_performance()
def write_pocket_pdb(
    pdb_id: str, pocket_predictions: List[str], num_pockets: int
) -> bool:
    """Write individual pocket PDB files with improved parsing"""
    try:
        # Validate input PDB file
        pdb_file_path = validate_file_exists(
            RESULTS_BASE_DIR / f"{pdb_id}.pdb", f"PDB file for {pdb_id}"
        )

        with open(pdb_file_path, "r") as f:
            pdb_lines = f.readlines()

        logger.info(f"Writing {num_pockets} pocket PDB files for {pdb_id}")

        # Process each pocket
        for pocket_idx in range(num_pockets):
            if pocket_idx >= len(pocket_predictions):
                logger.warning(
                    f"Not enough pocket predictions for {pdb_id}, pocket {pocket_idx + 1}"
                )
                continue

            try:
                # Parse pocket information from CSV line
                pocket_info = pocket_predictions[pocket_idx].strip().split(",")
                if len(pocket_info) < 10:
                    logger.warning(
                        f"Invalid pocket prediction format for {pdb_id}, pocket {pocket_idx + 1}"
                    )
                    continue

                # Extract residue information (typically in column 9)
                residue_info = pocket_info[9].strip()

                # Parse chain and residue numbers
                if not residue_info:
                    logger.warning(
                        f"No residue information for {pdb_id}, pocket {pocket_idx + 1}"
                    )
                    continue

                # Extract chain identifier and residue numbers
                chain_match = re.match(r"^([A-Z]_)", residue_info)
                if not chain_match:
                    logger.warning(
                        f"Could not parse chain from residue info: {residue_info}"
                    )
                    continue

                chain_prefix = chain_match.group(1)
                residue_numbers = []

                # Split by chain prefix and extract numbers
                parts = residue_info.split(chain_prefix)[1:]
                for part in parts:
                    try:
                        res_num = int(re.search(r"\d+", part).group())
                        residue_numbers.append(res_num)
                    except (AttributeError, ValueError):
                        continue

                residue_numbers.sort()

                if not residue_numbers:
                    logger.warning(
                        f"No valid residue numbers found for {pdb_id}, pocket {pocket_idx + 1}"
                    )
                    continue

                # Write pocket PDB file
                pocket_output_dir = ensure_directory(P2RANK_DIR / f"{pdb_id}_pockets")
                pocket_file_path = (
                    pocket_output_dir / f"{pdb_id}_pkt_{pocket_idx + 1}.pdb"
                )

                with open(pocket_file_path, "w") as pocket_file:
                    # Write pocket header (required for APoc)
                    pocket_file.write(
                        f"PKT        11    101    {pdb_id}_pkt_{pocket_idx + 1}\n"
                    )

                    atom_count = 0
                    # Write atoms for this pocket
                    for res_num in residue_numbers:
                        for line in pdb_lines[:-1]:  # Exclude END line
                            if line.startswith("ATOM"):
                                parts = line.split()
                                if len(parts) >= 6:
                                    try:
                                        # Handle different PDB formats
                                        if len(parts[4]) == 1:
                                            line_res_num = int(parts[5])
                                        else:
                                            line_res_num = int(
                                                re.sub(r"\D", "", parts[4])
                                            )

                                        if line_res_num == res_num:
                                            pocket_file.write(line)
                                            atom_count += 1
                                    except (ValueError, IndexError):
                                        continue

                    pocket_file.write("TER\n")

                logger.debug(
                    f"Created pocket file: {pocket_file_path} ({atom_count} atoms, {len(residue_numbers)} residues)"
                )

            except Exception as e:
                logger.error(
                    f"Error processing pocket {pocket_idx + 1} for {pdb_id}: {str(e)}"
                )
                continue

        return True

    except Exception as e:
        logger.error(f"Error in write_pocket_pdb for {pdb_id}: {str(e)}")
        raise ProteinAnalysisError(
            f"Failed to write pocket PDB files for {pdb_id}: {str(e)}"
        )


@log_performance()
def combine_pdb_file_pocket_files(pdb_id: str, num_pockets: int) -> bool:
    """Combine main PDB file with pocket files for APoc analysis"""
    try:
        # Ensure APoc directory exists
        apoc_dir = ensure_directory(APOC_DIR)

        # Copy main PDB file to APoc directory
        main_pdb_src = validate_file_exists(
            RESULTS_BASE_DIR / f"{pdb_id}.pdb", f"Main PDB file for {pdb_id}"
        )
        main_pdb_dst = apoc_dir / f"{pdb_id}.pdb"

        shutil.copy2(main_pdb_src, main_pdb_dst)

        # Append pocket files
        pockets_appended = 0
        with open(main_pdb_dst, "a") as combined_file:
            for pocket_idx in range(1, num_pockets + 1):
                pocket_file_path = (
                    P2RANK_DIR / f"{pdb_id}_pockets" / f"{pdb_id}_pkt_{pocket_idx}.pdb"
                )

                if pocket_file_path.exists():
                    with open(pocket_file_path, "r") as pocket_file:
                        combined_file.write("\n")
                        combined_file.write(pocket_file.read())
                        pockets_appended += 1
                else:
                    logger.warning(f"Pocket file not found: {pocket_file_path}")

        logger.info(
            f"Combined PDB file created for {pdb_id} with {pockets_appended}/{num_pockets} pockets"
        )
        return True

    except Exception as e:
        logger.error(f"Error combining PDB files for {pdb_id}: {str(e)}")
        raise ProteinAnalysisError(
            f"Failed to combine PDB files for {pdb_id}: {str(e)}"
        )


@log_analysis_step("Extract Pocket PDB Files")
def extractPocket_pdb(pdb_1_id: str, pdb_2_id: str) -> Tuple[int, int]:
    """Extract pocket PDB files from P2Rank results with improved error handling"""
    try:
        with log_context(tool="extractPocket", pdb_1=pdb_1_id, pdb_2=pdb_2_id):
            # Process both proteins
            pocket_counts = []

            for pdb_id in [pdb_1_id, pdb_2_id]:
                # Validate P2Rank results file
                predictions_file = validate_file_exists(
                    P2RANK_DIR / f"{pdb_id}_pockets" / f"{pdb_id}.pdb_predictions.csv",
                    f"P2Rank predictions for {pdb_id}",
                )

                # Read predictions
                with open(predictions_file, "r") as f:
                    lines = f.readlines()

                # Skip header and get pocket data
                pocket_data = [line for line in lines[1:] if line.strip()]
                num_pockets = len(pocket_data)

                if num_pockets == 0:
                    logger.warning(f"No pockets found for {pdb_id}")
                    pocket_counts.append(0)
                    continue

                # Write individual pocket PDB files
                write_pocket_pdb(pdb_id, pocket_data, num_pockets)

                # Combine files for APoc
                combine_pdb_file_pocket_files(pdb_id, num_pockets)

                pocket_counts.append(num_pockets)
                logger.info(f"Extracted {num_pockets} pockets for {pdb_id}")

            return tuple(pocket_counts)

    except Exception as e:
        logger.error(f"Error in extractPocket_pdb: {str(e)}")
        raise ProteinAnalysisError(f"Failed to extract pocket PDB files: {str(e)}")


@log_analysis_step("APoc Pocket Comparison")
def runApoc(pdb_1_id: str, pdb_2_id: str) -> bool:
    """Run APoc binding pocket comparison with improved error handling"""
    try:
        with log_context(tool="APoc", pdb_1=pdb_1_id, pdb_2=pdb_2_id):
            # Ensure APoc directory exists
            apoc_dir = ensure_directory(APOC_DIR)

            # Validate input files
            pdb_1_path = validate_file_exists(
                apoc_dir / f"{pdb_1_id}.pdb", f"Combined PDB file for {pdb_1_id}"
            )
            pdb_2_path = validate_file_exists(
                apoc_dir / f"{pdb_2_id}.pdb", f"Combined PDB file for {pdb_2_id}"
            )

            # Prepare command
            output_path = apoc_dir / f"{pdb_1_id}_{pdb_2_id}_pocket_compare_results.txt"
            command = f"apoc {pdb_1_path} {pdb_2_path} -pvol 50 -plen 5 > {output_path}"

            # Run APoc
            logger.info(f"Running APoc comparison: {pdb_1_id} vs {pdb_2_id}")
            returncode, stdout, stderr = run_command(
                command, timeout=600
            )  # 10 minute timeout

            if returncode != 0:
                logger.warning(
                    f"APoc returned non-zero exit code {returncode}, but this may be normal"
                )

            # Verify output file was created and has content
            if not output_path.exists() or output_path.stat().st_size == 0:
                raise ExternalToolError(
                    f"APoc did not generate valid output: {output_path}"
                )

            logger.info("APoc completed successfully")
            return True

    except (FileNotFoundError, ExternalToolError) as e:
        logger.error(f"APoc error: {str(e)}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error in APoc: {str(e)}")
        raise ExternalToolError(f"Unexpected APoc error: {str(e)}")


@log_analysis_step("Parse Pocket Comparison Results")
def parse_pocket_comparison_result(pdb_1_id: str, pdb_2_id: str) -> bool:
    """Parse APoc results with improved error handling and validation"""
    try:
        with log_context(tool="parse_APoc", pdb_1=pdb_1_id, pdb_2=pdb_2_id):
            # Validate input file
            result_file = validate_file_exists(
                APOC_DIR / f"{pdb_1_id}_{pdb_2_id}_pocket_compare_results.txt",
                "APoc results file",
            )

            # Read and clean lines
            with open(result_file, "r") as f:
                lines = [line.strip() for line in f if not line.isspace()]

            if not lines:
                raise ProteinAnalysisError("APoc results file is empty")

            logger.info(f"Parsing APoc results from {len(lines)} lines")

            # Initialize results structure
            columns = [
                "pocket_1",
                "pocket_2",
                "PS_score",
                "P_value",
                "RMSD",
                "Number_aligned_residues",
                "Seq_Identity",
            ]
            results = {col: [] for col in columns}

            # Parse results
            i = 0
            parsed_count = 0
            while i < len(lines):
                line = lines[i]

                if line.startswith(">>>>>>>>>>>>>>>>>>>>>>>>>") and "Pocket" in line:
                    try:
                        # Check if we have enough lines for a complete result
                        if i + 4 >= len(lines):
                            break

                        # Skip lines that indicate no results
                        if (
                            "The number of values to be sorted is not positive."
                            in lines[i + 1]
                        ):
                            i += 2
                            continue

                        # Extract pocket information
                        pocket_line_1 = lines[i + 1]
                        pocket_line_2 = lines[i + 2] if i + 2 < len(lines) else ""

                        # Parse pocket IDs
                        pocket_1_match = re.search(r"Pocket:\s*(\w+)", pocket_line_1)
                        pocket_2_match = re.search(r"Pocket:\s*(\w+)", pocket_line_2)

                        if not (pocket_1_match and pocket_2_match):
                            i += 1
                            continue

                        results["pocket_1"].append(pocket_1_match.group(1))
                        results["pocket_2"].append(pocket_2_match.group(1))

                        # Parse statistics
                        stats_line = lines[i + 3] if i + 3 < len(lines) else ""
                        residues_line = lines[i + 4] if i + 4 < len(lines) else ""
                        rmsd_line = lines[i + 5] if i + 5 < len(lines) else ""

                        # Extract PS score and P-value
                        ps_match = re.search(r"PS-score\s*=\s*([\d.-]+)", stats_line)
                        pval_match = re.search(r"P-value\s*=\s*([\d.e-]+)", stats_line)

                        results["PS_score"].append(
                            ps_match.group(1) if ps_match else "0.0"
                        )
                        results["P_value"].append(
                            pval_match.group(1) if pval_match else "1.0"
                        )

                        # Extract number of aligned residues
                        residues_match = re.search(r"=\s*(\d+)", residues_line)
                        results["Number_aligned_residues"].append(
                            residues_match.group(1) if residues_match else "0"
                        )

                        # Extract RMSD and sequence identity
                        rmsd_match = re.search(r"RMSD\s*=\s*([\d.-]+)", rmsd_line)
                        seq_id_match = re.search(
                            r"Seq identity\s*=\s*([\d.-]+)", rmsd_line
                        )

                        results["RMSD"].append(
                            rmsd_match.group(1) if rmsd_match else "0.0"
                        )
                        results["Seq_Identity"].append(
                            seq_id_match.group(1) if seq_id_match else "0.0"
                        )

                        parsed_count += 1
                        i += 6

                    except (IndexError, AttributeError) as e:
                        logger.warning(
                            f"Error parsing APoc result at line {i}: {str(e)}"
                        )
                        i += 1
                        continue
                else:
                    i += 1

            # Write results to CSV
            output_file = APOC_DIR / "Combine_pocket_comp_results.csv"

            with open(output_file, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(columns)

                # Write data rows
                for i in range(len(results["pocket_1"])):
                    row = [
                        results[col][i] if i < len(results[col]) else ""
                        for col in columns
                    ]
                    writer.writerow(row)

            logger.info(
                f"Parsed {parsed_count} pocket comparisons and saved to {output_file}"
            )
            return True

    except Exception as e:
        logger.error(f"Error parsing APoc results: {str(e)}")
        raise ProteinAnalysisError(f"Failed to parse APoc results: {str(e)}")


def intersection_len(list1: List, list2: List) -> int:
    """Calculate intersection length between two lists"""
    intersection = len(set(list1) & set(list2))
    logger.debug(
        f"Intersection calculation: {len(list1)} ∩ {len(list2)} = {intersection}"
    )
    return intersection


def score_overlap(position_pocket: List[int], position_3d_similar: List[int]) -> float:
    """Calculate overlap score between pocket and 3D similar region"""
    try:
        if not position_pocket:
            return 0.0

        intersection_count = intersection_len(position_pocket, position_3d_similar)
        score = intersection_count / len(position_pocket)

        logger.debug(
            f"Overlap score: {intersection_count}/{len(position_pocket)} = {score:.3f}"
        )
        return score

    except Exception as e:
        logger.warning(f"Error calculating overlap score: {str(e)}")
        return 0.0


@st.cache_data(show_spinner=False)
def get_pocket_overlap_score_table(pdb_id: str) -> pd.DataFrame:
    """Get pocket overlap score table with caching"""
    try:
        score_file = OVERLAP_DIR / f"{pdb_id}_overlap_scores.csv"
        validate_file_exists(score_file, f"Overlap scores for {pdb_id}")
        df = pd.read_csv(score_file)
        logger.debug(f"Loaded overlap scores for {pdb_id}: {df.shape}")
        return df
    except Exception as e:
        logger.error(f"Error loading overlap scores for {pdb_id}: {str(e)}")
        return pd.DataFrame()


@st.cache_data(show_spinner=False)
def get_combine_pockets_comparison_table() -> pd.DataFrame:
    """Get combined pocket comparison table with caching"""
    try:
        comparison_file = APOC_DIR / "Combine_pocket_comp_results.csv"
        validate_file_exists(comparison_file, "Combined pocket comparison results")
        df = pd.read_csv(comparison_file)
        logger.debug(f"Loaded pocket comparison table: {df.shape}")
        return df
    except Exception as e:
        logger.error(f"Error loading pocket comparison table: {str(e)}")
        return pd.DataFrame()


@log_analysis_step("Calculate Pocket Overlap Scores")
def pocket_3D_similar_overlap_score(
    pdb_1_id: str, pdb_2_id: str, num_pockets: int
) -> bool:
    """Calculate overlap scores between pockets and 3D similar regions"""
    try:
        with log_context(
            tool="overlap_score",
            pdb_1=pdb_1_id,
            pdb_2=pdb_2_id,
            num_pockets=num_pockets,
        ):
            # Ensure output directory
            overlap_dir = ensure_directory(OVERLAP_DIR)

            # Load 3D similar region
            flex_file = validate_file_exists(
                FATCAT_DIR / f"{pdb_1_id}_flex_with_{pdb_2_id}.pdb",
                f"3D similar region for {pdb_1_id}",
            )

            with open(flex_file, "r") as f:
                flex_content = f.read()

            similar_residues = get_residue_id(flex_content)
            logger.info(f"Found {len(similar_residues)} residues in 3D similar region")

            # Initialize results
            columns = [
                "3D_Structurally_Similar_Part",
                f"{pdb_1_id}_Pocket_Number",
                "Pocket_and_3D_Similar_Part_Overlap_Score",
            ]
            results = {col: [] for col in columns}

            # Calculate overlap for each pocket
            for pocket_num in range(1, num_pockets + 1):
                try:
                    pocket_file = (
                        P2RANK_DIR
                        / f"{pdb_1_id}_pockets"
                        / f"{pdb_1_id}_pkt_{pocket_num}.pdb"
                    )

                    if not pocket_file.exists():
                        logger.warning(f"Pocket file not found: {pocket_file}")
                        continue

                    with open(pocket_file, "r") as f:
                        pocket_content = f.read()

                    pocket_residues = get_residue_id(pocket_content)
                    overlap_score = score_overlap(pocket_residues, similar_residues)

                    results["3D_Structurally_Similar_Part"].append(
                        f"{pdb_1_id}_compared_with_{pdb_2_id}"
                    )
                    results[f"{pdb_1_id}_Pocket_Number"].append(pocket_num)
                    results["Pocket_and_3D_Similar_Part_Overlap_Score"].append(
                        overlap_score
                    )

                    logger.debug(
                        f"Pocket {pocket_num}: {len(pocket_residues)} residues, overlap score: {overlap_score:.3f}"
                    )

                except Exception as e:
                    logger.warning(
                        f"Error processing pocket {pocket_num} for {pdb_1_id}: {str(e)}"
                    )
                    continue

            # Write results
            output_file = overlap_dir / f"{pdb_1_id}_overlap_scores.csv"

            with open(output_file, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(columns)

                for i in range(len(results[columns[0]])):
                    row = [results[col][i] for col in columns]
                    writer.writerow(row)

            logger.info(
                f"Calculated overlap scores for {len(results[columns[0]])} pockets of {pdb_1_id}"
            )
            return True

    except Exception as e:
        logger.error(f"Error calculating overlap scores for {pdb_1_id}: {str(e)}")
        raise ProteinAnalysisError(f"Failed to calculate overlap scores: {str(e)}")


@log_analysis_step("PDB2PQR Processing")
def run_pockets_pdb2pqr(pdb_id: str, num_pockets: int) -> bool:
    """Run PDB2PQR for pocket electrostatic analysis"""
    try:
        with log_context(tool="PDB2PQR", pdb_id=pdb_id, num_pockets=num_pockets):
            # Ensure output directory
            output_dir = ensure_directory(APBS_DIR / pdb_id)

            logger.info(f"Running PDB2PQR for {num_pockets} pockets of {pdb_id}")
            processed_count = 0

            # Process each pocket
            for pocket_num in range(1, num_pockets + 1):
                pocket_file = validate_file_exists(
                    P2RANK_DIR / f"{pdb_id}_pockets" / f"{pdb_id}_pkt_{pocket_num}.pdb",
                    f"Pocket {pocket_num} PDB file for {pdb_id}",
                )

                output_file = output_dir / f"{pdb_id}_pkt_{pocket_num}_processed.pqr"
                command = (
                    f"pdb2pqr --ff=AMBER --with-ph=7.0 {pocket_file} {output_file}"
                )

                logger.debug(f"Running PDB2PQR for {pdb_id} pocket {pocket_num}")
                returncode, stdout, stderr = run_command(command, timeout=300)

                if returncode != 0:
                    logger.warning(
                        f"PDB2PQR failed for {pdb_id} pocket {pocket_num}: {stderr}"
                    )
                    continue

                # Verify output
                if not output_file.exists() or output_file.stat().st_size == 0:
                    logger.warning(
                        f"PDB2PQR did not generate valid output for {pdb_id} pocket {pocket_num}"
                    )
                else:
                    processed_count += 1

            logger.info(
                f"PDB2PQR completed for {pdb_id}: {processed_count}/{num_pockets} pockets processed"
            )
            return True

    except Exception as e:
        logger.error(f"Error in PDB2PQR for {pdb_id}: {str(e)}")
        raise ExternalToolError(f"PDB2PQR failed for {pdb_id}: {str(e)}")


@log_analysis_step("APBS Electrostatic Calculation")
def run_pockets_apbs(pdb_id: str, num_pockets: int) -> bool:
    """Run APBS for electrostatic potential calculation"""
    try:
        with log_context(tool="APBS", pdb_id=pdb_id, num_pockets=num_pockets):
            # Ensure output directory
            output_dir = ensure_directory(APBS_DIR / pdb_id)

            logger.info(f"Running APBS for {num_pockets} pockets of {pdb_id}")
            processed_count = 0

            # Process each pocket
            for pocket_num in range(1, num_pockets + 1):
                pqr_file = output_dir / f"{pdb_id}_pkt_{pocket_num}_processed.pqr"
                dx_output_file = output_dir / f"{pdb_id}_pkt_{pocket_num}_processed_pot"

                if not pqr_file.exists():
                    logger.warning(
                        f"PQR file not found for {pdb_id} pocket {pocket_num}"
                    )
                    continue

                # Create APBS input file
                input_file_content = f"""
read
    mol pqr {pqr_file}
end
elec name solv 
    mg-auto
    dime 97 97 97
    cglen 80.0 80.0 80.0
    fglen 60.0 60.0 60.0
    cgcent mol 1
    fgcent mol 1
    mol 1
    lpbe
    bcfl sdh
    srfm smol
    chgm spl2
    ion charge +1 conc 0.150 radius 2.0
    ion charge -1 conc 0.150 radius 2.0
    pdie 2.0
    sdie 78.54
    srad 1.4
    swin 0.3
    temp 298.15
    sdens 10.0
    calcenergy total
    write pot dx {dx_output_file}
end
"""

                input_file = output_dir / f"{pdb_id}_pkt_{pocket_num}_apbs.in"
                with open(input_file, "w") as f:
                    f.write(input_file_content)

                # Run APBS
                command = f"apbs {input_file}"
                logger.debug(f"Running APBS for {pdb_id} pocket {pocket_num}")

                returncode, stdout, stderr = run_command(command, cwd="./", timeout=600)

                if returncode != 0:
                    logger.warning(
                        f"APBS failed for {pdb_id} pocket {pocket_num}: {stderr}"
                    )
                    continue

                # Verify output
                expected_output = (
                    output_dir / f"{pdb_id}_pkt_{pocket_num}_processed_pot.dx"
                )
                if not expected_output.exists():
                    logger.warning(
                        f"APBS did not generate expected output for {pdb_id} pocket {pocket_num}"
                    )
                else:
                    processed_count += 1

            logger.info(
                f"APBS completed for {pdb_id}: {processed_count}/{num_pockets} pockets processed"
            )
            return True

    except Exception as e:
        logger.error(f"Error in APBS for {pdb_id}: {str(e)}")
        raise ExternalToolError(f"APBS failed for {pdb_id}: {str(e)}")


@log_performance()
def parse_dx_file(dx_file_path: Union[str, Path]) -> np.ndarray:
    """Parse DX file and return potential values"""
    try:
        dx_file_path = Path(dx_file_path)
        validate_file_exists(dx_file_path, "DX file")

        logger.debug(f"Parsing DX file: {dx_file_path}")
        pot_grid = Grid(str(dx_file_path))

        logger.debug(f"DX grid shape: {pot_grid.grid.shape}")
        return pot_grid.grid

    except Exception as e:
        logger.error(f"Error parsing DX file {dx_file_path}: {str(e)}")
        raise ProteinAnalysisError(f"Failed to parse DX file: {str(e)}")


def parse_dx_file_for_view(dx_file_name: str, main_pdb_id: str) -> str:
    """Parse DX file and return content for visualization"""
    try:
        dx_file = APBS_DIR / main_pdb_id / dx_file_name
        dx_file_path = Path(dx_file)
        validate_file_exists(dx_file_path, "DX file")

        with open(dx_file_path, "r") as f:
            dx_content = f.read()

        logger.debug(
            f"Read DX file for visualization: {dx_file_path} ({len(dx_content)} chars)"
        )
        return dx_content

    except Exception as e:
        logger.error(f"Error parsing DX file {dx_file_path}: {str(e)}")
        raise ProteinAnalysisError(f"Failed to parse DX file: {str(e)}")


@log_analysis_step("Compare Pocket Electrostatic Potentials")
def compare_pockets_potentials(
    pdb_1_id: str, pdb_2_id: str, pdb_1_num_pockets: int, pdb_2_num_pockets: int
) -> pd.DataFrame:
    """Compare electrostatic potentials between pockets"""
    try:
        with log_context(tool="compare_potentials", pdb_1=pdb_1_id, pdb_2=pdb_2_id):
            # Collect all pocket potentials
            pocket_potentials = {}

            logger.info(
                f"Comparing electrostatic potentials: {pdb_1_id} ({pdb_1_num_pockets} pockets) vs {pdb_2_id} ({pdb_2_num_pockets} pockets)"
            )

            # Process protein 1 pockets
            for pocket_num in range(1, pdb_1_num_pockets + 1):
                pocket_name = f"{pdb_1_id}_pkt_{pocket_num}"
                dx_file = APBS_DIR / pdb_1_id / f"{pocket_name}_processed_pot.dx"

                if dx_file.exists():
                    try:
                        pot_values = parse_dx_file(dx_file)
                        pocket_potentials[pocket_name] = pot_values
                    except Exception as e:
                        logger.warning(
                            f"Failed to load potential for {pocket_name}: {str(e)}"
                        )

            # Process protein 2 pockets
            for pocket_num in range(1, pdb_2_num_pockets + 1):
                pocket_name = f"{pdb_2_id}_pkt_{pocket_num}"
                dx_file = APBS_DIR / pdb_2_id / f"{pocket_name}_processed_pot.dx"

                if dx_file.exists():
                    try:
                        pot_values = parse_dx_file(dx_file)
                        pocket_potentials[pocket_name] = pot_values
                    except Exception as e:
                        logger.warning(
                            f"Failed to load potential for {pocket_name}: {str(e)}"
                        )

            if not pocket_potentials:
                raise ProteinAnalysisError("No pocket potentials could be loaded")

            logger.info(f"Loaded {len(pocket_potentials)} pocket potentials")

            # Calculate similarity matrix
            pocket_names = list(pocket_potentials.keys())
            n_pockets = len(pocket_names)
            similarity_matrix = np.zeros((n_pockets, n_pockets))

            for i, pocket_1 in enumerate(pocket_names):
                for j, pocket_2 in enumerate(pocket_names):
                    if i == j:
                        similarity_matrix[i, j] = 1.0
                    elif i > j:
                        similarity_matrix[i, j] = similarity_matrix[j, i]
                    else:
                        try:
                            pot_1 = pocket_potentials[pocket_1].flatten()
                            pot_2 = pocket_potentials[pocket_2].flatten()

                            # Calculate correlation coefficient
                            correlation_matrix = np.corrcoef(pot_1, pot_2)
                            similarity_score = correlation_matrix[0, 1]

                            # Handle NaN values
                            if np.isnan(similarity_score):
                                similarity_score = 0.0

                            similarity_matrix[i, j] = similarity_score

                        except Exception as e:
                            logger.warning(
                                f"Error calculating similarity between {pocket_1} and {pocket_2}: {str(e)}"
                            )
                            similarity_matrix[i, j] = 0.0

            # Create DataFrame
            similarity_df = pd.DataFrame(
                similarity_matrix, index=pocket_names, columns=pocket_names
            )

            # Save results
            output_file = (
                APBS_DIR / f"{pdb_1_id}_{pdb_2_id}_pocket_corr_similarity_matrix.csv"
            )
            similarity_df.to_csv(output_file)

            logger.info(
                f"Calculated similarity matrix for {len(pocket_names)} pockets and saved to {output_file}"
            )
            return similarity_df

    except Exception as e:
        logger.error(f"Error comparing pocket potentials: {str(e)}")
        raise ProteinAnalysisError(f"Failed to compare pocket potentials: {str(e)}")


@log_performance()
def heatmap_pocket_corr(
    similarity_matrix: pd.DataFrame, pdb_1_id: str, pdb_2_id: str
) -> plt.Figure:
    """Create heatmap of pocket correlation matrix"""
    try:
        logger.info(f"Creating heatmap for pocket correlation matrix")

        # Create figure
        fig, ax = plt.subplots(figsize=(12, 10))

        # Define colormap
        cmap = LinearSegmentedColormap.from_list(
            "electrostatic", ["red", "white", "blue"], N=256
        )

        # Create heatmap
        heatmap = sns.heatmap(
            similarity_matrix,
            cmap=cmap,
            vmin=-1,
            vmax=1,
            center=0,
            fmt=".3f",
            linewidths=0.5,
            square=True,
            ax=ax,
            cbar_kws={"label": "Electrostatic Potential Similarity", "shrink": 0.8},
        )

        # Customize plot
        ax.set_title(
            f"Pocket Electrostatic Potential Similarity\n{pdb_1_id} vs {pdb_2_id}",
            fontsize=14,
            fontweight="bold",
            pad=20,
        )

        ax.set_xlabel("Pockets", fontsize=12)
        ax.set_ylabel("Pockets", fontsize=12)

        # Rotate labels for better readability
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

        plt.tight_layout()

        # Save figure
        output_file = APBS_DIR / f"{pdb_1_id}_{pdb_2_id}_pocket_corr_heatmap.png"
        plt.savefig(output_file, dpi=300, bbox_inches="tight")

        logger.info(f"Heatmap saved to {output_file}")
        return fig

    except Exception as e:
        logger.error(f"Error creating heatmap: {str(e)}")
        raise ProteinAnalysisError(f"Failed to create heatmap: {str(e)}")


@log_performance()
def cluster_pockets(
    similarity_matrix: pd.DataFrame,
    pdb_1_id: str,
    pdb_2_id: str,
    method: str = "ward",
    optimize_ordering: bool = True,
) -> plt.Figure:
    """Create dendrogram for pocket clustering"""
    try:
        logger.info(f"Creating dendrogram for pocket clustering using {method} method")

        # Validate input
        if similarity_matrix.empty:
            raise ProteinAnalysisError("Similarity matrix is empty")

        # Convert similarity to distance matrix
        distance_matrix = 1 - similarity_matrix.values

        # Handle invalid values
        distance_matrix = np.nan_to_num(
            distance_matrix, nan=1.0, posinf=1.0, neginf=0.0
        )

        # Ensure matrix is symmetric
        distance_matrix = (distance_matrix + distance_matrix.T) / 2

        # Convert to condensed distance matrix
        try:
            condensed_dist = squareform(distance_matrix, checks=False)
        except Exception as e:
            logger.warning(f"Error creating condensed distance matrix: {str(e)}")
            # Fallback: use upper triangle
            n = distance_matrix.shape[0]
            condensed_dist = []
            for i in range(n):
                for j in range(i + 1, n):
                    condensed_dist.append(distance_matrix[i, j])
            condensed_dist = np.array(condensed_dist)

        # Perform hierarchical clustering
        try:
            linkage_matrix = linkage(condensed_dist, method=method)

            # Apply optimal leaf ordering if requested
            if optimize_ordering:
                try:
                    linkage_matrix = optimal_leaf_ordering(
                        linkage_matrix, condensed_dist
                    )
                except Exception as e:
                    logger.warning(f"Optimal leaf ordering failed: {str(e)}")
        except Exception as e:
            logger.warning(f"Clustering failed with method {method}: {str(e)}")
            # Fallback to single linkage
            linkage_matrix = linkage(condensed_dist, method="single")

        # Create dendrogram
        fig, ax = plt.subplots(figsize=(14, 8))

        # Set color palette
        mpl.rcParams["axes.prop_cycle"] = mpl.cycler(color=plt.cm.tab20.colors)

        dendrogram(
            linkage_matrix,
            labels=similarity_matrix.index.tolist(),
            leaf_rotation=90,
            leaf_font_size=10,
            ax=ax,
            color_threshold=None,
            above_threshold_color="gray",
        )

        # Customize plot
        ax.set_title(
            f"Hierarchical Clustering of Pocket Electrostatic Similarity\n{pdb_1_id} vs {pdb_2_id}",
            fontsize=14,
            fontweight="bold",
            pad=20,
        )

        ax.set_xlabel("Pockets", fontsize=12)
        ax.set_ylabel("Distance (1 - Similarity)", fontsize=12)

        # Add method annotation
        ax.text(
            0.02,
            0.98,
            f"Method: {method}",
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )

        # Add grid
        ax.yaxis.grid(True, linestyle="--", alpha=0.7)
        ax.set_axisbelow(True)

        plt.tight_layout()

        # Save figure
        output_file = APBS_DIR / f"{pdb_1_id}_{pdb_2_id}_pocket_clustering.png"
        plt.savefig(output_file, dpi=300, bbox_inches="tight")

        logger.info(f"Dendrogram saved to {output_file}")
        return fig

    except Exception as e:
        logger.error(f"Error creating dendrogram: {str(e)}")
        raise ProteinAnalysisError(f"Failed to create dendrogram: {str(e)}")


###### Utility functions for data export and analysis


def export_analysis_summary(
    pdb_1_id: str, pdb_2_id: str, pdb_1_num_pockets: int, pdb_2_num_pockets: int
) -> dict:
    """Export comprehensive analysis summary"""
    try:
        summary = {
            "proteins": {
                "protein_1": {"id": pdb_1_id, "num_pockets": pdb_1_num_pockets},
                "protein_2": {"id": pdb_2_id, "num_pockets": pdb_2_num_pockets},
            },
            "analysis_folders": {
                "structural_alignment": str(FATCAT_DIR),
                "pocket_comparison": str(APOC_DIR),
                "similarity_analysis": str(APBS_DIR),
                "overlap_scores": str(OVERLAP_DIR),
                "pocket_prediction": str(P2RANK_DIR),
            },
            "timestamp": datetime.now().isoformat(),
            "version": "1.0.0",
        }

        logger.debug(f"Created analysis summary for {pdb_1_id} vs {pdb_2_id}")
        return summary

    except Exception as e:
        logger.error(f"Error creating analysis summary: {str(e)}")
        return {}


@log_performance()
def add_folder_to_zip(zip_file, folder_path, folder_name_in_zip):
    """Recursively add a folder and all its contents to a ZIP file"""
    folder_path = Path(folder_path)

    if not folder_path.exists():
        logger.warning(f"Folder does not exist: {folder_path}")
        return False

    files_added = 0
    # Add all files and subdirectories recursively
    for root, dirs, files in os.walk(folder_path):
        root_path = Path(root)

        # Calculate relative path from the original folder
        rel_path = root_path.relative_to(folder_path)

        # Create directory structure in ZIP
        if str(rel_path) != ".":
            zip_path = f"{folder_name_in_zip}/{rel_path}"
        else:
            zip_path = folder_name_in_zip

        # Add files
        for file in files:
            file_path = root_path / file
            if str(rel_path) != ".":
                archive_path = f"{folder_name_in_zip}/{rel_path}/{file}"
            else:
                archive_path = f"{folder_name_in_zip}/{file}"

            try:
                zip_file.write(file_path, archive_path)
                files_added += 1
            except Exception as e:
                logger.warning(f"Could not add {file_path} to ZIP: {str(e)}")

    logger.debug(f"Added {files_added} files from {folder_path} to ZIP")
    return True


@log_performance()
def create_results_package(pdb_1_id, pdb_2_id, pdb_1_num_pocket, pdb_2_num_pocket):
    """Create a ZIP of the entire results directory"""
    try:
        with st.spinner("📦 Creating results package..."):
            logger.info(f"Creating results package for {pdb_1_id} vs {pdb_2_id}")

            summary = export_analysis_summary(
                pdb_1_id, pdb_2_id, pdb_1_num_pocket, pdb_2_num_pocket
            )

            zip_buffer = io.BytesIO()

            with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
                # Add summary at root level
                summary_json = json.dumps(summary, indent=2)
                zip_file.writestr("analysis_summary.json", summary_json)
                logger.debug("Added analysis summary to ZIP")

                # Add the entire results directory
                if RESULTS_BASE_DIR.exists():
                    add_folder_to_zip(zip_file, RESULTS_BASE_DIR, "results")
                    st.success("✅ Results package created successfully!")
                    logger.info(f"Results package created successfully")
                else:
                    st.error("❌ Results directory not found")
                    logger.error(f"Results directory not found: {RESULTS_BASE_DIR}")
                    return None

            zip_buffer.seek(0)
            package_size = len(zip_buffer.getvalue())
            logger.info(f"Results package size: {package_size / 1024 / 1024:.2f} MB")
            return zip_buffer.getvalue()

    except Exception as e:
        st.error(f"❌ Error creating results package: {str(e)}")
        logger.error(f"Error creating results package: {str(e)}", exc_info=True)
        return None


# Log module import
logger.info("protFuncs module loaded successfully")
