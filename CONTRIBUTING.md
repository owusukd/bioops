# Contributing to BioOps

Thank you for your interest in contributing to BioOps! This document provides guidelines and instructions for contributing to the project.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [How to Contribute](#how-to-contribute)
- [Pull Request Process](#pull-request-process)
- [Style Guidelines](#style-guidelines)
- [Testing](#testing)
- [Documentation](#documentation)
- [Issue Reporting](#issue-reporting)

## Code of Conduct

By participating in this project, you agree to abide by our Code of Conduct:

- Be respectful and inclusive
- Welcome newcomers and help them get started
- Focus on constructive criticism
- Respect differing viewpoints and experiences
- Show empathy towards other community members

## Getting Started

1. **Fork the Repository**
   ```bash
   # Fork on GitHub, then:
   git clone https://github.com/owusukd/bioops.git
   cd bioops
   ```

2. **Create a Development Branch**
   ```bash
   git checkout -b feature/your-feature-name
   # or
   git checkout -b fix/issue-description
   ```

3. **Set Up Development Environment**
   ```bash
   # Create virtual environment
   python -m venv venv-dev
   source venv-dev/bin/activate  # On Windows: venv-dev\Scripts\activate
   
   # Install dependencies
   pip install -r requirements.txt
   pip install -r requirements-dev.txt
   ```

## Development Setup

### Prerequisites

- Python 3.12+
- Docker and Docker Compose (optional)
- Git
- Make (optional, for using Makefile commands)

### Installing Development Dependencies

```bash
pip install -r requirements-dev.txt
```

This includes:
- `pytest` - Testing framework
- `pytest-cov` - Coverage reporting
- `black` - Code formatting
- `flake8` - Linting
- `mypy` - Type checking
- `pre-commit` - Git hooks

### Setting Up Pre-commit Hooks

```bash
pre-commit install
```

This will run formatting and linting checks before each commit.

## How to Contribute

### Types of Contributions

1. **Bug Fixes**
   - Fix bugs reported in GitHub Issues
   - Add tests to prevent regression

2. **New Features**
   - Discuss the feature in an issue first
   - Implement with tests and documentation

3. **Documentation**
   - Improve existing documentation
   - Add examples and tutorials
   - Fix typos and clarify instructions

4. **Performance Improvements**
   - Optimize existing code
   - Add benchmarks to demonstrate improvements

5. **Test Coverage**
   - Add missing tests
   - Improve test quality

### Contribution Workflow

1. **Check Existing Issues**
   - Look for issues tagged with `good first issue` or `help wanted`
   - Comment on the issue to indicate you're working on it

2. **Create/Update Your Fork**
   ```bash
   git remote add upstream https://github.com/owusukd/bioops.git
   git fetch upstream
   git checkout main
   git merge upstream/main
   ```

3. **Make Your Changes**
   - Write clean, readable code
   - Follow the style guidelines
   - Add tests for new functionality
   - Update documentation as needed

4. **Test Your Changes**
   ```bash
   # Run tests
   make test
   
   # Check code style
   make lint
   
   # Format code
   make format
   ```

5. **Commit Your Changes**
   ```bash
   git add .
   git commit -m "Brief description of changes"
   ```

   Commit message format:
   ```
   type: subject
   
   body (optional)
   
   footer (optional)
   ```
   
   Types: `feat`, `fix`, `docs`, `style`, `refactor`, `test`, `chore`

6. **Push to Your Fork**
   ```bash
   git push origin feature/your-feature-name
   ```

## Pull Request Process

1. **Create a Pull Request**
   - Use a clear, descriptive title
   - Reference any related issues
   - Describe what changes you made and why
   - Include screenshots for UI changes

2. **Pull Request Template**
   ```markdown
   ## Description
   Brief description of changes
   
   ## Related Issue
   Fixes #(issue number)
   
   ## Type of Change
   - [ ] Bug fix
   - [ ] New feature
   - [ ] Documentation update
   - [ ] Performance improvement
   
   ## Testing
   - [ ] Tests pass locally
   - [ ] New tests added
   - [ ] Documentation updated
   
   ## Screenshots (if applicable)
   ```

3. **Review Process**
   - Maintainers will review your PR
   - Address any requested changes
   - Once approved, your PR will be merged

## Style Guidelines

### Python Code Style

We use [PEP 8](https://pep8.org/) with the following specifications:

- Line length: 120 characters
- Use `black` for formatting
- Use type hints where possible
- Docstrings for all public functions/classes

Example:
```python
from typing import List, Optional

def analyze_proteins(
    pdb_id_1: str,
    pdb_id_2: str,
    chain_ids: Optional[List[str]] = None
) -> dict:
    """
    Analyze two protein structures.
    
    Args:
        pdb_id_1: First protein PDB identifier
        pdb_id_2: Second protein PDB identifier
        chain_ids: Optional list of chain identifiers
        
    Returns:
        Dictionary containing analysis results
        
    Raises:
        ValueError: If PDB IDs are invalid
    """
    # Implementation
    pass
```

### Import Organization

```python
# Standard library imports
import os
import sys
from pathlib import Path

# Third-party imports
import numpy as np
import pandas as pd
import streamlit as st

# Local imports
from bioops import config
from protFuncs import analyze_structure
```

### Error Handling

```python
try:
    result = risky_operation()
except SpecificError as e:
    logger.error(f"Operation failed: {e}")
    raise AnalysisError(f"Could not complete analysis: {e}") from e
```

## Testing

### Writing Tests

- Place tests in the `tests/` directory
- Mirror the source code structure
- Use descriptive test names
- Include both positive and negative test cases

Example:
```python
# tests/test_protFuncs.py
import pytest
from protFuncs import validate_pdb_id

class TestValidation:
    def test_valid_pdb_id(self):
        """Test validation of correct PDB ID format"""
        assert validate_pdb_id("1ABC") == True
        assert validate_pdb_id("5CXV") == True
    
    def test_invalid_pdb_id(self):
        """Test validation rejects invalid formats"""
        assert validate_pdb_id("ABC") == False
        assert validate_pdb_id("12345") == False
        assert validate_pdb_id("") == False
    
    def test_pdb_id_case_insensitive(self):
        """Test PDB ID validation is case-insensitive"""
        assert validate_pdb_id("1abc") == True
        assert validate_pdb_id("1ABC") == True
```

### Running Tests

```bash
# Run all tests
make test

# Run specific test file
pytest tests/test_protFuncs.py

# Run with coverage
pytest --cov=bioops --cov=protFuncs --cov-report=html

# Run only unit tests
make test-unit

# Run only integration tests
make test-integration
```

## Documentation

### Code Documentation

- Add docstrings to all public functions and classes
- Use Google-style docstrings
- Include type hints
- Provide examples in docstrings when helpful

### User Documentation

- Update README.md for significant changes
- Add/update entries in docs/ directory
- Include screenshots for UI changes
- Provide clear examples

### API Documentation

- Document all public APIs
- Include request/response examples
- Note any breaking changes

## Issue Reporting

### Bug Reports

Include:
- Clear description of the bug
- Steps to reproduce
- Expected behavior
- Actual behavior
- System information (OS, Python version, etc.)
- Error messages and stack traces
- Screenshots if applicable

### Feature Requests

Include:
- Clear description of the feature
- Use case / motivation
- Proposed implementation (optional)
- Alternative solutions considered

### Good First Issues

If you're new to the project, look for issues labeled:
- `good first issue` - Simple fixes suitable for beginners
- `documentation` - Documentation improvements
- `help wanted` - Issues where we need community help

## Questions?

If you have questions about contributing:
1. Check existing documentation
2. Search closed issues
3. Ask in an issue or discussion
4. Contact the maintainers

Thank you for contributing to BioOps! 🧬
