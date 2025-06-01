# BioOps - Protein Analysis Tool

<div align="center">
  <img src="logo/BioOps_Main_Logo.png" alt="BioOps Logo" width="400"/>
  
  [![Python](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/)
  [![Streamlit](https://img.shields.io/badge/streamlit-1.44-red.svg)](https://streamlit.io/)
  [![Docker](https://img.shields.io/badge/docker-ready-blue.svg)](https://www.docker.com/)
  [![License: CC BY-NC 4.0](https://licensebuttons.net/l/by-nc/4.0/88x31.png)](https://creativecommons.org/licenses/by-nc/4.0/)
</div>

## 🧬 Overview

BioOps is a comprehensive protein analysis platform specializing in Protein analysis. It integrates multiple state-of-the-art computational tools to provide detailed structural and electrostatic analysis of protein binding pockets.

### Key Features

- **3D Structural Comparison** using FATCAT algorithm
- **Binding Pocket Prediction** with P2Rank
- **Pocket Similarity Analysis** using APoc
- **Electrostatic Potential Calculation** with APBS/PDB2PQR
- **Interactive 3D Visualization** powered by py3Dmol
- **Comprehensive Results Export** with detailed analysis reports

## 🚀 Quick Start

### Using Docker (Recommended)

```bash
# Pull image from dockerhub
sudo docker pull owusukd/bioops

# Run with Docker
sudo docker run -d -p 8501:8501 --name bioops owusukd/bioops:latest

# Access at http://localhost:8501
```

### Using Make

```
BioOps Development and Deployment Commands
==========================================
Development:
  make install      - Install dependencies and tools
  make dev          - Run in development mode
  make test         - Run tests
  make lint         - Run linting
  make format       - Format code
  make docs         - Build documentation (not implemented)

Docker:
  make build        - Build Docker image
  make run          - Run with Docker Compose
  make stop         - Stop all containers
  make clean        - Clean up containers and volumes
  make prod         - Run in production mode
  make monitoring   - Run with monitoring stack (still in the works)

Maintenance:
  make logs         - Show application logs
  make shell        - Open shell in container
  make backup       - Backup analysis results
  make restore      - Restore from backup
```

### Manual Installation

```bash
# Clone the repository
git clone https://github.com/owusukd/bioops.git
cd bioops

# Run the installation script
chmod +x install_tools
./install_tools

# Activate the virtual environment
source .bioops/bin/activate

# Run the application
streamlit run bioops.py
```

## 📋 Requirements

### System Requirements
- Ubuntu 20.04+ (or compatible Linux distribution)
- Python 3.12+
- Java 17+ (for P2Rank)
- 8GB RAM minimum (16GB recommended)
- 10GB free disk space

### External Tools
The following tools are automatically installed by the setup script:
- **FATCAT** - Flexible structure alignment
- **P2Rank** - Machine learning-based pocket prediction
- **APoc** - Pocket structural comparison
- **APBS** - Electrostatic calculations
- **PDB2PQR** - PDB file preparation for APBS

## 🔧 Configuration

### Environment Variables

Create a `.env` file in the project root:

```env
# Logging Configuration
LOG_LEVEL=INFO
LOG_FORMAT=json
LOG_FILE_PATH=logs/bioops.log
LOG_MAX_SIZE=10485760  # 10MB
LOG_BACKUP_COUNT=5

# Analysis Configuration
ANALYSIS_TIMEOUT=900  # 15 minutes
MAX_FILE_SIZE=52428800  # 50MB
TEMP_DIR=tempDownloadDir

# Performance Configuration
MAX_WORKERS=8
CACHE_ENABLED=true
CACHE_TTL=3600  # 1 hour
```

### Logging Configuration

BioOps uses a comprehensive logging system. Configure logging in `config/logging.yaml`:

```yaml
version: 1
formatters:
  default:
    format: '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
  json:
    class: pythonjsonlogger.jsonlogger.JsonFormatter
    format: '%(asctime)s %(name)s %(levelname)s %(message)s'

handlers:
  console:
    class: logging.StreamHandler
    level: INFO
    formatter: default
    stream: ext://sys.stdout
  
  file:
    class: logging.handlers.RotatingFileHandler
    level: DEBUG
    formatter: json
    filename: logs/bioops.log
    maxBytes: 10485760  # 10MB
    backupCount: 5

loggers:
  bioops:
    level: DEBUG
    handlers: [console, file]
    propagate: false
  
  protFuncs:
    level: DEBUG
    handlers: [console, file]
    propagate: false

root:
  level: INFO
  handlers: [console, file]
```

## 📊 Usage

### Input Methods

1. **PDB ID Entry**
   - Enter 4-character PDB identifiers
   - Optionally specify chain IDs
   - Files are automatically downloaded from RCSB

2. **File Upload**
   - Upload PDB files directly
   - Maximum file size: 50MB
   - Supported format: .pdb

### Analysis Workflow

1. **Select Input Method** - Choose between PDB ID or file upload
2. **Provide Structures** - Enter IDs or upload files for two proteins
3. **Specify Chains** - Optional, defaults to chain 'A'
4. **Run Analysis** - Click "Start Analysis" to begin processing
5. **View Results** - Explore interactive visualizations and download reports

### Results Include

- **Structure Analysis Tab**
  - 3D protein visualizations
  - Structural alignment results
  - Pocket overlap scores
  - APoc pocket comparison metrics

- **Electrostatic Analysis Tab**
  - Similarity heatmaps
  - Hierarchical clustering dendrograms
  - Individual pocket potential visualizations

- **Summary Tab**
  - Analysis statistics
  - Downloadable results (CSV, ZIP package)

## 🏗️ Architecture

```
bioops/
├── bioops.py              # Main Streamlit application  
├── protFuncs.py           # Core analysis functions
├── install_tools          # Tool installation script
├── requirements.txt       # Python dependencies
├── Dockerfile            # Docker container definition
├── docker-compose.yml    # Docker Compose configuration
├── config/
│   ├── logging.yaml      # Logging configuration
│   └── analysis.yaml     # Analysis parameters
├── logs/                 # Application logs
├── tempDownloadDir/      # Temporary analysis files
│   └── results/          # Analysis results
│       ├── fatcat/       # FATCAT alignment results
│       ├── p2rank/       # P2Rank pocket predictions
│       ├── apoc/         # APoc comparison results
│       ├── apbs/         # APBS electrostatic calculations
│       └── overlap_score/ # Overlap score calculations
└── logo/                 # Application logos
```

## 🐛 Troubleshooting

### Common Issues

1. **Java Not Found**
   ```bash
   export JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
   export PATH=$JAVA_HOME/bin:$PATH
   ```

2. **Tool Command Not Found**
   ```bash
   source ~/.bashrc
   # Or manually add tool paths to PATH
   ```

3. **Permission Denied**
   ```bash
   chmod +x install_tools
   chmod +x tools/*/bin/*
   ```

4. **Memory Issues**
   - Increase Docker memory allocation
   - Reduce number of concurrent analyses
   - Process smaller proteins

### Log Files

Check logs for detailed error information:
- Application logs: `logs/bioops.log`
- Streamlit logs: `logs/streamlit.log`
- Tool-specific logs in `tempDownloadDir/results/*/`

## 🤝 Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Development Setup

```bash
# Clone the repository
git clone https://github.com/owusukd/bioops.git
cd bioops

# Create development environment
python -m venv venv-dev
source venv-dev/bin/activate

# Install development dependencies
pip install -r requirements-dev.txt

# Run tests
pytest tests/

# Run linting
flake8 bioops.py protFuncs.py
black --check bioops.py protFuncs.py
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

BioOps integrates the following tools:
- [FATCAT](https://fatcat.godziklab.org/) - Ye Y. & Godzik A. (2003)
- [P2Rank](https://github.com/rdk/p2rank) - Krivák R. & Hoksza D. (2018)
- [APoc](https://doi.org/10.1093/bioinformatics/btt024) - Gao M. & Skolnick J. (2013)
- [APBS](https://www.poissonboltzmann.org/) - Baker N.A. et al. (2001)
- [PDB2PQR](https://www.poissonboltzmann.org/) - Dolinsky T.J. et al. (2007)

## 📧 Contact

- **Developer**: Kwabena O. Dankwah
- **Email**: owusukd@yahoo.com
- **GitHub**: [@owusukd](https://github.com/owusukd)
- **Issues**: [GitHub Issues](https://github.com/owusukd/bioops/issues)

## 🗺️ Roadmap

- [ ] Support for AlphaFold structures
- [ ] Batch processing capabilities
- [ ] REST API implementation
- [ ] Machine learning-based pocket classification
- [ ] Integration with molecular dynamics simulations
- [ ] Support for protein-protein interaction analysis

---

<div align="center">
  Made with ❤️ for the bioinformatics community
</div>
