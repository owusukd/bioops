SHELL := /bin/bash
# BioOps Makefile
# Simplifies common development and deployment tasks

.PHONY: help build run stop clean test lint format docs install dev prod monitoring

# Default target
help:
	@echo "BioOps Development and Deployment Commands"
	@echo "=========================================="
	@echo "Development:"
	@echo "  make install      - Install dependencies and tools"
	@echo "  make dev          - Run in development mode"
	@echo "  make test         - Run tests"
	@echo "  make lint         - Run linting"
	@echo "  make format       - Format code"
	@echo "  make docs         - Build documentation"
	@echo ""
	@echo "Docker:"
	@echo "  make build        - Build Docker image"
	@echo "  make run          - Run with Docker Compose"
	@echo "  make stop         - Stop all containers"
	@echo "  make clean        - Clean up containers and volumes"
	@echo "  make prod         - Run in production mode"
	@echo "  make monitoring   - Run with monitoring stack"
	@echo ""
	@echo "Maintenance:"
	@echo "  make logs         - Show application logs"
	@echo "  make shell        - Open shell in container"
	@echo "  make backup       - Backup analysis results"
	@echo "  make restore      - Restore from backup"

# Installation
install:
	@echo "Installing BioOps dependencies..."
	chmod +x install_tools
	./install_tools
	@echo "Installation complete!"

# Development
dev:
	@echo "Starting BioOps in development mode..."
	source .bioops/bin/activate && streamlit run bioops.py --server.port=8501 --server.enableCORS=false --server.enableXsrfProtection=false

# Testing
test:
	@echo "Running tests..."
	source .bioops/bin/activate && python -m pytest tests/ -v --cov=bioops --cov=protFuncs --cov-report=html

test-unit:
	@echo "Running unit tests..."
	source .bioops/bin/activate && python -m pytest tests/unit/ -v

test-integration:
	@echo "Running integration tests..."
	source .bioops/bin/activate && python -m pytest tests/integration/ -v

# Code quality
lint:
	@echo "Running linting..."
	source .bioops/bin/activate && flake8 bioops.py protFuncs.py --max-line-length=120 --ignore=E203,W503
	source .bioops/bin/activate && mypy bioops.py protFuncs.py --ignore-missing-imports

format:
	@echo "Formatting code..."
	source .bioops/bin/activate && black bioops.py protFuncs.py
	source .bioops/bin/activate && isort bioops.py protFuncs.py

# Documentation
docs:
	@echo "Building documentation..."
	source .bioops/bin/activate && mkdocs build

docs-serve:
	@echo "Serving documentation..."
	source .bioops/bin/activate && mkdocs serve

# Docker operations
build:
	@echo "Building Docker image..."
	docker-compose build

run:
	@echo "Starting BioOps with Docker Compose..."
	docker-compose up -d
	@echo "BioOps is running at http://localhost:8501"

stop:
	@echo "Stopping all containers..."
	docker-compose down

clean:
	@echo "Cleaning up containers and volumes..."
	docker-compose down -v --remove-orphans
	rm -rf tempDownloadDir/results/*
	rm -rf logs/*

# Production deployment
prod:
	@echo "Starting BioOps in production mode..."
	docker-compose --profile production up -d

monitoring:
	@echo "Starting BioOps with monitoring stack..."
	docker-compose --profile monitoring up -d
	@echo "Services:"
	@echo "  - BioOps: http://localhost:8501"
	@echo "  - Grafana: http://localhost:3000"
	@echo "  - Prometheus: http://localhost:9090"

# Logs and debugging
logs:
	docker-compose logs -f bioops

logs-tail:
	docker-compose logs --tail=100 bioops

shell:
	docker-compose exec bioops /bin/bash

# Backup and restore
backup:
	@echo "Creating backup..."
	mkdir -p backups
	tar -czf backups/bioops-backup-$$(date +%Y%m%d-%H%M%S).tar.gz tempDownloadDir/results logs
	@echo "Backup created in backups/"

restore:
	@echo "Available backups:"
	@ls -1 backups/*.tar.gz
	@echo ""
	@echo "To restore, run: tar -xzf backups/[backup-file] -C ."

# Environment setup
env:
	@echo "Setting up environment..."
	cp .env.example .env
	@echo "Environment file created. Please edit .env with your settings."

# Health check
health:
	@echo "Checking service health..."
	@curl -f http://localhost:8501/_stcore/health || echo "Service is not healthy"

# Performance monitoring
perf:
	@echo "Checking performance metrics..."
	@curl -s http://localhost:9090/metrics | grep -E "bioops_|process_" | head -20

# Database migrations (for future use)
migrate:
	@echo "Running database migrations..."
	@echo "No migrations to run (database not yet implemented)"

# Security scanning
security:
	@echo "Running security scan..."
	source .bioops/bin/activate && pip-audit
	source .bioops/bin/activate && bandit -r bioops.py protFuncs.py

# Docker image optimization
optimize:
	@echo "Analyzing Docker image size..."
	docker images bioops:latest
	@echo ""
	@echo "Running dive for detailed analysis..."
	dive bioops:latest

# Clean Python cache
clean-cache:
	@echo "Cleaning Python cache..."
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete
	rm -rf .pytest_cache
	rm -rf .mypy_cache
	rm -rf htmlcov
	rm -rf .coverage

# Full clean
clean-all: clean clean-cache
	@echo "Cleaning all generated files..."
	rm -rf .bioops
	rm -rf tools/
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info

# Version management
version:
	@echo "BioOps version information:"
	@grep -E "^__version__|^VERSION" bioops.py || echo "__version__ = '1.0.0'"
	@echo ""
	@echo "Tool versions:"
	@echo "Python: $$(python --version)"
	@echo "Streamlit: $$(streamlit --version)"
	@docker-compose exec bioops prank -version 2>/dev/null || echo "P2Rank: Not running"
	@docker-compose exec bioops FATCAT -help 2>&1 | head -1 || echo "FATCAT: Not running"
