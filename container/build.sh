#!/bin/bash

# Build script for DNA Barcoding Analysis Container

set -e

echo "===================================================="
echo "Building DNA Barcoding Analysis Container"
echo "===================================================="
echo ""

# Check if Docker is running
if ! docker info > /dev/null 2>&1; then
    echo "ERROR: Docker is not running!"
    echo "Please start Docker Desktop and try again."
    exit 1
fi

# Build the container
echo "Building container..."
docker build -t dna-barcoding:latest container/

# Check the size
echo ""
echo "===================================================="
echo "Build complete! Checking container size..."
echo "===================================================="
docker images dna-barcoding:latest

# Get the size in MB
SIZE=$(docker images dna-barcoding:latest --format "{{.Size}}")
echo ""
echo "Container size: $SIZE"
echo ""

# Verify tools are installed
echo "===================================================="
echo "Verifying installed tools..."
echo "===================================================="

echo "Checking MAFFT..."
docker run --rm dna-barcoding:latest bash -c "mafft --version" 2>&1 | head -n 1

echo "Checking IQ-TREE..."
docker run --rm dna-barcoding:latest bash -c "iqtree --version" 2>&1 | head -n 1

echo "Checking BLAST..."
docker run --rm dna-barcoding:latest bash -c "blastn -version" 2>&1 | head -n 1

echo "Checking Python..."
docker run --rm dna-barcoding:latest bash -c "python --version"

echo "Checking BioPython..."
docker run --rm dna-barcoding:latest bash -c "python -c 'import Bio; print(f\"BioPython {Bio.__version__}\")'"

echo "Checking R..."
docker run --rm dna-barcoding:latest bash -c "R --version" 2>&1 | head -n 1

echo "Checking ape package..."
docker run --rm dna-barcoding:latest bash -c "R -e 'packageVersion(\"ape\")' 2>&1" | grep -A 1 "packageVersion"

echo ""
echo "===================================================="
echo "Build successful!"
echo "===================================================="
echo ""
echo "To run the container:"
echo "  docker run -it -v \$(pwd):/workspace dna-barcoding:latest"
echo ""
echo "Or use docker-compose:"
echo "  docker-compose up -d"
echo "  docker exec -it dna-barcoding zsh"
echo ""
