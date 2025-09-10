#!/bin/bash
# build_dss_container.sh - Build the DSS methylation analysis container

echo "Building DSS Methylation Analysis Container..."
echo "============================================="

# Build the Docker image
docker build --platform linux/amd64 -t dss-methylation:latest .

if [ $? -eq 0 ]; then
    echo "✅ Container built successfully!"
    echo "Image: dss-methylation:latest"
    echo ""
    echo "To run the container, use:"
    echo "./run_dss_container.sh"
else
    echo "❌ Container build failed!"
    exit 1
fi
