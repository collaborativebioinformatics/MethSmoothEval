#!/bin/bash
# run_dss_container.sh - Run the DSS methylation analysis container

# Configuration
CONTAINER_NAME="dss-methylation-lab"
IMAGE_NAME="dss-methylation:latest2"
HOST_PORT=8888
CONTAINER_PORT=8888

# Create local directories for persistent storage
mkdir -p ./data ./results ./scripts

echo "Starting DSS Methylation Analysis Container..."
echo "============================================="

# Check if container is already running
if [ $(docker ps -q -f name=$CONTAINER_NAME) ]; then
    echo "⚠️  Container $CONTAINER_NAME is already running"
    echo "Stopping existing container..."
    docker stop $CONTAINER_NAME
    docker rm $CONTAINER_NAME
fi

# Run the container with volume mounts
docker run -d \
    --name $CONTAINER_NAME \
    -p $HOST_PORT:$CONTAINER_PORT \
    -v $(pwd)/data:/workspace/data \
    -v $(pwd)/results:/workspace/results \
    -v $(pwd)/scripts:/workspace/scripts \
    $IMAGE_NAME

if [ $? -eq 0 ]; then
    echo "✅ Container started successfully!"
    echo ""
    echo "Container Details:"
    echo "=================="
    echo "Name: $CONTAINER_NAME"
    echo "Port: $HOST_PORT"
    echo "JupyterLab URL: http://localhost:$HOST_PORT"
    echo ""
    echo "Mounted Volumes:"
    echo "- ./data -> /workspace/data"
    echo "- ./results -> /workspace/results"
    echo "- ./scripts -> /workspace/scripts"
    echo ""
    echo "Useful Commands:"
    echo "- View logs: docker logs $CONTAINER_NAME"
    echo "- Enter container: docker exec -it $CONTAINER_NAME bash"
    echo "- Stop container: docker stop $CONTAINER_NAME"
    echo ""
    echo "Waiting for JupyterLab to start..."
    sleep 5
    
    # Check if JupyterLab is responding
    if curl -s http://localhost:$HOST_PORT > /dev/null; then
        echo "🚀 JupyterLab is ready at: http://localhost:$HOST_PORT"
    else
        echo "⏳ JupyterLab is still starting up. Check logs with: docker logs $CONTAINER_NAME"
    fi
else
    echo "❌ Failed to start container!"
    exit 1
fi

