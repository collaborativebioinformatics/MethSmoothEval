#!/bin/bash
# manage_container.sh - Container management utilities

CONTAINER_NAME="dss-methylation-lab"

case "$1" in
    "start")
        echo "Starting container..."
        docker start $CONTAINER_NAME
        ;;
    "stop")
        echo "Stopping container..."
        docker stop $CONTAINER_NAME
        ;;
    "restart")
        echo "Restarting container..."
        docker restart $CONTAINER_NAME
        ;;
    "logs")
        echo "Showing container logs..."
        docker logs -f $CONTAINER_NAME
        ;;
    "shell")
        echo "Entering container shell..."
        docker exec -it $CONTAINER_NAME bash
        ;;
    "status")
        echo "Container status:"
        docker ps -f name=$CONTAINER_NAME
        ;;
    "clean")
        echo "Cleaning up container and image..."
        docker stop $CONTAINER_NAME 2>/dev/null
        docker rm $CONTAINER_NAME 2>/dev/null
        docker rmi dss-methylation:latest 2>/dev/null
        echo "Cleanup complete!"
        ;;
    *)
        echo "Usage: $0 {start|stop|restart|logs|shell|status|clean}"
        echo ""
        echo "Commands:"
        echo "  start   - Start the container"
        echo "  stop    - Stop the container"
        echo "  restart - Restart the container"
        echo "  logs    - Show container logs"
        echo "  shell   - Enter container shell"
        echo "  status  - Show container status"
        echo "  clean   - Remove container and image"
        ;;
esac