#!/bin/bash
docker run -it --rm \
    --name dss-methylation-terminal \
    -v $(pwd)/data:/workspace/data \
    -v $(pwd)/results:/workspace/results \
    -v $(pwd)/scripts:/workspace/scripts \
    -v $(pwd)/notebooks:/workspace/notebooks \
    redndgreen8/methsmootheval:latest bash
