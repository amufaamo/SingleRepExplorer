#!/bin/bash
echo "Building amufaamo/singlerepexplorer:v2.0.0 and latest..."
docker build -t amufaamo/singlerepexplorer:v2.0.0 -t amufaamo/singlerepexplorer:latest .

echo "Pushing amufaamo/singlerepexplorer:v2.0.0..."
docker push amufaamo/singlerepexplorer:v2.0.0

echo "Pushing amufaamo/singlerepexplorer:latest..."
docker push amufaamo/singlerepexplorer:latest

echo "Done!"
