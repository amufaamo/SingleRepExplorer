#!/bin/bash
echo "Building amufaamo/singlerepexplorer:v1.1.15 and latest..."
docker build -t amufaamo/singlerepexplorer:v1.1.15 -t amufaamo/singlerepexplorer:latest .

echo "Pushing amufaamo/singlerepexplorer:v1.1.15..."
docker push amufaamo/singlerepexplorer:v1.1.15

echo "Pushing amufaamo/singlerepexplorer:latest..."
docker push amufaamo/singlerepexplorer:latest

echo "Done!"
