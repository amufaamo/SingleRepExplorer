#!/bin/bash
echo "Building amufaamo/singlerepexplorer:v1.1.17 and latest..."
docker build -t amufaamo/singlerepexplorer:v1.1.17 -t amufaamo/singlerepexplorer:latest .

echo "Pushing amufaamo/singlerepexplorer:v1.1.17..."
docker push amufaamo/singlerepexplorer:v1.1.17

echo "Pushing amufaamo/singlerepexplorer:latest..."
docker push amufaamo/singlerepexplorer:latest

echo "Done!"
