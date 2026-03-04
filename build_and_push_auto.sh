#!/bin/bash
# ==============================================================================
# SingleRepExplorer Automatic Build & Push Script
# Extracts version from app.R and pushes tags to Docker Hub
# ==============================================================================

# 1. app.Rのタイトル行からバージョン番号を自動抽出 (例: v1.1.17)
# 例: title = "SingleRepExplorer v1.1.17: Web Application...
VERSION=$(grep 'title = "SingleRepExplorer v' app/app.R | sed -n 's/.*SingleRepExplorer \([vV][0-9]*\.[0-9]*\.[0-9]*\):.*/\1/p')

# バージョンが取得できなかった場合のフォールバック（手動入力）
if [ -z "$VERSION" ]; then
    echo "Warning: Automatically extracting version from app/app.R failed."
    read -p "Please enter the version manually (e.g., v1.1.17): " VERSION
fi

# 念のための確認プロンプト
echo "Detected Version: ${VERSION}"
read -p "Are you sure you want to build and push this version? (y/n): " confirm
if [ "$confirm" != "y" ]; then
    echo "Build aborted."
    exit 1
fi

echo "=================================================="
echo "Building amufaamo/singlerepexplorer:${VERSION} and latest..."
echo "=================================================="

# Docker ビルド
docker build -t amufaamo/singlerepexplorer:${VERSION} -t amufaamo/singlerepexplorer:latest .

# ビルド成功確認
if [ $? -eq 0 ]; then
    echo "Build successful! Proceeding to push..."
    
    echo "Pushing amufaamo/singlerepexplorer:${VERSION}..."
    docker push amufaamo/singlerepexplorer:${VERSION}
    
    echo "Pushing amufaamo/singlerepexplorer:latest..."
    docker push amufaamo/singlerepexplorer:latest
    
    echo "Done! The image has been pushed to Docker Hub."
else
    echo "Docker build failed. Push aborted."
    exit 1
fi
