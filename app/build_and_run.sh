#!/bin/bash

# SingleRepExplorer Docker ビルド・実行スクリプト

echo "=== SingleRepExplorer Docker Setup Script ==="

# アーキテクチャの検出
ARCH=$(uname -m)
echo "Detected architecture: $ARCH"

# Docker Composeで使用するDockerfile名を環境変数としてエクスポート
if [[ "$ARCH" == "arm64" ]]; then
    echo "Using ARM64 optimized build (Dockerfile.arm64)..."
    export DOCKERFILE_NAME="Dockerfile.arm64"
    # 注意: Dockerfile.arm64 が存在しない場合、このスクリプトはエラーになります。
    # 必要に応じて `cp Dockerfile Dockerfile.arm64` を実行するか、
    # このスクリプトの分岐を修正して常に 'Dockerfile' を使うようにしてください。
else
    echo "Using standard x86_64 build (Dockerfile)..."
    export DOCKERFILE_NAME="Dockerfile"
fi

# 既存のコンテナを停止・削除 (依存関係のない孤立したコンテナも削除)
echo "Stopping and removing existing containers..."
docker-compose down --remove-orphans 2>/dev/null || true # 存在しない場合にエラーにならないようにする

# Docker imageの削除（オプション）
# ユーザーに確認 (-r オプションで Enter キーを不要に、-n 1 で1文字だけ読み込む)
read -p "Remove existing Docker images for this project? (y/N): " -n 1 -r
echo # 改行を追加

if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Removing existing images..."
    docker-compose down --rmi all 2>/dev/null || true # composeで管理されるイメージを削除
fi

echo "Building and running with ${DOCKERFILE_NAME}..."
# docker-compose up でビルドと起動 (-d でバックグラウンド実行)
docker-compose up --build -d

echo "==================================================================="
echo "SingleRepExplorer Shiny App should be starting up..."
echo "Access URL: http://localhost:3838/"
echo "To view logs: docker-compose logs -f singlerepexplorer"
echo "To stop: docker-compose down"
echo "==================================================================="
