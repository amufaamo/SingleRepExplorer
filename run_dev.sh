#!/bin/bash
# 開発用高速起動スクリプト（ビルドなし）
# docker build をスキップし、既存のローカルイメージ or Docker Hub イメージで即起動
#
# Usage:
#   ./run_dev.sh              # latest イメージで起動
#   ./run_dev.sh v1.0.0       # 指定バージョンで起動
#   ./run_dev.sh --pull       # Docker Hub から pull してから起動
set -euo pipefail

IMAGE_NAME="amufaamo/singlerepexplorer"
TAG="latest"
DO_PULL=false

for arg in "$@"; do
    case "$arg" in
        --pull) DO_PULL=true ;;
        v*) TAG="$arg" ;;
        *) echo "Unknown option: $arg"; exit 1 ;;
    esac
done

IMAGE="${IMAGE_NAME}:${TAG}"

echo "=== SingleRepExplorer 高速起動 ==="
echo "  Image: ${IMAGE}"
echo ""

# pull オプションが指定された場合のみ Docker Hub から取得
if $DO_PULL; then
    echo ">>> Pulling ${IMAGE} from Docker Hub..."
    docker pull "${IMAGE}"
fi

# 既存コンテナを停止・削除（エラーは無視）
docker rm -f singlerepexplorer 2>/dev/null || true

# コンテナ起動（ビルドなし・app ディレクトリをマウントして即反映）
docker run -d \
    --name singlerepexplorer \
    -p 3838:3838 \
    -v "$(pwd)/app:/srv/shiny-server" \
    -e SHINY_LOG_LEVEL=INFO \
    --memory=8g \
    --restart unless-stopped \
    "${IMAGE}"

echo ""
echo "==================================================================="
echo "  起動完了!"
echo "  URL     : http://localhost:3838/"
echo "  ログ確認: docker logs -f singlerepexplorer"
echo "  停止    : docker stop singlerepexplorer"
echo "==================================================================="
