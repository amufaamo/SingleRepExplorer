#!/bin/bash
# ==============================================================================
# SingleRepExplorer Docker ビルド・実行スクリプト
#
# Usage:
#   ./build_and_run.sh [VERSION] [OPTIONS]
#
# Arguments:
#   VERSION       イメージにつけるバージョンタグ (例: v1.0.0)
#                 省略時は app/app.R から自動抽出
#
# Options:
#   --pull        ローカルビルドせず Docker Hub からイメージを pull して起動
#   --no-cache    Docker ビルド時にキャッシュを使わない
#   -h, --help    このヘルプを表示する
#
# Examples:
#   ./build_and_run.sh                    # バージョン自動抽出してビルド・起動
#   ./build_and_run.sh v1.0.0             # v1.0.0 タグでビルド・起動
#   ./build_and_run.sh v1.0.0 --pull      # Docker Hub から v1.0.0 を pull して起動
#   ./build_and_run.sh v1.0.0 --no-cache  # キャッシュなしでビルド・起動
# ==============================================================================

set -euo pipefail

IMAGE_NAME="amufaamo/singlerepexplorer"
VERSION=""
DO_PULL=false
NO_CACHE=""

# --- 引数パース ---
for arg in "$@"; do
    case "$arg" in
        --pull)
            DO_PULL=true
            ;;
        --no-cache)
            NO_CACHE="--no-cache"
            ;;
        -h|--help)
            sed -n '/^# Usage:/,/^# =====/ { /^# =====/d; s/^# \{0,1\}//; p }' "$0"
            exit 0
            ;;
        -*)
            echo "Unknown option: $arg" >&2
            echo "Run '$0 --help' for usage." >&2
            exit 1
            ;;
        *)
            VERSION="$arg"
            ;;
    esac
done

# --- バージョン解決 ---
if [ -z "$VERSION" ]; then
    echo "No version specified. Trying to extract from app/app.R..."
    VERSION=$(grep 'SingleRepExplorer v' app/app.R \
        | sed -n 's/.*SingleRepExplorer v\([0-9]*\.[0-9]*\.[0-9]*\).*/\1/p' \
        | head -1)

    if [ -z "$VERSION" ]; then
        echo "Warning: Could not extract version automatically."
        read -rp "Please enter the version manually (e.g., v1.0.0): " VERSION
    fi
fi

# v が付いていなければ付ける
[[ "$VERSION" == v* ]] || VERSION="v${VERSION}"

FULL_IMAGE="${IMAGE_NAME}:${VERSION}"

echo "=== SingleRepExplorer Docker Setup Script ==="
echo "  Version : ${VERSION}"
$DO_PULL && echo "  Mode    : pull from Docker Hub" || echo "  Mode    : local build"
echo ""

# --- アーキテクチャ検出（ローカルビルド時のみ使用） ---
ARCH=$(uname -m)
if [[ "$ARCH" == "arm64" ]]; then
    export DOCKERFILE_NAME="Dockerfile.arm64"
    echo "Detected architecture: arm64 → using Dockerfile.arm64"
else
    export DOCKERFILE_NAME="Dockerfile"
    echo "Detected architecture: x86_64 → using Dockerfile"
fi

# --- 既存コンテナを停止・削除 ---
echo ""
echo "Stopping existing containers..."
docker-compose down --remove-orphans 2>/dev/null || true

# --- イメージ削除の確認（ローカルビルド時のみ） ---
if ! $DO_PULL; then
    read -rp "Remove existing local images for this project? (y/N): " -n 1
    echo
    if [[ "$REPLY" =~ ^[Yy]$ ]]; then
        echo "Removing existing images..."
        docker-compose down --rmi all 2>/dev/null || true
    fi
fi

# --- pull または ローカルビルド ---
if $DO_PULL; then
    echo ""
    echo ">>> Pulling ${FULL_IMAGE} from Docker Hub..."
    docker pull "${FULL_IMAGE}"

    # docker-compose が使う IMAGE を上書き
    export APP_IMAGE="${FULL_IMAGE}"

    echo ""
    echo ">>> Starting container with pulled image..."
    # pull モードでは build なしで起動（docker-compose.yml の image: を環境変数で上書き）
    APP_IMAGE="${FULL_IMAGE}" docker run -d \
        --name singlerepexplorer \
        -p 3838:3838 \
        -v "$(pwd)/app:/srv/shiny-server" \
        -e SHINY_LOG_LEVEL=INFO \
        --memory=8g \
        --memory-reservation=4g \
        --restart unless-stopped \
        "${FULL_IMAGE}"
else
    echo ""
    echo ">>> Building ${FULL_IMAGE} locally..."
    # shellcheck disable=SC2086
    docker build ${NO_CACHE} \
        -f "${DOCKERFILE_NAME}" \
        -t "${FULL_IMAGE}" \
        -t "${IMAGE_NAME}:latest" \
        .

    echo ""
    echo ">>> Starting container..."
    docker-compose up -d
fi

echo ""
echo "==================================================================="
echo "  SingleRepExplorer ${VERSION} is starting up!"
echo "  Access URL : http://localhost:3838/"
echo "  View logs  : docker logs -f singlerepexplorer"
echo "             : (or) docker-compose logs -f singlerepexplorer"
echo "  Stop       : docker-compose down"
echo "==================================================================="
