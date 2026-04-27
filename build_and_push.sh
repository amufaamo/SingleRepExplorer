#!/bin/bash
# ==============================================================================
# SingleRepExplorer Build & Push Script
#
# Usage:
#   ./build_and_push.sh [VERSION] [OPTIONS]
#
# Arguments:
#   VERSION      バージョン番号 (例: v1.0.0)。省略時は app/app.R から自動抽出。
#
# Options:
#   --no-latest  latest タグをビルド・プッシュしない
#   --no-push    ビルドのみ行い、Docker Hub へのプッシュをスキップする
#   -h, --help   このヘルプを表示する
#
# Examples:
#   ./build_and_push.sh v1.0.0
#   ./build_and_push.sh               # app.R からバージョン自動抽出
#   ./build_and_push.sh v1.0.0 --no-latest
#   ./build_and_push.sh v1.0.0 --no-push
# ==============================================================================

set -euo pipefail

IMAGE_NAME="amufaamo/singlerepexplorer"
PUSH_LATEST=true
DO_PUSH=true
VERSION=""

# --- 引数パース ---
for arg in "$@"; do
    case "$arg" in
        --no-latest)
            PUSH_LATEST=false
            ;;
        --no-push)
            DO_PUSH=false
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
    # Match either tags$title("SingleRepExplorer vX.Y.Z ...") or
    # title = "SingleRepExplorer vX.Y.Z ..." styles.
    VERSION=$(grep -E 'SingleRepExplorer v[0-9]+\.[0-9]+\.[0-9]+' app/app.R \
        | sed -nE 's/.*SingleRepExplorer (v[0-9]+\.[0-9]+\.[0-9]+).*/\1/p' \
        | head -1)

    if [ -z "$VERSION" ]; then
        echo "Warning: Could not extract version automatically."
        read -rp "Please enter the version manually (e.g., v1.0.0): " VERSION
    else
        echo "Detected version: $VERSION"
    fi
fi

# vが付いていなければ付ける
[[ "$VERSION" == v* ]] || VERSION="v${VERSION}"

# --- 確認プロンプト ---
echo ""
echo "=============================================="
echo "  Image : ${IMAGE_NAME}:${VERSION}"
$PUSH_LATEST && echo "  Also  : ${IMAGE_NAME}:latest"
$DO_PUSH     && echo "  Action: build + push" || echo "  Action: build only (no push)"
echo "=============================================="
read -rp "Proceed? (y/N): " confirm
if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
    echo "Aborted."
    exit 0
fi

# --- ビルド ---
BUILD_TAGS="-t ${IMAGE_NAME}:${VERSION}"
$PUSH_LATEST && BUILD_TAGS="$BUILD_TAGS -t ${IMAGE_NAME}:latest"

echo ""
echo ">>> Building ${IMAGE_NAME}:${VERSION} ..."
# shellcheck disable=SC2086
docker build $BUILD_TAGS .

echo "Build successful!"

# --- プッシュ ---
if $DO_PUSH; then
    echo ""
    echo ">>> Pushing ${IMAGE_NAME}:${VERSION} ..."
    docker push "${IMAGE_NAME}:${VERSION}"

    if $PUSH_LATEST; then
        echo ">>> Pushing ${IMAGE_NAME}:latest ..."
        docker push "${IMAGE_NAME}:latest"
    fi

    echo ""
    echo "Done! Image pushed to Docker Hub."
else
    echo ""
    echo "Done! (push skipped)"
fi
