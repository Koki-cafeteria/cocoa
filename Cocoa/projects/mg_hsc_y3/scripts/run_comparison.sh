#!/bin/bash
# Kaiser RSD vs Velocileptors 比較を自動実行するスクリプト

set -e  # エラーが発生したら停止

# カラー設定
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Kaiser RSD vs Velocileptors 比較実行${NC}"
echo -e "${BLUE}========================================${NC}"

# Cocoaディレクトリに移動
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COCOA_DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
cd "${COCOA_DIR}"

echo -e "\n${BLUE}作業ディレクトリ: ${COCOA_DIR}${NC}"

# Cocoaの環境を読み込む
if [ ! -f start_cocoa.sh ]; then
    echo -e "${RED}❌ start_cocoa.sh が見つかりません${NC}"
    exit 1
fi

echo -e "${YELLOW}📦 Cocoa環境を読み込んでいます...${NC}"
source start_cocoa.sh

# 出力ファイルのパス
KAISER_OUTPUT="projects/mg_hsc_y3/data/mg_hsc_y3_theory_kaiser_rsd_boss_dr12.modelvector"
VELO_OUTPUT="projects/mg_hsc_y3/data/velocileptors_pk_multipoles_mu10.txt"

# Kaiser RSDの実行確認
NEED_KAISER=false
if [ ! -f "${KAISER_OUTPUT}" ]; then
    echo -e "${YELLOW}⚠️  Kaiser RSD の出力ファイルが見つかりません${NC}"
    NEED_KAISER=true
else
    echo -e "${GREEN}✅ Kaiser RSD の出力ファイルが見つかりました${NC}"
    echo -e "   ${KAISER_OUTPUT}"
    read -p "Kaiser RSDを再実行しますか？ (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        NEED_KAISER=true
    fi
fi

# Velocileptorsの実行確認
NEED_VELO=false
if [ ! -f "${VELO_OUTPUT}" ]; then
    echo -e "${YELLOW}⚠️  Velocileptors の出力ファイルが見つかりません${NC}"
    NEED_VELO=true
else
    echo -e "${GREEN}✅ Velocileptors の出力ファイルが見つかりました${NC}"
    echo -e "   ${VELO_OUTPUT}"
    read -p "Velocileptorsを再実行しますか？ (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        NEED_VELO=true
    fi
fi

# Kaiser RSDの実行
if [ "${NEED_KAISER}" = true ]; then
    echo -e "\n${BLUE}========================================${NC}"
    echo -e "${BLUE}Kaiser RSD を実行中...${NC}"
    echo -e "${BLUE}========================================${NC}"
    cobaya-run projects/mg_hsc_y3/EXAMPLE_EVALUATE_kaiser_RSD.yaml
    echo -e "${GREEN}✅ Kaiser RSD 完了${NC}"
else
    echo -e "\n${GREEN}ℹ️  Kaiser RSD はスキップします${NC}"
fi

# Velocileptorsの実行
if [ "${NEED_VELO}" = true ]; then
    echo -e "\n${BLUE}========================================${NC}"
    echo -e "${BLUE}Velocileptors を実行中...${NC}"
    echo -e "${BLUE}========================================${NC}"
    cobaya-run projects/mg_hsc_y3/EXAMPLE_EVALUATE_velocileptors.yaml
    echo -e "${GREEN}✅ Velocileptors 完了${NC}"
else
    echo -e "\n${GREEN}ℹ️  Velocileptors はスキップします${NC}"
fi

# 比較プロットの作成
echo -e "\n${BLUE}========================================${NC}"
echo -e "${BLUE}比較プロットを作成中...${NC}"
echo -e "${BLUE}========================================${NC}"

python projects/mg_hsc_y3/scripts/compare_kaiser_velocileptors.py

echo -e "\n${GREEN}========================================${NC}"
echo -e "${GREEN}✅ すべての処理が完了しました！${NC}"
echo -e "${GREEN}========================================${NC}"

echo -e "\n${BLUE}📊 出力ファイル:${NC}"
echo -e "   • Kaiser RSD: ${KAISER_OUTPUT}"
echo -e "   • Velocileptors: ${VELO_OUTPUT}"
echo -e "   • 比較プロット1: projects/mg_hsc_y3/data/comparison_kaiser_velocileptors_single_z.png"
echo -e "   • 比較プロット2: projects/mg_hsc_y3/data/comparison_kaiser_velocileptors_multiple_z.png"

echo -e "\n${BLUE}📖 詳細はガイドを参照:${NC}"
echo -e "   projects/mg_hsc_y3/KAISER_VELOCILEPTORS_COMPARISON_GUIDE.md"



