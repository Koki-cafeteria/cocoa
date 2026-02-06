#!/bin/bash
# BOSS DR12 CMASS 共分散行列データのダウンロードスクリプト
#
# BOSS DR12のデータリリースは通常以下の場所にあります：
# - SDSS Data Release 12: https://www.sdss.org/dr12/
# - BOSS DR12 Power Spectrum: https://www.sdss.org/dr12/data_access/
#
# 共分散行列ファイルは通常、パワースペクトルデータと同じディレクトリにあります。
# ファイル名の例：
# - Covariance_Matrix_cmass_ngc_v5.txt
# - Covariance_Matrix_cmass_sgc_v5.txt
# - Covariance_Matrix_cmass_combined_v5.txt

echo "BOSS DR12 CMASS 共分散行列データのダウンロード"
echo "================================================"
echo ""
echo "BOSS DR12のデータリリースから共分散行列を取得してください："
echo ""
echo "1. SDSS Data Release 12 のウェブサイトにアクセス："
echo "   https://www.sdss.org/dr12/"
echo ""
echo "2. BOSS DR12 Power Spectrum データのページに移動"
echo ""
echo "3. 共分散行列ファイルをダウンロード（ファイル名の例）："
echo "   - Covariance_Matrix_cmass_ngc_v5.txt (NGCサンプル用)"
echo "   - Covariance_Matrix_cmass_sgc_v5.txt (SGCサンプル用)"
echo "   - Covariance_Matrix_cmass_combined_v5.txt (結合サンプル用)"
echo ""
echo "4. ダウンロードしたファイルをこのディレクトリに配置："
echo "   $(pwd)"
echo ""
echo "5. 以下のコマンドで変換："
echo "   python convert_boss_dr12_covariance.py \\"
echo "       Covariance_Matrix_cmass_ngc_v5.txt \\"
echo "       boss_dr12_cmass_covariance.txt \\"
echo "       --nk 40 \\"
echo "       --format mono_quad"
echo ""
echo "注意: ファイル名や形式はデータリリースによって異なる場合があります。"
echo "      実際のファイル形式に合わせて --format オプションを調整してください。"

