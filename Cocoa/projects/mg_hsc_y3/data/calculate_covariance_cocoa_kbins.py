#!/usr/bin/env python3
"""
モックデータをCOCOAのkビンに補間してから共分散行列を計算するスクリプト

使用方法:
    python calculate_covariance_cocoa_kbins.py

入力:
    - patchy_mocks_cmass_ngc/Power_Spectrum_cmass_ngc_v5_Patchy_*.txt: モックデータファイル

出力:
    - boss_dr12_cmass_cocoa_kbins_covariance.txt: 補間後の共分散行列
    - boss_dr12_cmass_cocoa_kbins_datavector_mean.txt: 平均データベクトル
"""

import numpy as np
from scipy.interpolate import interp1d
import tarfile
import os
import glob
import sys

# 設定
TAR_FILE_PATH = "./boss dr12/Power_Spectrum_cmass_ngc_v5_Patchy.tar.gz"
EXTRACT_DIR = "./patchy_mocks_cmass_ngc"
K_MAX_FIT = 0.395040  # BOSS DR12の最大k値
OUTPUT_COV_FILE = "boss_dr12_cmass_cocoa_kbins_covariance.txt"
OUTPUT_MEAN_DATA_FILE = "boss_dr12_cmass_cocoa_kbins_datavector_mean.txt"

# COCOAのkビン設定
k_min = 0.007424
k_max = 0.395040
Nk = 40

def calculate_cocoa_k_bins(k_min, k_max, Nk):
    """COCOAの対数等間隔kビンを計算"""
    log_dk = (np.log10(k_max) - np.log10(k_min)) / (Nk - 1)
    k_bins = np.array([10**(np.log10(k_min) + i * log_dk) for i in range(Nk)])
    return k_bins

def interpolate_power_spectrum(k_original, P_original, k_target, kind='cubic'):
    """パワースペクトルを補間"""
    if np.any(P_original < 0):
        # 負の値がある場合は線形スケールで補間
        interp_func = interp1d(k_original, P_original, kind=kind,
                              bounds_error=False, fill_value=np.nan)
        P_interp = interp_func(k_target)
    else:
        # 正の値のみの場合は対数スケールで補間
        log_k_orig = np.log10(k_original)
        log_P_orig = np.log10(P_original)
        interp_func = interp1d(log_k_orig, log_P_orig, kind=kind,
                              bounds_error=False, fill_value=np.nan)
        log_P_interp = interp_func(np.log10(k_target))
        P_interp = 10**log_P_interp
    return P_interp

print("="*60)
print("モックデータからCOCOA kビン用の共分散行列を計算")
print("="*60)

# COCOAのkビンを計算
k_cocoa = calculate_cocoa_k_bins(k_min, k_max, Nk)
print(f"\nCOCOA k bins: {len(k_cocoa)} bins")
print(f"k range: {k_cocoa[0]:.6f} - {k_cocoa[-1]:.6f} [h/Mpc]")

# モックファイルを読み込む
file_pattern = os.path.join(EXTRACT_DIR, "Power_Spectrum_cmass_ngc_v5_Patchy_*.txt")
mock_files = sorted(glob.glob(file_pattern))

if not mock_files:
    print(f"\nError: '{EXTRACT_DIR}' 内に '{file_pattern}' に一致するファイルが見つかりません。")
    print("-> モックデータが展開されているか確認してください。")
    sys.exit(1)

print(f"\n{len(mock_files)} 個のモックファイルを発見。")

all_mock_vectors = []

for i, filepath in enumerate(mock_files):
    try:
        # BOSS DR12形式: k-centerbin, k-eff, P0, P2, P4, N_modes, P0_shotnoise
        data = np.loadtxt(filepath, comments='#', usecols=(1, 2, 3))  # k-eff, P0, P2
        
        k_eff = data[:, 0]
        p0 = data[:, 1]
        p2 = data[:, 2]
        
        # k_maxでフィルタリング
        mask = k_eff <= K_MAX_FIT
        k_eff_filtered = k_eff[mask]
        p0_filtered = p0[mask]
        p2_filtered = p2[mask]
        
        # COCOAのkビンに補間
        p0_interp = interpolate_power_spectrum(k_eff_filtered, p0_filtered, k_cocoa, kind='cubic')
        p2_interp = interpolate_power_spectrum(k_eff_filtered, p2_filtered, k_cocoa, kind='cubic')
        
        # NaNチェック
        if np.any(np.isnan(p0_interp)) or np.any(np.isnan(p2_interp)):
            print(f"Warning: NaN values in mock {i+1}, skipping...")
            continue
        
        # Cocoa形式 [P0(k1), ..., P0(kNk), P2(k1), ..., P2(kNk)] に並べ替え
        cocoa_vector = np.empty(2 * Nk)
        cocoa_vector[:Nk] = p0_interp
        cocoa_vector[Nk:] = p2_interp
        
        all_mock_vectors.append(cocoa_vector)
        
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        continue

if not all_mock_vectors:
    print("Error: 有効なモックデータがありません。")
    sys.exit(1)

# モック行列を作成
mock_matrix = np.array(all_mock_vectors)
print(f"\nモック行列の形状 (N_real, N_bins): {mock_matrix.shape}")

# 共分散行列を計算
print("\n共分散行列を計算中...")
cov_matrix = np.cov(mock_matrix, rowvar=False)
print(f"共分散行列の形状 (N_bins, N_bins): {cov_matrix.shape}")

# 平均データベクトル
mean_data_vector = np.mean(mock_matrix, axis=0)
print(f"平均データベクトルの形状 (N_bins,): {mean_data_vector.shape}")

# 共分散行列をCOCOA形式で保存
print(f"\n共分散行列を保存: {OUTPUT_COV_FILE}")
with open(OUTPUT_COV_FILE, 'w') as f:
    for i in range(cov_matrix.shape[0]):
        for j in range(cov_matrix.shape[1]):
            f.write(f"{i}  {j}  {cov_matrix[i, j]:.6e}\n")

# 平均データベクトルを保存
print(f"平均データベクトルを保存: {OUTPUT_MEAN_DATA_FILE}")
with open(OUTPUT_MEAN_DATA_FILE, 'w') as f:
    for i, value in enumerate(mean_data_vector):
        f.write(f"{i}  {value:.6e}\n")

print("\n" + "="*60)
print("完了!")
print("="*60)
print(f"\n生成されたファイル:")
print(f"  - {OUTPUT_COV_FILE}")
print(f"  - {OUTPUT_MEAN_DATA_FILE}")














