#!/usr/bin/env python3
"""
モックデータから共分散行列を作成するコード（修正版）

修正点:
1. データベクトルを塊形式 [P0(k1), ..., P0(kNk), P2(k1), ..., P2(kNk)] に変更
2. 共分散行列を直接Cocoa形式で出力
"""

import numpy as np
import tarfile
import os
import glob
import sys

print(f"Numpy version: {np.__version__}")

# --- 設定項目 ---
TAR_FILE_PATH = "./boss dr12/Power_Spectrum_cmass_ngc_v5_Patchy.tar.gz"
EXTRACT_DIR = "./patchy_mocks_cmass_ngc"
K_MAX_FIT = 0.395040  # 全k範囲を使用（40ビン、データベクトルと一致させる）
OUTPUT_COV_FILE = "kaiser_rsd_covariance.txt"
OUTPUT_MEAN_DATA_FILE = "kaiser_rsd_datavector_mean.txt"

print(f"'{TAR_FILE_PATH}' を '{EXTRACT_DIR}' に展開します...")

try:
    with tarfile.open(TAR_FILE_PATH, 'r:gz') as tar:
        tar.extractall(path=EXTRACT_DIR)
    print("展開完了。")
except FileNotFoundError:
    print(f"エラー: '{TAR_FILE_PATH}' が見つかりません。", file=sys.stderr)
except Exception as e:
    print(f"展開エラー: {e}", file=sys.stderr)

print("モックファイルの読み込みを開始...")

file_pattern = os.path.join(EXTRACT_DIR, "Power_Spectrum_cmass_ngc_v5_Patchy_*.txt")
mock_files = sorted(glob.glob(file_pattern))

if not mock_files:
    print(f"エラー: '{EXTRACT_DIR}' 内に '{file_pattern}' に一致するファイルが見つかりません。", file=sys.stderr)
    print("-> `EXTRACT_DIR` の中身を確認し、ファイルパスが正しいか確認してください。")
else:
    print(f"{len(mock_files)} 個のモックファイルを発見。")
    all_mock_vectors = []
    k_bins_filter_indices = None
    k_bins_count = 0
    
    for i, filepath in enumerate(mock_files):
        try:
            # README Section 3.1 形式:
            # 0:k-center, 1:k-eff, 2:P0, 3:P2, 4:P4, 5:N_modes, 6:P0_shotnoise
            # 使うのは 1:k-eff, 2:P0, 3:P2
            data = np.loadtxt(filepath, comments='#', usecols=(1, 2, 3))
        except Exception as e:
            print(f"読み込みエラー {filepath}: {e}", file=sys.stderr)
            continue
        
        # 最初のファイルで k<=K_MAX_FIT のインデックスを取得
        if i == 0:
            k_eff_values = data[:, 0]
            k_bins_filter_indices = np.where(k_eff_values <= K_MAX_FIT)[0]
            k_bins_count = len(k_bins_filter_indices)
            
            if k_bins_count == 0:
                print(f"エラー: k <= {K_MAX_FIT} のkビンがありません。", file=sys.stderr)
                break
            
            print(f"k <= {K_MAX_FIT} h/Mpc の {k_bins_count} ビンを使用します。")
            print(f"(データベクトルの全長は 2 * {k_bins_count} = {2 * k_bins_count} になります)")
        
        # k_max でフィルタリング
        p0 = data[k_bins_filter_indices, 1]
        p2 = data[k_bins_filter_indices, 2]
        
        # 【修正】Cocoa 形式 [P0(k1), ..., P0(kNk), P2(k1), ..., P2(kNk)] に並べ替え
        # まずすべてのモノポール、その後にすべてのクアドルポール
        block_vector = np.concatenate([p0, p2])
        
        all_mock_vectors.append(block_vector)
    
    # (N_real, N_bins) の形状を持つ Numpy 配列を作成
    if all_mock_vectors:
        mock_matrix = np.array(all_mock_vectors)
        print(f"モック行列の形状 (N_real, N_bins): {mock_matrix.shape}")
    else:
        print("データベクトルが作成されませんでした。")

if 'mock_matrix' in locals() and mock_matrix.size > 0:
    print("共分散行列を計算中...")
    
    # np.cov はデフォルトで (N-1) で割るため、Eq. 7 と一致します。
    # rowvar=False は、入力の形状が (N_real, N_bins) であることを示し、
    # (N_bins, N_bins) の共分散行列を返します。
    cov_matrix = np.cov(mock_matrix, rowvar=False)
    
    print(f"共分散行列の形状 (N_bins, N_bins): {cov_matrix.shape}")
    
    # 参考: モックの平均データベクトル (Eq. 8 の D_tilde)
    mean_data_vector = np.mean(mock_matrix, axis=0)
    print(f"平均データベクトルの形状 (N_bins,): {mean_data_vector.shape}")
    
    # Cocoa形式で共分散行列を保存
    print(f"\nCocoa形式で共分散行列を '{OUTPUT_COV_FILE}' に保存中...")
    n = cov_matrix.shape[0]
    with open(OUTPUT_COV_FILE, 'w') as f:
        for i in range(n):
            for j in range(n):
                f.write(f"{i}  {j}  {cov_matrix[i, j]:.6e}\n")
    print(f"共分散行列を保存しました ({n}x{n} = {n*n} エントリ)")
    
    # 平均データベクトルも保存（参考用）
    print(f"\n平均データベクトルを '{OUTPUT_MEAN_DATA_FILE}' に保存中...")
    with open(OUTPUT_MEAN_DATA_FILE, 'w') as f:
        for idx, value in enumerate(mean_data_vector):
            f.write(f"{idx}  {value:.6e}\n")
    print(f"平均データベクトルを保存しました ({len(mean_data_vector)} エントリ)")
    
    print("\n" + "="*60)
    print("処理完了!")
    print("="*60)
    print(f"  - 共分散行列: {OUTPUT_COV_FILE}")
    print(f"  - 平均データベクトル: {OUTPUT_MEAN_DATA_FILE}")
    print(f"  - kビン数: {k_bins_count}")
    print(f"  - データベクトルサイズ: {2 * k_bins_count}")
    
else:
    print("モック行列が作成されていないため、計算をスキップします。")

