#!/usr/bin/env python3
"""
BOSS DR12のパワースペクトルデータをKaiser RSD用のデータファイル形式に変換するスクリプト

使用方法:
    python convert_boss_dr12_data.py input_file.txt output_prefix

入力ファイル形式:
    BOSS DR12のPower_Spectrum_cmass_ngc_v5形式
    カラム: k-centerbin, k-eff, Monopole-Pshotnoise, Quadrupole, Hexadecapole, number of modes, Pshotnoise
"""

import numpy as np
import sys
import os

def read_boss_dr12_file(filename):
    """
    BOSS DR12のパワースペクトルファイルを読み込む
    
    Returns:
        k_eff: 実効k値 [h/Mpc]
        P0: モノポール（ショットノイズを引いた後）[(Mpc/h)^3]
        P2: クアドルポール [(Mpc/h)^3]
        Pshotnoise: ショットノイズ [(Mpc/h)^3]
        n_modes: モード数
    """
    k_eff = []
    P0 = []
    P2 = []
    Pshotnoise = []
    n_modes = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            # コメント行をスキップ
            if line.startswith('#') or len(line) == 0:
                continue
            
            # データ行をパース
            parts = line.split()
            if len(parts) >= 7:
                k_eff.append(float(parts[1]))  # k-eff
                P0.append(float(parts[2]))     # Monopole-Pshotnoise
                P2.append(float(parts[3]))     # Quadrupole
                n_modes.append(int(parts[5]))  # number of modes
                Pshotnoise.append(float(parts[6]))  # Pshotnoise
    
    return np.array(k_eff), np.array(P0), np.array(P2), np.array(Pshotnoise), np.array(n_modes)


def create_data_files(k_eff, P0, P2, output_prefix, z_bin_index=0):
    """
    Kaiser RSD用のデータファイルを作成
    
    Parameters:
        k_eff: k値の配列 [h/Mpc]
        P0: モノポールの配列 [(Mpc/h)^3]
        P2: クアドルポールの配列 [(Mpc/h)^3]
        output_prefix: 出力ファイルのプレフィックス
        z_bin_index: 赤方偏移ビンのインデックス（単一zビンの場合は0）
    """
    Nk = len(k_eff)
    Nz = 1  # BOSS DR12は単一の赤方偏移ビン
    
    # データベクトルのサイズ: 2 * Nk * Nz (モノポール + クアドルポール)
    N_data = 2 * Nk * Nz
    
    print(f"Creating data files for BOSS DR12:")
    print(f"  - Number of k bins: {Nk}")
    print(f"  - Number of z bins: {Nz}")
    print(f"  - Total data points: {N_data}")
    
    # --- 1. データベクトルファイルの作成 ---
    # 構造: [P_0(k1,z1), ..., P_0(kNk,z1), P_2(k1,z1), ..., P_2(kNk,z1)]
    data_filename = f"{output_prefix}_datavector.txt"
    print(f"\nGenerating {data_filename}...")
    
    indices = np.arange(N_data)
    values = np.zeros(N_data)
    
    # モノポール部分（最初のNk個）
    values[0:Nk] = P0
    
    # クアドルポール部分（次のNk個）
    values[Nk:2*Nk] = P2
    
    data_vector = np.column_stack((indices, values))
    np.savetxt(data_filename, data_vector, fmt='%d  %.6e')
    print(f"  - Monopole range: indices 0-{Nk-1}")
    print(f"  - Quadrupole range: indices {Nk}-{2*Nk-1}")
    
    # --- 2. マスクファイルの作成 ---
    # 全てのデータ点を使用（値は1）
    mask_filename = f"{output_prefix}_mask.txt"
    print(f"\nGenerating {mask_filename}...")
    mask_values = np.ones(N_data, dtype=int)
    mask_file_content = np.column_stack((indices, mask_values))
    np.savetxt(mask_filename, mask_file_content, fmt='%d  %d')
    
    # --- 3. k値の情報を保存（参考用） ---
    k_info_filename = f"{output_prefix}_k_bins.txt"
    print(f"\nGenerating {k_info_filename} (for reference)...")
    k_info = np.column_stack((np.arange(Nk), k_eff))
    np.savetxt(k_info_filename, k_info, fmt='%d  %.6e', 
               header='index  k_eff[h/Mpc]')
    
    # --- 4. 共分散行列ファイルについて ---
    # 注意: BOSS DR12のデータには共分散行列の情報が含まれていないため、
    # ユーザーが別途準備する必要があります
    print(f"\nNote: Covariance matrix file ({output_prefix}_covariance.txt) must be created separately.")
    print(f"      The covariance matrix should be {N_data}x{N_data} in size.")
    print(f"      Format: [row_index] [col_index] [covariance_value]")
    
    print("\n" + "="*60)
    print("Files created successfully!")
    print("="*60)
    print(f"  - Data vector: {data_filename}")
    print(f"  - Mask: {mask_filename}")
    print(f"  - k bins info: {k_info_filename}")
    print(f"\nNext steps:")
    print(f"  1. Create covariance matrix file: {output_prefix}_covariance.txt")
    print(f"  2. Update .dataset file to set:")
    print(f"     - Nk = {Nk}")
    print(f"     - Nz = {Nz}")
    print(f"     - k_min = {k_eff[0]:.6f}")
    print(f"     - k_max = {k_eff[-1]:.6f}")
    print(f"     - z_min, z_max: Set according to BOSS DR12 redshift range")


def main():
    if len(sys.argv) < 3:
        print("Usage: python convert_boss_dr12_data.py <input_file> <output_prefix>")
        print("\nExample:")
        print("  python convert_boss_dr12_data.py Power_Spectrum_cmass_ngc_v5.txt boss_dr12_cmass")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_prefix = sys.argv[2]
    
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
    
    print(f"Reading BOSS DR12 data from: {input_file}")
    k_eff, P0, P2, Pshotnoise, n_modes = read_boss_dr12_file(input_file)
    
    print(f"\nLoaded {len(k_eff)} k bins")
    print(f"  k range: {k_eff[0]:.6f} - {k_eff[-1]:.6f} [h/Mpc]")
    print(f"  P0 range: {P0.min():.2e} - {P0.max():.2e} [(Mpc/h)^3]")
    print(f"  P2 range: {P2.min():.2e} - {P2.max():.2e} [(Mpc/h)^3]")
    
    create_data_files(k_eff, P0, P2, output_prefix)


if __name__ == "__main__":
    main()

