#!/usr/bin/env python3
"""
BOSS DR12のデータをCOCOAの対数等間隔kビンに補間して変換するスクリプト

使用方法:
    python convert_boss_dr12_to_cocoa_kbins.py

入力:
    - Power_Spectrum_cmass_ngc_v5.txt: BOSS DR12の元データ
    - boss_dr12_cmass_covariance.txt: 既存の共分散行列（モックデータから計算済み）

出力:
    - boss_dr12_cmass_datavector_cocoa_kbins.txt: 補間後のデータベクトル
    - boss_dr12_cmass_covariance_cocoa_kbins.txt: 補間後の共分散行列（モックデータから再計算）
    - boss_dr12_cmass_mask_cocoa_kbins.txt: マスクファイル
"""

import numpy as np
from scipy.interpolate import interp1d
import sys
import os

def read_boss_dr12_file(filename):
    """BOSS DR12のパワースペクトルファイルを読み込む"""
    k_eff = []
    P0 = []
    P2 = []
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#') or len(line.strip()) == 0:
                continue
            parts = line.split()
            if len(parts) >= 4:
                k_eff.append(float(parts[1]))  # k-eff
                P0.append(float(parts[2]))     # Monopole-Pshotnoise
                P2.append(float(parts[3]))     # Quadrupole
    
    return np.array(k_eff), np.array(P0), np.array(P2)


def calculate_cocoa_k_bins(k_min, k_max, Nk):
    """COCOAの対数等間隔kビンを計算"""
    log_dk = (np.log10(k_max) - np.log10(k_min)) / (Nk - 1)
    k_bins = np.array([10**(np.log10(k_min) + i * log_dk) for i in range(Nk)])
    return k_bins


def interpolate_power_spectrum(k_original, P_original, k_target, kind='cubic'):
    """
    パワースペクトルを補間
    
    Parameters:
        k_original: 元のk値
        P_original: 元のパワースペクトル値
        k_target: 補間先のk値
        kind: 補間方法 ('linear', 'cubic', etc.)
    
    Returns:
        P_interpolated: 補間後のパワースペクトル値
    """
    # 対数スケールで補間（パワースペクトルは通常対数スケールで補間する）
    # ただし、負の値がある場合は線形スケールで補間
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


def main():
    # 設定
    input_file = "Power_Spectrum_cmass_ngc_v5.txt"
    output_prefix = "boss_dr12_cmass_cocoa_kbins"
    
    # COCOAのkビン設定（.datasetファイルから）
    k_min = 0.007424
    k_max = 0.395040
    Nk = 40
    Nz = 1
    
    print("="*60)
    print("BOSS DR12データをCOCOAのkビンに補間")
    print("="*60)
    
    # 1. BOSS DR12の元データを読み込み
    print(f"\n1. Reading BOSS DR12 data from: {input_file}")
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
    
    k_eff, P0_eff, P2_eff = read_boss_dr12_file(input_file)
    print(f"   Loaded {len(k_eff)} k bins")
    print(f"   k range: {k_eff[0]:.6f} - {k_eff[-1]:.6f} [h/Mpc]")
    print(f"   P0 range: {P0_eff.min():.2e} - {P0_eff.max():.2e} [(Mpc/h)^3]")
    print(f"   P2 range: {P2_eff.min():.2e} - {P2_eff.max():.2e} [(Mpc/h)^3]")
    
    # 2. COCOAの対数等間隔kビンを計算
    print(f"\n2. Calculating COCOA k bins (logarithmically spaced)")
    k_cocoa = calculate_cocoa_k_bins(k_min, k_max, Nk)
    print(f"   COCOA k bins: {len(k_cocoa)} bins")
    print(f"   k range: {k_cocoa[0]:.6f} - {k_cocoa[-1]:.6f} [h/Mpc]")
    
    # 3. データを補間
    print(f"\n3. Interpolating P0 and P2 to COCOA k bins")
    P0_interp = interpolate_power_spectrum(k_eff, P0_eff, k_cocoa, kind='cubic')
    P2_interp = interpolate_power_spectrum(k_eff, P2_eff, k_cocoa, kind='cubic')
    
    # 補間結果の確認
    valid_mask_p0 = ~np.isnan(P0_interp)
    valid_mask_p2 = ~np.isnan(P2_interp)
    print(f"   Valid P0 bins: {np.sum(valid_mask_p0)}/{len(P0_interp)}")
    print(f"   Valid P2 bins: {np.sum(valid_mask_p2)}/{len(P2_interp)}")
    
    if np.any(np.isnan(P0_interp)) or np.any(np.isnan(P2_interp)):
        print("   Warning: Some bins could not be interpolated (NaN values)")
        print("   These will be masked in the mask file")
    
    # 4. データベクトルファイルを作成
    print(f"\n4. Creating data vector file")
    N_data = 2 * Nk * Nz
    indices = np.arange(N_data)
    values = np.zeros(N_data)
    
    # モノポール部分（最初のNk個）
    values[0:Nk] = P0_interp
    # クアドルポール部分（次のNk個）
    values[Nk:2*Nk] = P2_interp
    
    data_vector_file = f"{output_prefix}_datavector.txt"
    data_vector = np.column_stack((indices, values))
    np.savetxt(data_vector_file, data_vector, fmt='%d  %.6e')
    print(f"   Saved: {data_vector_file}")
    
    # 5. マスクファイルを作成（NaNの場合は0、それ以外は1）
    print(f"\n5. Creating mask file")
    mask_values = np.ones(N_data, dtype=int)
    mask_values[0:Nk] = valid_mask_p0.astype(int)
    mask_values[Nk:2*Nk] = valid_mask_p2.astype(int)
    
    mask_file = f"{output_prefix}_mask.txt"
    mask_data = np.column_stack((indices, mask_values))
    np.savetxt(mask_file, mask_data, fmt='%d  %d')
    print(f"   Saved: {mask_file}")
    print(f"   Masked bins: {np.sum(mask_values == 0)}/{N_data}")
    
    # 6. kビン情報を保存（参考用）
    k_info_file = f"{output_prefix}_k_bins.txt"
    k_info = np.column_stack((np.arange(Nk), k_cocoa))
    np.savetxt(k_info_file, k_info, fmt='%d  %.6e', 
               header='index  k_cocoa[h/Mpc]')
    print(f"\n6. Saved k bins info: {k_info_file}")
    
    print("\n" + "="*60)
    print("Data conversion completed!")
    print("="*60)
    print(f"\nNext steps:")
    print(f"  1. Recalculate covariance matrix from mock data")
    print(f"     (Use fix_mock_covariance_code.py with COCOA k bins)")
    print(f"  2. Update .dataset file to use:")
    print(f"     - data_file = {data_vector_file}")
    print(f"     - cov_file = {output_prefix}_covariance.txt")
    print(f"     - mask_file = {mask_file}")
    print(f"     - k_min = {k_min:.6f}")
    print(f"     - k_max = {k_max:.6f}")
    print(f"     - Nk = {Nk}")


if __name__ == "__main__":
    main()














