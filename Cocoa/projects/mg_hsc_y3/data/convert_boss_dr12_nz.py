#!/usr/bin/env python3
"""
BOSS DR12のn(z)ファイルをCocoa形式に変換するスクリプト

BOSS DR12形式:
  z <nobs> <wc nobs> <wc wfkp nobs>

Cocoa形式（単一zビンの場合）:
  z n(z)

使用方法:
    python convert_boss_dr12_nz.py input_file.txt output_file.nz
"""

import numpy as np
import sys
import os

def read_boss_dr12_nz(filename):
    """
    BOSS DR12のn(z)ファイルを読み込む
    
    Returns:
        z: 赤方偏移の配列
        nz: n(z)の配列（通常は<nobs>列を使用、または<wc wfkp nobs>を使用）
    """
    z = []
    nz = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            # コメント行をスキップ
            if line.startswith('#') or len(line) == 0:
                continue
            
            # データ行をパース
            parts = line.split()
            if len(parts) >= 4:
                z.append(float(parts[0]))  # z
                # <wc wfkp nobs>を使用（重み付き密度、通常はこれを使用）
                nz.append(float(parts[3]))  # <wc wfkp nobs>
    
    return np.array(z), np.array(nz)


def write_cocoa_nz(z, nz, output_file, n_bins=1):
    """
    Cocoa形式でn(z)ファイルを書き込む
    
    Parameters:
        z: 赤方偏移の配列
        nz: n(z)の配列（単一ビンの場合）
        output_file: 出力ファイル名
        n_bins: トモグラフィービンの数（BOSS DR12は1）
    """
    print(f"Writing Cocoa format n(z) file: {output_file}")
    print(f"  - Number of z bins: {n_bins}")
    print(f"  - Number of z points: {len(z)}")
    print(f"  - z range: {z[0]:.6f} - {z[-1]:.6f}")
    
    with open(output_file, 'w') as f:
        for i in range(len(z)):
            # Cocoa形式: z n(z)_bin0 n(z)_bin1 ... (単一ビンの場合は2列)
            if n_bins == 1:
                f.write(f"{z[i]:.6f}  {nz[i]:.12e}\n")
            else:
                # 複数ビンの場合は、各ビンに対応するn(z)を並べる
                # ここでは単一ビンのみ対応
                f.write(f"{z[i]:.6f}  {nz[i]:.12e}\n")
    
    print(f"Successfully wrote {len(z)} lines.")


def main():
    if len(sys.argv) < 3:
        print("Usage: python convert_boss_dr12_nz.py <input_file> <output_file>")
        print("\nExample:")
        print("  python convert_boss_dr12_nz.py \\")
        print("      'boss dr12/Density_galaxies_cmass_ngc_v5.txt' \\")
        print("      boss_dr12_cmass.nz")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
    
    print(f"Reading BOSS DR12 n(z) file from: {input_file}")
    z, nz = read_boss_dr12_nz(input_file)
    
    print(f"\nLoaded {len(z)} z points")
    print(f"  z range: {z[0]:.6f} - {z[-1]:.6f}")
    print(f"  n(z) range: {nz.min():.6e} - {nz.max():.6e}")
    print(f"  Total n(z) integral: {np.trapz(nz, z):.6e}")
    
    write_cocoa_nz(z, nz, output_file, n_bins=1)
    
    print("\n" + "="*60)
    print("Conversion completed successfully!")
    print("="*60)
    print(f"Output file: {output_file}")
    print(f"\nNext step: Update .dataset file to use this n(z) file")


if __name__ == "__main__":
    main()

