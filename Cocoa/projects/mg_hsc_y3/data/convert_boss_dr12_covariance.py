#!/usr/bin/env python3
"""
BOSS DR12の共分散行列をKaiser RSD用の形式に変換するスクリプト

BOSS DR12のデータリリースでは、共分散行列は通常以下の形式で提供されます：
1. モノポールとクアドルポールの共分散行列が統合された形式
2. モノポール、クアドルポール、ヘキサデカポールを含む形式

使用方法:
    python convert_boss_dr12_covariance.py <input_cov_file> <output_file> [--format FORMAT]
    
    FORMATオプション:
        - 'full': モノポール+クアドルポール+ヘキサデカポールを含む完全な共分散行列
        - 'mono_quad': モノポールとクアドルポールのみを含む共分散行列
        - 'separate': モノポールとクアドルポールが別々のファイルとして提供される場合
"""

import numpy as np
import sys
import os
import argparse

def read_covariance_matrix_full(filename, nk):
    """
    モノポール+クアドルポール+ヘキサデカポールを含む完全な共分散行列を読み込む
    
    BOSS DR12の形式: 3*nk x 3*nk の行列
    順序: [P0(k1), ..., P0(nk), P2(k1), ..., P2(nk), P4(k1), ..., P4(nk)]
    
    Returns:
        cov_matrix: 共分散行列 (3*nk x 3*nk)
    """
    # ファイルを読み込む（コメント行をスキップ）
    data = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or len(line) == 0:
                continue
            # 数値行をパース
            try:
                row = [float(x) for x in line.split()]
                data.append(row)
            except ValueError:
                continue
    
    cov_matrix = np.array(data)
    
    # サイズチェック
    expected_size = 3 * nk
    if cov_matrix.shape[0] != expected_size or cov_matrix.shape[1] != expected_size:
        print(f"Warning: Expected {expected_size}x{expected_size} matrix, got {cov_matrix.shape}")
        print(f"Assuming the matrix contains only monopole and quadrupole (2*nk x 2*nk)")
        if cov_matrix.shape[0] == 2 * nk and cov_matrix.shape[1] == 2 * nk:
            return cov_matrix
    
    # モノポールとクアドルポールの部分のみを抽出
    # インデックス: 0-nk-1 (P0), nk-2*nk-1 (P2)
    mono_quad_cov = np.zeros((2 * nk, 2 * nk))
    mono_quad_cov[0:nk, 0:nk] = cov_matrix[0:nk, 0:nk]  # P0-P0
    mono_quad_cov[0:nk, nk:2*nk] = cov_matrix[0:nk, nk:2*nk]  # P0-P2
    mono_quad_cov[nk:2*nk, 0:nk] = cov_matrix[nk:2*nk, 0:nk]  # P2-P0
    mono_quad_cov[nk:2*nk, nk:2*nk] = cov_matrix[nk:2*nk, nk:2*nk]  # P2-P2
    
    return mono_quad_cov


def read_covariance_matrix_mono_quad(filename, nk):
    """
    モノポールとクアドルポールのみを含む共分散行列を読み込む
    
    BOSS DR12の形式: 2*nk x 2*nk の行列
    順序: [P0(k1), ..., P0(nk), P2(k1), ..., P2(nk)]
    
    Returns:
        cov_matrix: 共分散行列 (2*nk x 2*nk)
    """
    data = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or len(line) == 0:
                continue
            try:
                row = [float(x) for x in line.split()]
                data.append(row)
            except ValueError:
                continue
    
    cov_matrix = np.array(data)
    
    expected_size = 2 * nk
    if cov_matrix.shape[0] != expected_size or cov_matrix.shape[1] != expected_size:
        raise ValueError(f"Expected {expected_size}x{expected_size} matrix, got {cov_matrix.shape}")
    
    return cov_matrix


def read_covariance_matrix_separate(mono_file, quad_file, nk):
    """
    モノポールとクアドルポールが別々のファイルとして提供される場合
    
    Returns:
        cov_matrix: 統合された共分散行列 (2*nk x 2*nk)
    """
    # モノポールの共分散行列を読み込む
    mono_cov = read_covariance_matrix_mono_quad(mono_file, nk)
    # クアドルポールの共分散行列を読み込む
    quad_cov = read_covariance_matrix_mono_quad(quad_file, nk)
    
    # 統合（クロス項は0と仮定、または別途提供される場合がある）
    # ここでは簡易的に対角ブロックのみを使用
    cov_matrix = np.zeros((2 * nk, 2 * nk))
    cov_matrix[0:nk, 0:nk] = mono_cov[0:nk, 0:nk]
    cov_matrix[nk:2*nk, nk:2*nk] = quad_cov[nk:2*nk, nk:2*nk]
    
    print("Warning: Cross-correlation between monopole and quadrupole is set to zero.")
    print("If cross-correlation data is available, please provide it separately.")
    
    return cov_matrix


def write_covariance_matrix(cov_matrix, output_file):
    """
    共分散行列をKaiser RSD用の形式で書き込む
    
    形式: [row_index] [col_index] [covariance_value]
    """
    n = cov_matrix.shape[0]
    print(f"Writing {n}x{n} covariance matrix to {output_file}...")
    
    with open(output_file, 'w') as f:
        for i in range(n):
            for j in range(n):
                f.write(f"{i}  {j}  {cov_matrix[i, j]:.6e}\n")
    
    print(f"Successfully wrote {n*n} covariance entries.")


def main():
    parser = argparse.ArgumentParser(
        description='Convert BOSS DR12 covariance matrix to Kaiser RSD format'
    )
    parser.add_argument('input_file', help='Input covariance matrix file')
    parser.add_argument('output_file', help='Output covariance matrix file')
    parser.add_argument('--nk', type=int, required=True, help='Number of k bins')
    parser.add_argument('--format', choices=['full', 'mono_quad', 'separate'], 
                       default='mono_quad',
                       help='Format of input covariance matrix')
    parser.add_argument('--quad-file', help='Quadrupole covariance file (for separate format)')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' not found.")
        sys.exit(1)
    
    print(f"Reading BOSS DR12 covariance matrix from: {args.input_file}")
    print(f"Format: {args.format}")
    print(f"Number of k bins: {args.nk}")
    
    try:
        if args.format == 'full':
            cov_matrix = read_covariance_matrix_full(args.input_file, args.nk)
        elif args.format == 'mono_quad':
            cov_matrix = read_covariance_matrix_mono_quad(args.input_file, args.nk)
        elif args.format == 'separate':
            if not args.quad_file:
                print("Error: --quad-file is required for separate format")
                sys.exit(1)
            cov_matrix = read_covariance_matrix_separate(args.input_file, args.quad_file, args.nk)
        
        print(f"\nCovariance matrix shape: {cov_matrix.shape}")
        print(f"  - Expected: ({2*args.nk}, {2*args.nk})")
        print(f"  - Monopole block: [0:{args.nk}, 0:{args.nk}]")
        print(f"  - Quadrupole block: [{args.nk}:{2*args.nk}, {args.nk}:{2*args.nk}]")
        
        write_covariance_matrix(cov_matrix, args.output_file)
        
        print("\n" + "="*60)
        print("Conversion completed successfully!")
        print("="*60)
        print(f"Output file: {args.output_file}")
        print(f"\nNext step: Update .dataset file to use this covariance matrix")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()

