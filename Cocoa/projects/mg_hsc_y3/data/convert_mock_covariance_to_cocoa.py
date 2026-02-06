#!/usr/bin/env python3
"""
モックデータから作成した共分散行列をCocoa形式に変換するスクリプト

モックデータの形式: [P0(k1), P2(k1), P0(k2), P2(k2), ...] (交互)
Cocoa形式: [P0(k1), ..., P0(kNk), P2(k1), ..., P2(kNk)] (塊)

使用方法:
    python convert_mock_covariance_to_cocoa.py \
        input_covariance.npy \
        output_covariance.txt \
        --nk N_K_BINS
"""

import numpy as np
import sys
import argparse

def convert_interleaved_to_block_format(cov_interleaved, nk):
    """
    交互形式の共分散行列を塊形式に変換
    
    Parameters:
        cov_interleaved: 交互形式の共分散行列 (2*nk x 2*nk)
                        形式: [P0(k1), P2(k1), P0(k2), P2(k2), ...]
        nk: kビンの数
    
    Returns:
        cov_block: 塊形式の共分散行列 (2*nk x 2*nk)
                  形式: [P0(k1), ..., P0(kNk), P2(k1), ..., P2(kNk)]
    """
    if cov_interleaved.shape[0] != 2 * nk or cov_interleaved.shape[1] != 2 * nk:
        raise ValueError(f"Expected (2*nk, 2*nk) = ({2*nk}, {2*nk}) matrix, got {cov_interleaved.shape}")
    
    # 変換行列を作成
    # 交互形式のインデックス: [0, 1, 2, 3, ..., 2*nk-2, 2*nk-1]
    # 塊形式のインデックス: [0, 2, 4, ..., 2*nk-2, 1, 3, 5, ..., 2*nk-1]
    
    # インデックスのマッピング
    # 交互形式 -> 塊形式
    # P0(k1) [0] -> [0]
    # P2(k1) [1] -> [nk]
    # P0(k2) [2] -> [1]
    # P2(k2) [3] -> [nk+1]
    # ...
    
    cov_block = np.zeros_like(cov_interleaved)
    
    # モノポール部分 (0 to nk-1)
    for i in range(nk):
        old_idx_i = 2 * i  # 交互形式でのP0(ki)のインデックス
        new_idx_i = i      # 塊形式でのP0(ki)のインデックス
        
        for j in range(nk):
            old_idx_j = 2 * j
            new_idx_j = j
            cov_block[new_idx_i, new_idx_j] = cov_interleaved[old_idx_i, old_idx_j]
    
    # クアドルポール部分 (nk to 2*nk-1)
    for i in range(nk):
        old_idx_i = 2 * i + 1  # 交互形式でのP2(ki)のインデックス
        new_idx_i = nk + i     # 塊形式でのP2(ki)のインデックス
        
        for j in range(nk):
            old_idx_j = 2 * j + 1
            new_idx_j = nk + j
            cov_block[new_idx_i, new_idx_j] = cov_interleaved[old_idx_i, old_idx_j]
    
    # クロス項 (P0-P2)
    for i in range(nk):
        old_idx_i = 2 * i      # 交互形式でのP0(ki)
        new_idx_i = i         # 塊形式でのP0(ki)
        
        for j in range(nk):
            old_idx_j = 2 * j + 1  # 交互形式でのP2(kj)
            new_idx_j = nk + j     # 塊形式でのP2(kj)
            cov_block[new_idx_i, new_idx_j] = cov_interleaved[old_idx_i, old_idx_j]
            cov_block[new_idx_j, new_idx_i] = cov_interleaved[old_idx_j, old_idx_i]
    
    return cov_block


def write_cocoa_covariance(cov_matrix, output_file):
    """
    Cocoa形式で共分散行列を書き込む
    
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
        description='Convert mock covariance matrix from interleaved to Cocoa block format'
    )
    parser.add_argument('input_file', help='Input covariance matrix file (.npy or .txt)')
    parser.add_argument('output_file', help='Output covariance matrix file (Cocoa format)')
    parser.add_argument('--nk', type=int, required=True, help='Number of k bins')
    parser.add_argument('--format', choices=['npy', 'txt'], default='npy',
                       help='Input file format (default: npy)')
    
    args = parser.parse_args()
    
    print(f"Reading covariance matrix from: {args.input_file}")
    print(f"Format: {args.format}")
    print(f"Number of k bins: {args.nk}")
    
    # 共分散行列を読み込む
    try:
        if args.format == 'npy':
            cov_interleaved = np.load(args.input_file)
        else:
            cov_interleaved = np.loadtxt(args.input_file)
        
        print(f"Input covariance matrix shape: {cov_interleaved.shape}")
        
        # 交互形式から塊形式に変換
        print("Converting from interleaved format to block format...")
        cov_block = convert_interleaved_to_block_format(cov_interleaved, args.nk)
        
        print(f"Output covariance matrix shape: {cov_block.shape}")
        print(f"  - Monopole block: [0:{args.nk}, 0:{args.nk}]")
        print(f"  - Quadrupole block: [{args.nk}:{2*args.nk}, {args.nk}:{2*args.nk}]")
        
        # Cocoa形式で書き込む
        write_cocoa_covariance(cov_block, args.output_file)
        
        print("\n" + "="*60)
        print("Conversion completed successfully!")
        print("="*60)
        print(f"Output file: {args.output_file}")
        print(f"\nNext step: Update .dataset file to use this covariance matrix")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

