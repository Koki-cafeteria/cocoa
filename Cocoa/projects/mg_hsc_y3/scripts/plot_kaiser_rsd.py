#!/usr/bin/env python3
"""
Kaiser RSD理論データベクトルからパワースペクトルをプロット
縦軸: P(k), 横軸: k
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# データファイルのパス
theory_file = "./projects/mg_hsc_y3/data/mg_hsc_y3_theory_kaiser_rsd.modelvector"
dataset_file = "./projects/mg_hsc_y3/data/kaiser_rsd_data.dataset"

# データセットファイルからビニング情報を読み込む
def read_dataset_config(dataset_file):
    """データセットファイルから設定を読み込む"""
    config = {}
    with open(dataset_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()
                try:
                    # 数値に変換を試みる
                    if '.' in value:
                        config[key] = float(value)
                    else:
                        config[key] = int(value)
                except ValueError:
                    config[key] = value
    return config

# 理論データベクトルを読み込む
def load_theory_vector(theory_file):
    """理論データベクトルを読み込む"""
    data = np.loadtxt(theory_file)
    if data.ndim == 2:
        return data[:, 1]  # 2列目がデータ
    else:
        return data

# メイン処理
if __name__ == "__main__":
    # データセット設定を読み込む
    config = read_dataset_config(dataset_file)
    Nk = config.get('Nk', 100)
    k_min = config.get('k_min', 0.01)
    k_max = config.get('k_max', 1.0)
    Nz = config.get('Nz', 10)
    z_min = config.get('z_min', 0.1)
    z_max = config.get('z_max', 1.0)
    
    print(f"設定: Nk={Nk}, k_min={k_min}, k_max={k_max}, Nz={Nz}, z_min={z_min}, z_max={z_max}")
    
    # 理論データベクトルを読み込む
    theory_vector = load_theory_vector(theory_file)
    print(f"理論データベクトルのサイズ: {len(theory_vector)}")
    print(f"期待されるサイズ (Nk * Nz): {Nk * Nz}")
    
    # kとzのグリッドを作成
    k_bins = np.logspace(np.log10(k_min), np.log10(k_max), Nk)
    z_bins = np.linspace(z_min, z_max, Nz)
    
    # データベクトルを2次元配列に再構成 (k, z)
    # データベクトルの順序: index = i * Nz + j (i: k index, j: z index)
    Pkz = theory_vector.reshape(Nk, Nz)
    
    # プロット作成
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # 各zについてプロット
    colors = plt.cm.viridis(np.linspace(0, 1, Nz))
    for j, z in enumerate(z_bins):
        ax.loglog(k_bins, Pkz[:, j], 
                 color=colors[j], 
                 label=f'z = {z:.2f}',
                 linewidth=1.5)
    
    ax.set_xlabel('k [h/Mpc]', fontsize=14)
    ax.set_ylabel('P(k) [(Mpc/h)$^3$]', fontsize=14)
    ax.set_title('Kaiser RSD Power Spectrum', fontsize=16)
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, ncol=1)
    
    plt.tight_layout()
    
    # 保存
    output_file = "./projects/mg_hsc_y3/data/kaiser_rsd_power_spectrum.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"プロットを保存しました: {output_file}")
    
    # 表示
    plt.show()

