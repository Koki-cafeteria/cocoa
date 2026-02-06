import numpy as np

# --- 設定 ---
# .datasetファイルで指定した値と一致させる
Nk = 100
Nz = 10
# モノポール（ℓ=0）とクアドルポール（ℓ=2）の両方を含む
# データベクトルの構造（3x2ptと同様に塊になっている）:
# [P_0(k1,z1), P_0(k1,z2), ..., P_0(kNk,zNz), P_2(k1,z1), P_2(k1,z2), ..., P_2(kNk,zNz)]
# まずすべてのモノポール、その後にすべてのクアドルポール
N_data = 2 * Nk * Nz  # モノポール + クアドルポール

# --- ファイル名 ---
data_filename = "kaiser_rsd_datavector.txt"
cov_filename = "kaiser_rsd_covariance.txt"
mask_filename = "kaiser_rsd_mask.txt"

# --- 1. データベクトルファイルの作成 ---
# 書式: [インデックス] [値]
# テスト用なので、値はすべて0.0でOK
print(f"Generating {data_filename} with {N_data} data points (monopole + quadrupole)...")
indices = np.arange(N_data)
values = np.zeros(N_data)
data_vector = np.column_stack((indices, values))
np.savetxt(data_filename, data_vector, fmt='%d  %.6e')

# --- 2. マスクファイルの作成 ---
# 書式: [インデックス] [マスク値 (0 or 1)]
# 全てのデータ点を使用するため、値はすべて1にする
print(f"Generating {mask_filename}...")
mask_values = np.ones(N_data, dtype=int)
mask_file_content = np.column_stack((indices, mask_values))
np.savetxt(mask_filename, mask_file_content, fmt='%d  %d')

# --- 3. 共分散行列ファイルの作成 ---
# 書式: [行インデックス i] [列インデックス j] [共分散 C_ij]
# テスト用なので、相関のない対角行列を作成 (対角成分にのみ値を入れる)
# ゼロ割を避けるため、小さな値 (例: 1.0) を入れておく
print(f"Generating {cov_filename} ({N_data}x{N_data} matrix)...")
cov_values = np.ones(N_data) # 対角成分の値
covariance_matrix = np.column_stack((indices, indices, cov_values))
np.savetxt(cov_filename, covariance_matrix, fmt='%d  %d  %.6e')

print("\nAll dummy files have been generated successfully!")
print(f"  - Data vector: {N_data} points (monopole + quadrupole)")
print(f"  - Covariance matrix: {N_data}x{N_data}")
print(f"  - Mask: {N_data} points")
