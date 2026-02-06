import numpy as np

# --- 設定 ---
# 3x2ptのデータ点数
N_3x2pt = 240
# P(k,z)のデータ点数
N_pk = 1000
# 合計のデータ点数
N_total = N_3x2pt + N_pk

# --- 入力ファイル名 ---
# 3x2ptで使っていた既存のマスクファイル
original_3x2pt_mask_file = "mg_hsc_y3_M1_GGLOLAP0.05.mask"

# --- 出力ファイル名 ---
combined_mask_file = "combined_mask.txt"
combined_data_file = "combined_datavector.txt"
combined_cov_file = "combined_covariance.txt"


# --- 1. 結合マスクファイルの作成 ---
print(f"Generating {combined_mask_file}...")
# 既存の3x2ptマスクを読み込む
mask_3x2pt = np.loadtxt(original_3x2pt_mask_file, usecols=1, dtype=int)

# P(k)部分のマスクを作成 (全て1)
mask_pk = np.ones(N_pk, dtype=int)

# 2つのマスクを結合
combined_mask_values = np.concatenate([mask_3x2pt, mask_pk])

# ファイルに書き出す
indices = np.arange(N_total)
mask_output = np.column_stack((indices, combined_mask_values))
np.savetxt(combined_mask_file, mask_output, fmt='%d  %d')


# --- 2. 結合データベクトル（ダミー）の作成 ---
print(f"Generating {combined_data_file}...")
data_values = np.zeros(N_total)
data_output = np.column_stack((indices, data_values))
np.savetxt(combined_data_file, data_output, fmt='%d  %.6e')


# --- 3. 結合共分散行列（ダミー）の作成 ---
print(f"Generating {combined_cov_file}...")
cov_values = np.ones(N_total) # 対角成分に1を入れる
cov_output = np.column_stack((indices, indices, cov_values))
np.savetxt(combined_cov_file, cov_output, fmt='%d  %d  %.6e')


print(f"\nSuccessfully generated combined files for {N_total} data points!")
