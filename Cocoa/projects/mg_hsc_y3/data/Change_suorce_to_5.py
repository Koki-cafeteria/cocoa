import numpy as np

# 入力ファイル名と出力ファイル名
input_filename = "mg_hsc_y3_source.nz"
output_filename = "mg_hsc_y3_source5bin.nz"

# データを読み込む
data = np.loadtxt(input_filename)

# 2列目のデータを取得
second_column = data[:, 1]

# 新しいデータを作成
# 元のデータ + 2列目を3~6列目にコピーしたデータ
new_data = np.hstack((data, second_column[:, None], second_column[:, None], second_column[:, None], second_column[:, None]))

# 新しいデータを書き込む
np.savetxt(output_filename, new_data, fmt="%.18e", delimiter=" ")

print(f"処理が完了しました。新しいデータは {output_filename} に保存されています。")