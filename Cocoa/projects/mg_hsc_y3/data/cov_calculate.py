# hsc_covariance.dat ファイルからデータを読み出す
data = []
with open('covariance_transformed.dat', 'r') as f:
    for line in f:
        # コメント行は無視する
        if line.startswith('#'):
            continue
        # 行をスペースで分割し、浮動小数点数に変換してリストに追加
        data.append([float(x) for x in line.strip().split()])

# hsc_y3_covariance.dat ファイルに特定の成分を書き込む
with open('mg_hsc_y3_cov_transformed', 'w') as f:
    for i in range(240):
        for j in range(240):
            # (i,j)成分のみを書き込む
            f.write(f"{i} {j} 0 0 0 0 0 0 {data[i][j]} 0\n")
