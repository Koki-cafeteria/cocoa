# Kaiser RSD データファイル形式

## 概要

Kaiser RSD（Redshift Space Distortion）の解析では、モノポール（ℓ=0）とクアドルポール（ℓ=2）の両方のパワースペクトルを計算します。

## データベクトルの構造

データベクトルは以下の順序で格納されます（3x2ptと同様に、各プローブが塊になっています）：

```
[P_0(k1,z1), P_0(k1,z2), ..., P_0(k1,zNz), 
 P_0(k2,z1), P_0(k2,z2), ..., P_0(k2,zNz),
 ...
 P_0(kNk,z1), P_0(kNk,z2), ..., P_0(kNk,zNz),
 P_2(k1,z1), P_2(k1,z2), ..., P_2(k1,zNz),
 P_2(k2,z1), P_2(k2,z2), ..., P_2(k2,zNz),
 ...
 P_2(kNk,z1), P_2(kNk,z2), ..., P_2(kNk,zNz)]
```

つまり、まずすべてのモノポール（ℓ=0）が格納され、その後にすべてのクアドルポール（ℓ=2）が格納されます。これは3x2ptのデータベクトル構造（cosmic shear、g-g lensing、galaxy clusteringがそれぞれ塊になっている）と同様の形式です。

データベクトルのサイズ: `2 * N_k * N_z`

## ファイル形式

### 1. データベクトルファイル (`kaiser_rsd_datavector.txt`)

形式: 2列（インデックス、値）

```
0  P_0(k1,z1)
1  P_0(k1,z2)
...
N_k*N_z-1  P_0(kNk,zNz)
N_k*N_z  P_2(k1,z1)
N_k*N_z+1  P_2(k1,z2)
...
2*N_k*N_z-1  P_2(kNk,zNz)
```

### 2. マスクファイル (`kaiser_rsd_mask.txt`)

形式: 2列（インデックス、マスク値）

- マスク値: 0（使用しない）または 1（使用する）

```
0  1
1  1
2  1
3  1
...
```

### 3. 共分散行列ファイル (`kaiser_rsd_covariance.txt`)

形式: 3列（行インデックス、列インデックス、共分散値）

```
0  0  C_00
0  1  C_01
0  2  C_02
...
1  0  C_10
1  1  C_11
...
```

### 4. データセット設定ファイル (`kaiser_rsd_data.dataset`)

INI形式の設定ファイル。以下のパラメータを設定：

- `Nk`: kビンの数
- `k_min`: 最小k値 [h/Mpc]
- `k_max`: 最大k値 [h/Mpc]
- `Nz`: 赤方偏移ビンの数
- `z_min`: 最小赤方偏移
- `z_max`: 最大赤方偏移

## SDSS/BOSS DR12などの実際のデータを使用する場合

### BOSS DR12データの変換

BOSS DR12のパワースペクトルデータ（`Power_Spectrum_cmass_ngc_v5`形式）を変換するには、`convert_boss_dr12_data.py`スクリプトを使用します：

```bash
cd /home/tanida/cocoa_v41/cocoa/Cocoa/projects/mg_hsc_y3/data
python convert_boss_dr12_data.py Power_Spectrum_cmass_ngc_v5.txt boss_dr12_cmass
```

このスクリプトは以下を生成します：
- `boss_dr12_cmass_datavector.txt`: データベクトルファイル
- `boss_dr12_cmass_mask.txt`: マスクファイル
- `boss_dr12_cmass_k_bins.txt`: kビン情報（参考用）

**注意**: 共分散行列ファイルは別途準備する必要があります。BOSS DR12のデータリリースから取得し、`convert_boss_dr12_covariance.py`スクリプトで変換してください。

### データセット設定ファイルの更新

BOSS DR12データを使用する場合、`.dataset`ファイルを以下のように更新します：

```ini
# kビンの設定（BOSS DR12のkビン数に合わせる）
Nk = 40  # 実際のkビン数に合わせて調整
k_min = 0.005  # 最小k値 [h/Mpc]
k_max = 0.395  # 最大k値 [h/Mpc]

# 赤方偏移ビンの設定（BOSS DR12 CMASSは単一のzビン）
Nz = 1
z_min = 0.43  # BOSS DR12 CMASSのz範囲に合わせて調整
z_max = 0.70

# データファイルのパス
data_file = boss_dr12_cmass_datavector.txt
cov_file = boss_dr12_cmass_covariance.txt  # 別途準備が必要
mask_file = boss_dr12_cmass_mask.txt
```

### 一般的なデータの準備手順

1. **モノポールとクアドルポールのパワースペクトルを取得**
   - SDSS/BOSSなどの観測データから、各kビンとzビンに対するP_0(k,z)とP_2(k,z)を取得
   - 単位: [Mpc/h]^3
   - 注意: ショットノイズを引いた後の値を使用する場合と、引く前の値を使用する場合がある

2. **データベクトルファイルの作成**
   - 上記の形式に従って、まずすべてのモノポール、その後にすべてのクアドルポールを並べる
   - インデックスは0から始まる連番

3. **共分散行列の作成**
   - モノポールとクアドルポールの間の相関も含む
   - サイズ: (2*N_k*N_z) × (2*N_k*N_z)
   - 形式: `[row_index] [col_index] [covariance_value]`

4. **マスクファイルの作成**
   - 使用するデータポイントを1、使用しないデータポイントを0に設定

### 注意事項

- データの単位（特にkの単位）が理論計算と一致していることを確認（通常は[h/Mpc]）
- 共分散行列は正定値である必要がある
- データベクトルと共分散行列のサイズが一致していることを確認
- BOSS DR12の場合は単一の赤方偏移ビン（Nz=1）になる

## 例: ダミーデータの生成

`generate_dummy_files.py`スクリプトを使用して、テスト用のダミーデータを生成できます：

```bash
cd /home/tanida/cocoa_v41/cocoa/Cocoa/projects/mg_hsc_y3/data
python generate_dummy_files.py
```

このスクリプトは、モノポールとクアドルポールの両方を含むデータファイルを生成します。

## 例: BOSS DR12データの変換

BOSS DR12のパワースペクトルデータを変換する例：

```bash
cd /home/tanida/cocoa_v41/cocoa/Cocoa/projects/mg_hsc_y3/data
python convert_boss_dr12_data.py Power_Spectrum_cmass_ngc_v5.txt boss_dr12_cmass
```

変換後、`.dataset`ファイルを更新して、生成されたデータファイルを指定します。

### BOSS DR12共分散行列の変換

BOSS DR12のデータリリースから共分散行列を取得したら、`convert_boss_dr12_covariance.py`スクリプトで変換します：

```bash
# モノポール+クアドルポール形式の場合（推奨）
python convert_boss_dr12_covariance.py \
    BOSS_DR12_covariance_matrix.txt \
    boss_dr12_cmass_covariance.txt \
    --nk 40 \
    --format mono_quad

# モノポール+クアドルポール+ヘキサデカポール形式の場合
python convert_boss_dr12_covariance.py \
    BOSS_DR12_full_covariance_matrix.txt \
    boss_dr12_cmass_covariance.txt \
    --nk 40 \
    --format full
```

**共分散行列の形式**:
- `mono_quad`: 2*nk x 2*nk の行列（モノポール+クアドルポール）
- `full`: 3*nk x 3*nk の行列（モノポール+クアドルポール+ヘキサデカポール）から必要な部分を抽出

**BOSS DR12データリリースの場所**:
- BOSS DR12の公式データリリースページから共分散行列ファイルをダウンロード
- 通常、パワースペクトルデータと同じディレクトリに含まれています

