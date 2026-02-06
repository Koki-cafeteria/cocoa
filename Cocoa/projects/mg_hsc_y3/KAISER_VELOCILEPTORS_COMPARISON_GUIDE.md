# Kaiser RSD vs Velocileptors 比較ガイド

## 概要

KaiserモデルとVelocileptersの計算結果（P0, P2）を比較するためのガイドです。

## 実行手順

### 1. Kaiser RSDの計算

```bash
cd /home/tanida/cocoa_v41/cocoa/Cocoa
source start_cocoa.sh
cobaya-run projects/mg_hsc_y3/EXAMPLE_EVALUATE_kaiser_RSD.yaml
```

**出力ファイル**:
- `projects/mg_hsc_y3/data/mg_hsc_y3_theory_kaiser_rsd_boss_dr12.modelvector`

**データ構造**:
- Monopole (P0) と Quadrupole (P2) が連続して並ぶ
- `[P0(k1,z1), ..., P0(kNk,zNz), P2(k1,z1), ..., P2(kNk,zNz)]`
- BOSS DR12 CMASSのk-binsを使用（Nk ≈ 24）
- 複数のz値（Nz=10, z=0.1-1.0）

### 2. Velocileptorsの計算

```bash
cobaya-run projects/mg_hsc_y3/EXAMPLE_EVALUATE_velocileptors.yaml
```

**出力ファイル**:
- `projects/mg_hsc_y3/data/velocileptors_pk_multipoles_mu10.txt`

**データ構造**:
- k値ごとに3つの多極モーメント（P0, P2, P4）
- `[k, P0, P2, P4]` の4列
- 単一のz値（z_eff=0.5）
- k範囲: 5e-3 - 0.3 h/Mpc, Nk=50

### 3. 比較プロットの作成

```bash
python projects/mg_hsc_y3/scripts/compare_kaiser_velocileptors.py
```

**出力ファイル**:
- `projects/mg_hsc_y3/data/comparison_kaiser_velocileptors_single_z.png`
  - 単一のz値での詳細な比較（P0, P2, P4, 相対差）
- `projects/mg_hsc_y3/data/comparison_kaiser_velocileptors_multiple_z.png`
  - Kaiser RSDの複数z値とVelocileptersの比較

## プロットの内容

### 単一z値での比較 (`comparison_kaiser_velocileptors_single_z.png`)

4つのサブプロット:
1. **P0 (Monopole)**: Kaiser RSD vs Velocileptors
2. **P2 (Quadrupole)**: Kaiser RSD vs Velocileptors
3. **P4 (Hexadecapole)**: Velocileptorsのみ（Kaiserは計算しない）
4. **相対差**: (Velocileptors - Kaiser) / Kaiser [%]

### 複数z値での比較 (`comparison_kaiser_velocileptors_multiple_z.png`)

2つのサブプロット:
1. **P0 (Monopole)**: Kaiser RSD（複数z）とVelocileptors（z_eff=0.5）
2. **P2 (Quadrupole)**: Kaiser RSD（複数z）とVelocileptors（z_eff=0.5）

## パラメータの対応関係

### 宇宙論パラメータ（両方で共通）

```yaml
mu0: 1.0           # 修正重力パラメータ
sigma0: 0.0        # 修正重力パラメータ
As_1e9: 2.2065     # 原始パワースペクトル振幅
ns: 0.9645         # スペクトル指数
H0: 67.089         # ハッブル定数
omegab: 0.049434   # バリオン密度
omegam: 0.3156     # 物質密度
```

### バイアスパラメータ

#### Kaiser RSD
- `hsc_B_KAISER_RSD_1: 2.0` (線形バイアスのみ)

#### Velocileptors（DESIスタイル）
- `velocileptors_b1plus1_sigma8: 1.36` → b1 = (1.36/0.8) - 1 = 0.7
- `velocileptors_b2_sigma8_sq: 0.32` → b2 = 0.32/0.64 = 0.5
- `velocileptors_bs_sigma8_sq: -0.192` → bs = -0.192/0.64 = -0.3
- `velocileptors_alpha0: 10.0` (カウンター項)
- `velocileptors_alpha2: 20.0` (カウンター項)
- `velocileptors_alpha4: -60.0` (カウンター項)
- `velocileptors_SN0: 0.0` (確率的項)
- `velocileptors_SN2: 0.0` (確率的項)

## 比較時の注意点

1. **バイアスパラメータの違い**:
   - KaiserはP線形バイアス（b1）のみ
   - Velocileptorsは高次バイアス（b2, bs）とカウンター項を含む
   - 公正な比較のため、Velocileptorsの高次項を0に設定することも可能

2. **k範囲の違い**:
   - Kaiser: BOSS DR12の実際のk-bins（~24点）
   - Velocileptors: 均等なk-bins（50点）
   - 比較スクリプトは自動的に補間します

3. **z値の違い**:
   - Kaiser: 複数のz値（Nz=10）
   - Velocileptors: 単一のz_eff（0.5）
   - 比較スクリプトは最も近いz値を選択します

## 結果の解釈

### 期待される違い

1. **低k領域（k < 0.1 h/Mpc）**:
   - 両者は比較的一致するはず（線形領域）
   - 差は主に線形バイアスの違いによる

2. **高k領域（k > 0.1 h/Mpc）**:
   - Velocileptorsは非線形効果を含むため、違いが大きくなる
   - Kaiserは単純な線形理論なので、小スケールで精度が悪化

3. **P2 (Quadrupole)**:
   - RSD効果の強さを示す
   - 高k領域での違いがP0よりも顕著

### 統計情報

比較スクリプトは以下を出力します:
- **平均相対差**: 全k範囲での平均的な差
- **RMS相対差**: 差の二乗平均平方根

## トラブルシューティング

### ファイルが見つからない場合

```bash
# Kaiser RSDを再実行
cobaya-run projects/mg_hsc_y3/EXAMPLE_EVALUATE_kaiser_RSD.yaml

# Velocileptorsを再実行
cobaya-run projects/mg_hsc_y3/EXAMPLE_EVALUATE_velocileptors.yaml
```

### パラメータを変更したい場合

1. YAMLファイルの`sampler.evaluate.override`セクションを編集
2. 再度`cobaya-run`を実行
3. 比較スクリプトを実行

### k範囲を揃えたい場合

**Option 1**: Velocileptorsのk範囲をKaiserに合わせる

`EXAMPLE_EVALUATE_velocileptors.yaml`を編集:
```yaml
velocileptors_k_min: 0.018  # Kaiserのk_min
velocileptors_k_max: 0.3    # Kaiserのk_max
velocileptors_nk: 24        # Kaiserのk点数
```

**Option 2**: Kaiserのk-binsファイルを直接使用

`.dataset`ファイルで`k_bins_file`を指定することで、
同じk-binsを使用できます。

## 参考文献

- Kaiser (1987): "Clustering in real space and in redshift space"
- Martin White et al. (2020): "Velocileptors: Lagrangian Perturbation Theory for galaxy clustering"
- DES Y3 3x2pt analysis (2022): バイアスパラメータの事前分布



