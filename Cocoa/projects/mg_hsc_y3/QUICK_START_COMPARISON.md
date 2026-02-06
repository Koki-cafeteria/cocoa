# Kaiser RSD vs Velocileptors 比較 - Quick Start

## 最速の実行方法

### 自動実行スクリプトを使用（推奨）

```bash
cd /home/tanida/cocoa_v41/cocoa/Cocoa
chmod +x projects/mg_hsc_y3/scripts/run_comparison.sh
./projects/mg_hsc_y3/scripts/run_comparison.sh
```

このスクリプトは以下を自動で実行します：
1. 既存の出力ファイルをチェック
2. 必要に応じてKaiser RSDを実行
3. 必要に応じてVelocileptersを実行
4. 比較プロットを作成

### 手動実行

```bash
cd /home/tanida/cocoa_v41/cocoa/Cocoa
source start_cocoa.sh

# 1. Kaiser RSDを実行
cobaya-run projects/mg_hsc_y3/EXAMPLE_EVALUATE_kaiser_RSD.yaml

# 2. Velocileptorsを実行
cobaya-run projects/mg_hsc_y3/EXAMPLE_EVALUATE_velocileptors.yaml

# 3. 比較プロットを作成
python projects/mg_hsc_y3/scripts/compare_kaiser_velocileptors.py
```

## 出力ファイル

比較プロットは以下に保存されます：

1. **単一z値での詳細比較**:
   ```
   projects/mg_hsc_y3/data/comparison_kaiser_velocileptors_single_z.png
   ```
   - P0 (Monopole) の比較
   - P2 (Quadrupole) の比較
   - P4 (Hexadecapole) - Velocileptorsのみ
   - 相対差のプロット

2. **複数z値での比較**:
   ```
   projects/mg_hsc_y3/data/comparison_kaiser_velocileptors_multiple_z.png
   ```
   - Kaiser RSDの複数のz値
   - Velocileptorsの単一z値（z_eff=0.5）

## 所要時間の目安

- Kaiser RSD: 約5-10分
- Velocileptors: 約3-5分
- プロット作成: 数秒

**合計**: 約10-15分

## 比較のポイント

### 低k領域（k < 0.1 h/Mpc）

- 両モデルは比較的一致するはず
- 線形領域ではKaiserモデルも精度が高い

### 高k領域（k > 0.1 h/Mpc）

- Velocileptorsは非線形効果を含む
- Kaiserは線形理論のため精度が低下
- **相対差が大きくなることが期待される**

### P2 (Quadrupole)

- RSD（赤方偏移空間歪み）の効果を示す
- 高k領域での違いがP0よりも顕著
- Velocileptorsの方が実際の観測データに近いはず

## パラメータのカスタマイズ

### 宇宙論パラメータを変更する場合

両方のYAMLファイルの`sampler.evaluate.override`セクションを編集：

```yaml
sampler:
  evaluate:
    N: 1
    override:
      mu0: 0.0      # 修正重力パラメータ（0.0 = GR）
      sigma0: 0.0   # 修正重力パラメータ
      As_1e9: 2.2065
      ns: 0.9645
      H0: 67.089
      omegab: 0.049434
      omegam: 0.3156
```

### バイアスパラメータを変更する場合

**Kaiser RSD** (`EXAMPLE_EVALUATE_kaiser_RSD.yaml`):
```yaml
sampler:
  evaluate:
    override:
      hsc_B_KAISER_RSD_1: 2.0  # 線形バイアス
```

**Velocileptors** (`EXAMPLE_EVALUATE_velocileptors.yaml`):
```yaml
sampler:
  evaluate:
    override:
      velocileptors_b1plus1_sigma8: 1.36   # (1+b1)*sigma8
      velocileptors_b2_sigma8_sq: 0.32     # b2*sigma8^2
      velocileptors_bs_sigma8_sq: -0.192   # bs*sigma8^2
      velocileptors_alpha0: 10.0           # カウンター項
      velocileptors_alpha2: 20.0
      velocileptors_alpha4: -60.0
      velocileptors_SN0: 0.0               # 確率的項
      velocileptors_SN2: 0.0
```

### 公正な比較のための設定

Kaiserと同じ線形バイアスのみを使用する場合：

```yaml
# Velocileptorsで高次項を0に設定
velocileptors_b2_sigma8_sq: 0.0
velocileptors_bs_sigma8_sq: 0.0
velocileptors_alpha0: 0.0
velocileptors_alpha2: 0.0
velocileptors_alpha4: 0.0
velocileptors_SN0: 0.0
velocileptors_SN2: 0.0
```

## トラブルシューティング

### エラー: `ModuleNotFoundError: No module named 'velocileptors'`

Velocileptorsがインストールされていません：

```bash
cd /home/tanida/cocoa_v41/cocoa/Cocoa
source start_cocoa.sh
./installation_scripts/setup_velocileptors.sh
./installation_scripts/compile_velocileptors.sh
```

### エラー: `File not found: kaiser_rsd_data.modelvector`

データファイルが見つかりません。YAMLファイルのパスを確認してください。

### プロットが表示されない

ヘッドレス環境（SSHなど）の場合、プロットは自動的に保存されます：

```bash
# 保存された画像を確認
ls -lh projects/mg_hsc_y3/data/comparison_*.png
```

## 次のステップ

1. **結果の解釈**:
   - 詳細は`KAISER_VELOCILEPTORS_COMPARISON_GUIDE.md`を参照

2. **パラメータスタディ**:
   - バイアスパラメータを変えて比較
   - 修正重力パラメータ（mu0, sigma0）の影響を調査

3. **実データとの比較**:
   - BOSS DR12 CMASSの観測データと比較
   - chi2を計算して最適なパラメータを求める

4. **論文用のプロット**:
   - `compare_kaiser_velocileptors.py`をカスタマイズ
   - フォントサイズ、カラーマップ、ラベルなどを調整

## 参考資料

- 詳細ガイド: `KAISER_VELOCILEPTORS_COMPARISON_GUIDE.md`
- Velocileptors設定: `VELOCILEPTORS_YAML_GUIDE.md`
- 比較スクリプト: `scripts/compare_kaiser_velocileptors.py`



