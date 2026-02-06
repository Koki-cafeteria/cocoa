# Velocileptors YAML設定ガイド

## 概要

velocileptorsのパワースペクトルを計算・プロットするためのYAMLファイルの書き方です。

## 基本的な設定

### 1. .datasetファイルでの設定（推奨）

velocileptorsの設定値（`k_min`, `k_max`, `nk`, `z_eff`）は、**`.dataset`ファイルから読み込むことができます**。
これは`theta_min_arcmin`などと同様の方法です。

#### .datasetファイルへの追加例

```ini
# 既存の設定...
data_file = mg_hsc_y3_theory12.modelvector
cov_file = mg_hsc_y3_cov_transformed
# ... その他の設定 ...

# velocileptors設定（追加）
velocileptors_k_min = 5e-3      # kの最小値 [h/Mpc]
velocileptors_k_max = 0.3       # kの最大値 [h/Mpc]
velocileptors_nk = 50           # k-binの数
velocileptors_z_eff = 0.5       # 有効赤方偏移
```

**優先順位:**
1. `.dataset`ファイルの値（最優先）
2. YAMLファイルの値（.datasetにない場合）
3. デフォルト値（どちらにもない場合）

### 2. YAMLファイルでの設定

```yaml
likelihood:
  mg_hsc_y3.mg_hsc_y3_2x2pt:  # 既存のクラスを使用（データファイルが必要）
    path: ./external_modules/data/mg_hsc_y3
    data_file: mg_hsc_y3_M1_GGL0.05.dataset  # 最小限のデータファイルでOK
    
    # velocileptors統合設定（必須）
    use_velocileptors: True
    velocileptors_path: ./external_modules/code/velocileptors
    
    # 以下の設定は.datasetファイルから読み込まれる（.datasetにない場合のみYAMLの値が使われる）
    velocileptors_k_min: 5e-3      # kの最小値 [h/Mpc]（.dataset優先）
    velocileptors_k_max: 0.3       # kの最大値 [h/Mpc]（.dataset優先）
    velocileptors_nk: 50           # k-binの数（.dataset優先）
    velocileptors_z_eff: 0.5       # 有効赤方偏移（.dataset優先）
    
    # velocileptors用のバイアスパラメータのデフォルト値
    velocileptors_bias_defaults:
      b1: 2.0      # 線形バイアス
      b2: 0.0      # 2次バイアス
      bs: 0.0      # シアーバイアス
      b3: 0.0      # 3次バイアス
      alpha0: 0.0  # カウンター項（monopole）
      alpha2: 0.0  # カウンター項（quadrupole）
      alpha4: 0.0  # カウンター項（hexadecapole）
      SN0: 0.0     # 確率的項（monopole）
      SN2: 0.0     # 確率的項（quadrupole）
```

### 3. 重要な設定項目の説明

#### `use_velocileptors`
- **必須**: `True`に設定するとvelocileptors計算が有効になります
- **デフォルト**: `False`（既存の統計量への影響なし）
- **設定場所**: YAMLファイルのみ（.datasetファイルからは読み込まれない）

#### `velocileptors_path`
- velocileptorsのパス（通常は`./external_modules/code/velocileptors`）
- **設定場所**: YAMLファイルのみ（.datasetファイルからは読み込まれない）

#### `velocileptors_k_min`, `velocileptors_k_max`, `velocileptors_nk`
- k-binの設定
- **設定場所**: `.dataset`ファイル（推奨）またはYAMLファイル
- **優先順位**: `.dataset`ファイル > YAMLファイル > デフォルト値
- **推奨値**:
  - `k_min`: `5e-3` (0.005) [h/Mpc]
  - `k_max`: `0.3` [h/Mpc]
  - `nk`: `50` (k-binの数)

#### `velocileptors_z_eff`
- 有効赤方偏移（データの平均z）
- **設定場所**: `.dataset`ファイル（推奨）またはYAMLファイル
- **優先順位**: `.dataset`ファイル > YAMLファイル > デフォルト値（0.5）
- **例**: BOSS DR12 CMASSの場合は`z_eff ≈ 0.57`

#### `velocileptors_bias_defaults`
- velocileptors用のバイアスパラメータのデフォルト値
- **パラメータの意味**:
  - `b1`: 線形バイアス（通常 1.5-3.0）
  - `b2`: 2次バイアス（通常 0.0-1.0）
  - `bs`: シアーバイアス（通常 0.0）
  - `b3`: 3次バイアス（通常 0.0）
  - `alpha0`, `alpha2`, `alpha4`: カウンター項（通常 0.0）
  - `SN0`, `SN2`: 確率的項（通常 0.0）

### 4. バイアスパラメータの指定方法

#### 方法1: `velocileptors_bias_defaults`で指定（推奨）
```yaml
likelihood:
  mg_hsc_y3.mg_hsc_y3_2x2pt:
    velocileptors_bias_defaults:
      b1: 2.0
      b2: 0.0
      # ... その他のパラメータ ...
```

#### 方法2: paramsセクションで指定（将来の実装）
```yaml
params:
  velocileptors_b1:
    value: 2.0
  velocileptors_b2:
    value: 0.0
  # ... その他のパラメータ ...
```

#### 方法3: sampler.overrideで指定（将来の実装）
```yaml
sampler:
  evaluate:
    override:
      velocileptors_b1: 2.0
      velocileptors_b2: 0.0
      # ... その他のパラメータ ...
```

### 5. 結果の取得方法

velocileptors計算結果は、`set_cosmo_related()`内で以下の属性に保存されます：

```python
self.velocileptors_kv   # k値 [h/Mpc]
self.velocileptors_p0   # Monopole (P0)
self.velocileptors_p2   # Quadrupole (P2)
self.velocileptors_p4   # Hexadecapole (P4)
self.velocileptors_z_eff # 有効赤方偏移
```

これらの値を取得するには、Likelihoodクラス内で`self.velocileptors_kv`等にアクセスするか、
専用のLikelihoodクラスを作成する必要があります。

## 完全なYAMLファイル例

`EXAMPLE_EVALUATE_velocileptors.yaml`を参照してください。

## .datasetファイルの例

`projects/mg_hsc_y3/data/velocileptors_example.dataset`を参照してください。

既存の.datasetファイルに以下の行を追加することで、velocileptors設定を指定できます：

```ini
# velocileptors設定（追加）
velocileptors_k_min = 5e-3
velocileptors_k_max = 0.3
velocileptors_nk = 50
velocileptors_z_eff = 0.5
```

## プロット用の専用Likelihoodクラス（オプション）

プロット専用のLikelihoodクラスを作成する場合：

```python
# mg_hsc_y3_velocileptors.py
from cobaya.likelihoods.mg_hsc_y3._cosmolike_prototype_base import _cosmolike_prototype_base
import numpy as np

class mg_hsc_y3_velocileptors(_cosmolike_prototype_base):
    def initialize(self):
        super().initialize(probe="velocileptors")
        # 最小限のデータファイルでOK（ダミーでも可）
    
    def internal_get_datavector(self, **params_values):
        self.set_cosmo_related()
        
        # velocileptors結果を取得
        if hasattr(self, 'velocileptors_kv'):
            kv = self.velocileptors_kv
            p0 = self.velocileptors_p0
            p2 = self.velocileptors_p2
            p4 = self.velocileptors_p4
            
            # データベクトルを構築（例: 多極子を結合）
            datavector = np.concatenate([p0, p2, p4])
        else:
            # velocileptors計算が失敗した場合
            datavector = np.array([])
        
        return datavector
```

## トラブルシューティング

### velocileptorsが計算されない場合

1. `use_velocileptors: True`が設定されているか確認
2. `velocileptors_path`が正しいか確認
3. プローブが`velocileptors`を必要とするか確認（現在は`kaiser_rsd`等のみ）
4. ログで警告メッセージを確認

### 計算時間が長い場合

- `velocileptors_nk`を減らす（例: 50 → 30）
- `velocileptors_k_max`を減らす（例: 0.3 → 0.2）

### バイアスパラメータを変更したい場合

- `velocileptors_bias_defaults`で設定を変更
- または、将来の実装でparamsセクションから指定可能になる予定

