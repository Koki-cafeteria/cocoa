# Velocileptors統合計画（オプション3: Likelihood内実装）

## 1. `_cosmolike_prototype_base.py`の役割

### 基底クラスとしての役割
`_cosmolike_prototype_base.py`は、**mg_hsc_y3プロジェクトの全Likelihoodクラスの共通基底クラス**です。

### 他のファイルとの関係

#### 継承関係
```
_cosmolike_prototype_base (基底クラス)
  ├── mg_hsc_y3_2x2pt.py      (2x2pt解析用)
  ├── mg_hsc_y3_3x2pt.py      (3x2pt解析用)
  ├── mg_hsc_y3_cosmic_shear.py (Cosmic Shear単独)
  ├── mg_hsc_y3_kaiser_RSD.py   (Kaiser RSD単独)
  ├── mg_hsc_y3_xi_gg.py        (Galaxy-Galaxy相関)
  └── mg_hsc_y3_xi_ggl.py       (Galaxy-Galaxy Lensing)
```

#### 各ファイルの役割

**`_cosmolike_prototype_base.py`（基底クラス）**
- **共通処理の実装**:
  - `initialize()`: CosmoLike C++インターフェースの初期化、データファイル読み込み、z/kグリッド設定
  - `get_requirements()`: Cobaya theory provider（MGCAMB）への要求定義
  - `set_cosmo_related()`: MGCAMBから宇宙論パラメータ・P(k)を取得し、CosmoLike C++に渡す
  - `set_lens_related()`: レンズ銀河のnuisanceパラメータ（バイアス、photo-z等）設定
  - `set_source_related()`: ソース銀河のnuisanceパラメータ（IA、shear calibration等）設定
  - `compute_logp()`: データベクトルからχ²を計算

**各Likelihoodクラス（例: `mg_hsc_y3_2x2pt.py`）**
- **プローブ固有の処理**:
  - `initialize()`: `probe="2x2pt"`などを指定して基底クラスを初期化
  - `internal_get_datavector()`: プローブ固有のデータベクトル計算
  - `logp()`, `get_datavector()`: Cobayaインターフェース

### 設計パターン
- **Template Method パターン**: 基底クラスが処理の骨組みを定義し、派生クラスが詳細を実装
- **DRY原則**: 共通処理（CosmoLike初期化、P(k)取得、キャッシュ管理）を1箇所に集約

---

## 2. Velocileptors統合実装計画

### 2.1 実装方針
- **Cobaya本体は変更しない**（オプション3）
- `_cosmolike_prototype_base.py`にvelocileptors統合機能を追加
- YAML設定でvelocileptorsの使用を有効化/無効化可能にする

### 2.2 実装箇所

#### A. `_cosmolike_prototype_base.py`への追加

**1. 初期化処理（`initialize()`メソッド内）**
```python
# velocileptors関連の設定を読み込み
self.use_velocileptors = getattr(self, 'use_velocileptors', False)
if self.use_velocileptors:
    # velocileptorsのパス設定
    self.velocileptors_path = getattr(self, 'velocileptors_path', 
        './external_modules/code/velocileptors')
    # k-bin設定（YAMLから読み込み可能）
    self.velocileptors_k_min = getattr(self, 'velocileptors_k_min', 5e-3)
    self.velocileptors_k_max = getattr(self, 'velocileptors_k_max', 0.3)
    self.velocileptors_nk = getattr(self, 'velocileptors_nk', 50)
    # バイアスパラメータのデフォルト値
    self.velocileptors_bias_defaults = getattr(self, 'velocileptors_bias_defaults', {})
```

**2. velocileptors計算メソッドの追加**
```python
def compute_velocileptors_pk(self, z_eff, bias_params=None):
    """
    Velocileptorsを使って赤方偏移空間パワースペクトルを計算
    
    Parameters:
    -----------
    z_eff : float
        有効赤方偏移
    bias_params : dict, optional
        バイアスパラメータ（b1, b2, bs, b3, alpha0, alpha2, alpha4, SN0, SN2等）
    
    Returns:
    --------
    kv : array
        k値 [h/Mpc]
    p0, p2, p4 : arrays
        多極子（monopole, quadrupole, hexadecapole）
    """
    if not self.use_velocileptors:
        return None, None, None, None
    
    # velocileptorsをインポート
    import sys
    if self.velocileptors_path not in sys.path:
        sys.path.insert(0, self.velocileptors_path)
    from velocileptors.LPT.moment_expansion_fftw import MomentExpansion
    
    # MGCAMBから線形P(k)を取得
    h = self.provider.get_param("H0")/100.0
    PKL = self.provider.get_Pk_interpolator(("delta_tot", "delta_tot"),
        nonlinear=False, extrap_kmax=self.extrap_kmax)
    
    # k-gridを準備（velocileptors用）
    k_bins = np.logspace(np.log10(self.velocileptors_k_min), 
                         np.log10(self.velocileptors_k_max), 
                         self.velocileptors_nk)
    
    # 線形P(k)を抽出（z=z_effでの値）
    klin = k_bins  # [h/Mpc]
    plin = PKL.P(z_eff, klin) * (h**3)  # CosmoLike形式から変換
    
    # 成長率f(z)を計算（MGCAMBから取得または近似式）
    # 注意: MGCAMBから直接f(z)を取得できない場合は、近似式を使用
    # f(z) ≈ Ω_m(z)^0.55 または CAMBdataから取得
    
    # MomentExpansionを初期化
    mome = MomentExpansion(klin, plin, threads=1,
                          cutoff=10, extrap_min=-4, extrap_max=3, jn=10,
                          nk=self.velocileptors_nk, 
                          kmin=self.velocileptors_k_min, 
                          kmax=self.velocileptors_k_max)
    
    # バイアスパラメータを準備（デフォルト値とマージ）
    if bias_params is None:
        bias_params = {}
    defaults = {
        'b1': 2.0, 'b2': 0.0, 'bs': 0.0, 'b3': 0.0,
        'alpha0': 0.0, 'alpha2': 0.0, 'alpha4': 0.0,
        'SN0': 0.0, 'SN2': 0.0
    }
    defaults.update(self.velocileptors_bias_defaults)
    defaults.update(bias_params)
    
    # パラメータベクトルを構築（reduced basis）
    biases = [defaults['b1'], defaults['b2'], defaults['bs'], defaults['b3']]
    cterms = [defaults['alpha0'], defaults['alpha2'], defaults['alpha4']]
    stoch = [defaults['SN0'], defaults['SN2']]
    pars = biases + cterms + stoch
    
    # 成長率f(z)を取得（簡易版: 近似式を使用）
    # TODO: MGCAMBから正確なf(z)を取得する方法を実装
    omegam_z = self.provider.get_param("omegam") * ((1 + z_eff)**3) / \
               (self.provider.get_param("omegam") * ((1 + z_eff)**3) + 
                (1 - self.provider.get_param("omegam")))
    f_z = omegam_z**0.55  # 簡易近似
    
    # 多極子を計算
    kv, p0, p2, p4 = mome.compute_redshift_space_power_multipoles(
        pars, z_eff, f_z, ngauss=4, reduced=True)
    
    return kv, p0, p2, p4
```

**3. `set_cosmo_related()`への統合**
```python
def set_cosmo_related(self):
    # 既存の処理（MGCAMBからP(k)取得、CosmoLike C++に渡す）
    # ... (既存コード) ...
    
    # velocileptors計算（オプション）
    if self.use_velocileptors:
        # 有効赤方偏移を取得（データから、または設定から）
        z_eff = getattr(self, 'velocileptors_z_eff', 0.5)
        
        # バイアスパラメータを取得（YAMLから、またはデフォルト値）
        bias_params = {}
        # TODO: YAMLからバイアスパラメータを読み込む処理
        
        # velocileptorsで計算
        kv, p0, p2, p4 = self.compute_velocileptors_pk(z_eff, bias_params)
        
        # 結果を保存（後で使用するため）
        self.velocileptors_kv = kv
        self.velocileptors_p0 = p0
        self.velocileptors_p2 = p2
        self.velocileptors_p4 = p4
        
        # CosmoLike C++に渡す（必要に応じて）
        # TODO: CosmoLike C++側にvelocileptors結果を渡すインターフェースを実装
```

#### B. YAML設定ファイルへの追加

**`EXAMPLE_EVALUATE32_mg_hsc_y3.yaml`の例:**
```yaml
likelihood:
  mg_hsc_y3.mg_hsc_y3_2x2pt:
    path: ./external_modules/data/mg_hsc_y3
    data_file: mg_hsc_y3_M1_GGL0.05.dataset
    # ... 既存の設定 ...
    
    # velocileptors統合設定
    use_velocileptors: True
    velocileptors_path: ./external_modules/code/velocileptors
    velocileptors_k_min: 5e-3
    velocileptors_k_max: 0.3
    velocileptors_nk: 50
    velocileptors_z_eff: 0.5  # データの有効赤方偏移
    velocileptors_bias_defaults:
      b1: 2.0
      b2: 0.0
      bs: 0.0
      b3: 0.0
      alpha0: 0.0
      alpha2: 0.0
      alpha4: 0.0
      SN0: 0.0
      SN2: 0.0
```

#### C. 新しいLikelihoodクラスの作成（オプション）

velocileptors専用のLikelihoodクラスを作成する場合:
```python
# mg_hsc_y3_velocileptors.py
from cobaya.likelihoods.mg_hsc_y3._cosmolike_prototype_base import _cosmolike_prototype_base
import numpy as np

class mg_hsc_y3_velocileptors(_cosmolike_prototype_base):
    def initialize(self):
        super().initialize(probe="velocileptors")
        # velocileptors専用の初期化
    
    def internal_get_datavector(self, **params_values):
        self.set_cosmo_related()
        
        # velocileptors結果を使用してデータベクトルを構築
        kv = self.velocileptors_kv
        p0 = self.velocileptors_p0
        p2 = self.velocileptors_p2
        p4 = self.velocileptors_p4
        
        # データベクトルを構築（例: 多極子を結合）
        datavector = np.concatenate([p0, p2, p4])
        
        return datavector
```

### 2.3 実装ステップ

1. **Phase 1: 基本統合**
   - `_cosmolike_prototype_base.py`に`compute_velocileptors_pk()`メソッドを追加
   - YAML設定でvelocileptorsを有効化できるようにする
   - テスト: 単一のz_effでvelocileptorsが動作することを確認

2. **Phase 2: 成長率f(z)の正確な取得**
   - MGCAMBからf(z)を取得する方法を実装
   - または、CAMBdataから直接計算

3. **Phase 3: バイアスパラメータの統合**
   - YAMLからバイアスパラメータを読み込む
   - または、既存の`hsc_B1_*`パラメータをvelocileptorsにマッピング

4. **Phase 4: CosmoLike C++との統合（オプション）**
   - CosmoLike C++側にvelocileptors結果を渡すインターフェースを実装
   - または、Python側で直接データベクトルを構築

### 2.4 既存統計量への影響分析

#### ⚠️ 重要な確認事項

**1. 計算時間への影響**

**影響なし（`use_velocileptors=False`の場合）:**
- `use_velocileptors`フラグが`False`（デフォルト）の場合、velocileptors計算は**完全にスキップ**されます
- 既存の3x2pt、2x2pt、cosmic_shear等の計算時間には**一切影響しません**
- 初期化時のオーバーヘッドも最小限（フラグチェックのみ）

**影響あり（`use_velocileptors=True`の場合）:**
- `set_cosmo_related()`内でvelocileptors計算が実行される
- 計算時間: 約0.1-1秒程度（k-bin数、z_eff数に依存）
- **全てのLikelihoodクラスで`set_cosmo_related()`が呼ばれるため、3x2pt等でも計算される**
- ただし、結果は使用されない（無駄な計算になる）

**推奨対応:**
```python
# set_cosmo_related()内での実装
def set_cosmo_related(self):
    # 既存の処理（常に実行）
    # ... (既存コード) ...
    
    # velocileptors計算（条件付き）
    if getattr(self, 'use_velocileptors', False):
        # プローブがvelocileptorsを必要とする場合のみ計算
        probe = getattr(self, 'probe', '')
        if 'velocileptors' in probe.lower() or probe == 'kaiser_rsd':
            kv, p0, p2, p4 = self.compute_velocileptors_pk(z_eff, bias_params)
            # 結果を保存
```

**2. メモリ使用量への影響**

**影響なし（`use_velocileptors=False`の場合）:**
- 追加のメモリ使用なし

**影響あり（`use_velocileptors=True`の場合）:**
- MomentExpansionオブジェクト: 約10-50MB（k-bin数に依存）
- 計算結果（kv, p0, p2, p4）: 約1-5MB
- **合計: 約10-55MBの追加メモリ使用**

**3. 既存のCosmoLike計算への影響**

**影響なし:**
- velocileptors計算は既存のCosmoLike C++計算とは**独立**
- `ci.init_linear_power_spectrum()`等の既存処理は変更なし
- 既存のキャッシュ機構も影響なし

**4. キャッシュ機構への影響**

**現状の問題:**
- `set_cache_alert()`は`chi`, `lnPL`, `lnPNL`, `cosmo`のみをチェック
- velocileptors結果はキャッシュされない（毎回再計算される可能性）

**推奨対応:**
```python
def set_cache_alert(self, chi, lnPL, lnPNL):
    # 既存のチェック
    cache_alert_1 = np.array_equal(self.do_cache_chi, chi)
    cache_alert_2 = np.array_equal(self.do_cache_lnPL, lnPL)
    cache_alert_3 = np.array_equal(self.do_cache_lnPNL, lnPNL)
    cache_alert_4 = np.array_equal(
        self.do_cache_cosmo,
        np.array([
            self.provider.get_param("omegam"),
            self.provider.get_param("H0")
        ])
    )
    
    # velocileptors結果のキャッシュチェック（オプション）
    cache_alert_5 = True
    if getattr(self, 'use_velocileptors', False):
        if hasattr(self, 'do_cache_velocileptors'):
            # velocileptors結果が変わっていないかチェック
            # （バイアスパラメータ等も含める必要がある）
            pass  # TODO: 実装
    
    return (cache_alert_1 and cache_alert_2 and cache_alert_3 and 
            cache_alert_4 and cache_alert_5 and not self.force_cache_false)
```

**5. エラーハンドリング**

**推奨実装:**
```python
def compute_velocileptors_pk(self, z_eff, bias_params=None):
    """VelocileptorsでP(k)多極子を計算"""
    if not getattr(self, 'use_velocileptors', False):
        return None, None, None, None
    
    try:
        import sys
        velocileptors_path = getattr(self, 'velocileptors_path', 
            './external_modules/code/velocileptors')
        if velocileptors_path not in sys.path:
            sys.path.insert(0, velocileptors_path)
        
        from velocileptors.LPT.moment_expansion_fftw import MomentExpansion
        # ... 計算処理 ...
        
    except ImportError as e:
        self.log.warning(f"Velocileptors import failed: {e}. Skipping velocileptors calculation.")
        return None, None, None, None
    except Exception as e:
        self.log.warning(f"Velocileptors calculation failed: {e}. Skipping velocileptors calculation.")
        return None, None, None, None
```

**6. プローブ別の最適化**

**推奨: プローブがvelocileptorsを必要とする場合のみ計算**
```python
def set_cosmo_related(self):
    # 既存の処理（常に実行）
    # ... (既存コード) ...
    
    # velocileptors計算（条件付き）
    if getattr(self, 'use_velocileptors', False):
        probe = getattr(self, 'probe', '')
        # velocileptorsを必要とするプローブのみ計算
        velocileptors_probes = ['velocileptors', 'kaiser_rsd', 'pk_multipoles']
        if any(vp in probe.lower() for vp in velocileptors_probes):
            kv, p0, p2, p4 = self.compute_velocileptors_pk(z_eff, bias_params)
            # 結果を保存
        else:
            # 3x2pt, 2x2pt等では計算しない（無駄を避ける）
            pass
```

### 2.5 注意事項

1. **単位の統一**
   - MGCAMB: k in [1/Mpc], P(k) in [Mpc^3]
   - Velocileptors: k in [h/Mpc], P(k) in [(Mpc/h)^3]
   - CosmoLike: k in [h/Mpc], P(k) in [(Mpc/h)^3]
   - 単位変換に注意

2. **キャッシュ管理**
   - velocileptors計算結果もキャッシュに含める
   - `set_cache_alert()`を拡張してvelocileptors結果もチェック
   - バイアスパラメータの変更もキャッシュ無効化の条件に含める

3. **エラーハンドリング**
   - velocileptorsのインポート失敗時は警告を出してスキップ
   - 計算失敗時も既存処理を継続（フォールバック）

4. **パフォーマンス最適化**
   - velocileptors計算は比較的重いため、キャッシュを活用
   - プローブがvelocileptorsを必要としない場合は計算をスキップ
   - 必要に応じて並列化（threadsパラメータ）

5. **後方互換性**
   - `use_velocileptors`のデフォルトは`False`
   - 既存のYAMLファイルは変更不要
   - 既存のLikelihoodクラスの動作は変更なし

---

## 3. 実装例（最小限）

最小限の実装例を示します:

```python
# _cosmolike_prototype_base.py への追加

def compute_velocileptors_pk(self, z_eff, bias_params=None):
    """VelocileptorsでP(k)多極子を計算"""
    if not getattr(self, 'use_velocileptors', False):
        return None, None, None, None
    
    import sys
    velocileptors_path = getattr(self, 'velocileptors_path', 
        './external_modules/code/velocileptors')
    if velocileptors_path not in sys.path:
        sys.path.insert(0, velocileptors_path)
    
    from velocileptors.LPT.moment_expansion_fftw import MomentExpansion
    
    # MGCAMBから線形P(k)を取得
    h = self.provider.get_param("H0")/100.0
    PKL = self.provider.get_Pk_interpolator(("delta_tot", "delta_tot"),
        nonlinear=False, extrap_kmax=self.extrap_kmax)
    
    k_min = getattr(self, 'velocileptors_k_min', 5e-3)
    k_max = getattr(self, 'velocileptors_k_max', 0.3)
    nk = getattr(self, 'velocileptors_nk', 50)
    k_bins = np.logspace(np.log10(k_min), np.log10(k_max), nk)
    
    klin = k_bins
    plin = PKL.P(z_eff, klin) * (h**3)
    
    mome = MomentExpansion(klin, plin, threads=1, nk=nk, 
                          kmin=k_min, kmax=k_max)
    
    if bias_params is None:
        bias_params = {'b1': 2.0, 'b2': 0.0, 'bs': 0.0, 'b3': 0.0,
                      'alpha0': 0.0, 'alpha2': 0.0, 'alpha4': 0.0,
                      'SN0': 0.0, 'SN2': 0.0}
    
    pars = [bias_params['b1'], bias_params['b2'], bias_params['bs'], bias_params['b3'],
            bias_params['alpha0'], bias_params['alpha2'], bias_params['alpha4'],
            bias_params['SN0'], bias_params['SN2']]
    
    omegam_z = self.provider.get_param("omegam") * ((1 + z_eff)**3) / \
               (self.provider.get_param("omegam") * ((1 + z_eff)**3) + 
                (1 - self.provider.get_param("omegam")))
    f_z = omegam_z**0.55
    
    kv, p0, p2, p4 = mome.compute_redshift_space_power_multipoles(
        pars, z_eff, f_z, ngauss=4, reduced=True)
    
    return kv, p0, p2, p4
```

---

## 4. 次のステップ

1. 最小限の実装を`_cosmolike_prototype_base.py`に追加
2. YAML設定ファイルでテスト
3. 必要に応じて機能を拡張（f(z)の正確な取得、バイアスパラメータ統合等）

