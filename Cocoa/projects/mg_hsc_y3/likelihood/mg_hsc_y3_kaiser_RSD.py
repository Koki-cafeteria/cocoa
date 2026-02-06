# from cobaya.likelihoods.mg_hsc_y3._cosmolike_prototype_base import _cosmolike_prototype_base
# import cosmolike_mg_hsc_y3_interface as ci
# import numpy as np
# from getdist import IniFile
# import os

# class mg_hsc_y3_kaiser_RSD(_cosmolike_prototype_base):
#   # 250801: initializeメソッドを修正
#   def initialize(self):
#     # 3x2ptとして基底クラスを初期化し、共通処理を任せる
#     super(mg_hsc_y3_kaiser_RSD, self).initialize(probe="3x2pt")
    
#     # .dataset ファイルを再度読み込み、kaiser_RSD用の設定を取得
#     dataset_file_path = os.path.join(self.path, self.data_file)
#     ini = IniFile(dataset_file_path)
    
#     # P(k,z)のビニングを設定
#     self.Nk = ini.int('Nk')
#     self.k_min = ini.float('k_min')
#     self.k_max = ini.float('k_max')
#     self.Nz = ini.int('Nz')
#     self.z_min = ini.float('z_min')
#     self.z_max = ini.float('z_max')
    
#     # C++側にk, zのビニング情報を渡す
#     ci.init_binning_kaiser_RSD(
#         Nk=self.Nk, k_min=self.k_min, k_max=self.k_max,
#         Nz=self.Nz, z_min=self.z_min, z_max=self.z_max
#     )
    
#     # 統合されたデータベクトルのサイズを計算
#     ci.init_data_vector_size_3x2pt_and_kaiser_RSD()

#   # logp, get_datavectorは変更なし
#   def logp(self, **params_values):
#     datavector = self.internal_get_datavector(**params_values)
#     return self.compute_logp(datavector)

#   def get_datavector(self, **params_values):        
#     datavector = self.internal_get_datavector(**params_values)
#     return np.array(datavector)

#   # internal_get_datavectorは、新しい統合関数を呼び出すように修正
#   def internal_get_datavector(self, **params_values):
#     self.set_cosmo_related()
#     self.set_lens_related(**params_values)
#     self.set_source_related(**params_values)
    
#     # 新しく作成した統合データベクトル計算関数を呼び出す
#     datavector = np.array(ci.compute_data_vector_3x2pt_and_kaiser_RSD_masked())

#     if self.print_datavector:
#       size = len(datavector)
#       out = np.zeros(shape=(size, 2))
#       out[:,0] = np.arange(0, size)
#       out[:,1] = datavector
#       fmt = '%d', '%1.8e'
#       np.savetxt(self.print_datavector_file, out, fmt = fmt)
                
#     return datavector
import os
import numpy as np
from getdist import IniFile
import cosmolike_mg_hsc_y3_interface as ci
from cobaya.likelihoods.mg_hsc_y3._cosmolike_prototype_base import _cosmolike_prototype_base

class mg_hsc_y3_kaiser_RSD(_cosmolike_prototype_base):

  def initialize(self):
    # .dataset ファイルから情報を読み込む
    dataset_file_path = os.path.join(self.path, self.data_file)
    ini = IniFile(dataset_file_path)
    
    # --- 1. C++側の状態をセットアップ ---
    ci.initial_setup()
    # Kaiser単独モード: 3x2ptは無効化
    ci.init_probes(possible_probes="kaiser_rsd")

    # IAモデルの設定（Kaiser単独では不要だが、エラー回避のため設定）
    self.IA_model = ini.int('IA_model', 0)
    self.IA_redshift_evolution = ini.int('IA_redshift_evolution', 0)
    ci.init_IA(ia_model=self.IA_model, ia_redshift_evolution=self.IA_redshift_evolution)

    # Kaiser単独ではsource sampleも初期化が必要（init_lens_sampleのチェック回避のため）
    ci.init_source_sample(
        filename=ini.relativeFileName('nz_source_file'),
        ntomo_bins=ini.int('source_ntomo'))
    
    # Kaiser単独ではレンズサンプルのみ必要（バイアス用）
    ci.init_lens_sample(
        filename=ini.relativeFileName('nz_lens_file'),
        ntomo_bins=ini.int('lens_ntomo'))
    
    # Kaiser用のk, zビニングを設定（k_bins_fileがあれば優先）
    self.Nz = ini.int('Nz')
    if ini.hasKey('k_bins_file'):
      k_bins_file = ini.relativeFileName('k_bins_file')
      k_bins = np.loadtxt(k_bins_file)
      if k_bins.ndim == 2:
        # 2列以上ある場合は最後の列をk配列として採用（例: [index, k_eff]）
        k_bins = k_bins[:, -1]
      elif k_bins.ndim != 1:
        raise RuntimeError("k_bins_file must contain a 1D array or NxM with k in last column")
      self.Nk = len(k_bins)
      ci.init_binning_kaiser_RSD_custom_kbins(
          k_bins_array=k_bins, Nz=self.Nz,
          z_min=ini.float('z_min'), z_max=ini.float('z_max'))
      print(f"[mg_hsc_y3_kaiser_RSD] Using custom k-bins from {k_bins_file} (Nk={self.Nk})")
    else:
      self.Nk = ini.int('Nk')
      ci.init_binning_kaiser_RSD(
          Nk=self.Nk, k_min=ini.float('k_min'), k_max=ini.float('k_max'),
          Nz=self.Nz, z_min=ini.float('z_min'), z_max=ini.float('z_max'))
    
    # --- 2. Kaiser単独のデータサイズを計算 ---
    ci.init_data_vector_size_kaiser_RSD()
    
    # --- 3. データ・共分散・マスクをロード ---
    # 注意: Kaiser単独用のデータファイルを使用
    self.data_file_name = ini.relativeFileName('data_file')
    self.cov_file = ini.relativeFileName('cov_file')
    self.mask_file = ini.relativeFileName('mask_file')
    
    # Kaiser単独用のデータロード関数を使用
    ci.init_data_kaiser_RSD(
        COV=self.cov_file, MASK=self.mask_file, DATA=self.data_file_name
    )
    
    # --- 4. Cobayaが必要とする属性を設定 ---
    self.data_vector = np.loadtxt(self.data_file_name)[:, 1]
    self.inv_cov = ci.get_inv_cov_masked()
    self.lens_ntomo = ini.int('lens_ntomo')
    # Kaiser単独ではsource_ntomoは不要だが、エラー回避のため設定
    self.source_ntomo = ini.int('source_ntomo', 1)
    # 基底クラスが持つ他の属性も念のため設定しておく
    self.use_baryon_pca = getattr(self, "use_baryon_pca", False)
    self.create_baryon_pca = getattr(self, "create_baryon_pca", False)
    
    # get_requirements()で必要な属性を設定
    self.z_interp_1D = np.linspace(0, 2.0, 1000)
    self.z_interp_1D = np.concatenate((self.z_interp_1D, np.linspace(2.0, 10.1, 200)), axis=0)
    self.z_interp_1D[0] = 0
    
    self.z_interp_2D = np.linspace(0, 2.0, 95)
    self.z_interp_2D = np.concatenate((self.z_interp_2D, np.linspace(2.0, 10, 5)), axis=0)
    self.z_interp_2D[0] = 0
    
    self.len_z_interp_2D = len(self.z_interp_2D)
    self.len_log10k_interp_2D = 1200
    self.log10k_interp_2D = np.linspace(-4.2, 2.0, self.len_log10k_interp_2D)
    self.k_interp_2D = np.power(10.0, self.log10k_interp_2D)
    self.len_k_interp_2D = len(self.k_interp_2D)
    self.len_pkz_interp_2D = self.len_log10k_interp_2D * self.len_z_interp_2D
    self.extrap_kmax = 2.5e2 * self.accuracyboost
    
    # キャッシュ用の配列
    self.do_cache_lnPL = np.zeros(self.len_pkz_interp_2D)
    self.do_cache_lnPNL = np.zeros(self.len_pkz_interp_2D)
    self.do_cache_chi = np.zeros(len(self.z_interp_1D))
    self.do_cache_cosmo = np.zeros(2)

  def get_requirements(self):
    """
    Override to add Kaiser RSD bias parameters
    
    YAMLで以下を設定してMG_paramsを制御：
    - use_mg_params: True  → MG_paramsを要求（修正重力理論を計算）
    - use_mg_params: False → MG_paramsを不要（標準宇宙論のみ）
    """
    req = super().get_requirements()
    
    # MG_paramsは条件付きで削除
    # use_mg_paramsフラグでcontrol（デフォルト: False = 不要）
    use_mg_params = getattr(self, 'use_mg_params', False)
    if not use_mg_params:
      req.pop("MG_params", None)
    
    # Add Kaiser RSD bias parameters to requirements
    Nz = getattr(self, 'Nz', 10)
    for i in range(Nz):
      req[f"hsc_B_KAISER_RSD_{i+1}"] = None
    return req
  
  def logp(self, **params_values):
    datavector = self.internal_get_datavector(**params_values)
    return self.compute_logp(datavector)

  def get_datavector(self, **params_values):        
    datavector = self.internal_get_datavector(**params_values)
    return np.array(datavector)

  # Kaiser単独用のデータベクトル計算
  def internal_get_datavector(self, **params_values):
    # ▼▼▼ MGパラメータとpower spectrumを設定（3x2ptと同様） ▼▼▼
    
    # MGパラメータを設定
    try:
      mg_params = self.provider.get_MG_params()
      z_table = mg_params["z"]
      mu_table = mg_params["mu"]
      Sigma_table = mg_params["Sigma"]
      ci.set_mg_params(z_table, mu_table, Sigma_table)
    except:
      pass
    
    h = self.provider.get_param("H0")/100.0
    
    # CAMBから線形・非線形パワースペクトルを取得
    PKL = self.provider.get_Pk_interpolator(("delta_tot", "delta_tot"),
      nonlinear=False, extrap_kmax=self.extrap_kmax)
    PKNL = self.provider.get_Pk_interpolator(("delta_tot", "delta_tot"),
      nonlinear=True, extrap_kmax=self.extrap_kmax)
    
    # パワースペクトルを準備
    lnPL = np.empty(self.len_pkz_interp_2D)
    lnPNL = np.empty(self.len_pkz_interp_2D)
    
    t1 = PKNL.logP(self.z_interp_2D, self.k_interp_2D).flatten()
    t2 = PKL.logP(self.z_interp_2D, self.k_interp_2D).flatten()
    
    # Cosmolikeはk in h/Mpcを要求
    log10k_interp_2D = self.log10k_interp_2D - np.log10(h)
    
    for i in range(self.len_z_interp_2D):
      lnPL[i::self.len_z_interp_2D] = t2[i*self.len_k_interp_2D:(i+1)*self.len_k_interp_2D]
    lnPL += np.log((h**3))
    
    for i in range(self.len_z_interp_2D):
      lnPNL[i::self.len_z_interp_2D] = t1[i*self.len_k_interp_2D:(i+1)*self.len_k_interp_2D]  
    lnPNL += np.log((h**3))
    
    # 共動距離を計算 - Mpc/hに変換
    chi = self.provider.get_comoving_radial_distance(self.z_interp_1D) * h
    
    # 宇宙論パラメータを設定
    ci.set_cosmological_parameters(
      omega_matter = self.provider.get_param("omegam"),
      hubble = self.provider.get_param("H0"),
      is_cached = False
    )
    
    # パワースペクトルをcosmolikeに渡す
    ci.init_linear_power_spectrum(log10k=log10k_interp_2D,
      z=self.z_interp_2D, lnP=lnPL)
    ci.init_non_linear_power_spectrum(log10k=log10k_interp_2D,
      z=self.z_interp_2D, lnP=lnPNL)
    
    # 成長関数を計算・設定
    G_growth = np.sqrt(PKL.P(self.z_interp_2D, 0.0005) / PKL.P(0, 0.0005))
    G_growth = G_growth * (1 + self.z_interp_2D) / G_growth[len(G_growth)-1]
    ci.init_growth(z=self.z_interp_2D, G=G_growth)
    
    # デバッグ：f(z)の値を確認（Velocileptorsと比較するため）
    # Kaiser RSDで使用するz値（通常は最初のビンの有効赤方偏移）
    try:
      from scipy.interpolate import interp1d
      # like.z_binsから有効赤方偏移を取得（最初のビン）
      z_kaiser = getattr(self, 'z_min', 0.5893)  # デフォルトは0.5893
      # MGCAMBから直接f(z)を取得して比較
      fsigma8_array = self.provider.get_fsigma8(self.z_interp_2D)
      sigma8_array = self.provider.get_sigma8_z(self.z_interp_2D)
      fsigma8_z = interp1d(self.z_interp_2D, fsigma8_array, kind='linear', 
                          fill_value='extrapolate', assume_sorted=True)(z_kaiser)
      sigma8_z = interp1d(self.z_interp_2D, sigma8_array, kind='linear', 
                         fill_value='extrapolate', assume_sorted=True)(z_kaiser)
      f_z_mgcamb = fsigma8_z / sigma8_z
      print(f'[mg_hsc_y3.mg_hsc_y3_kaiser_RSD] Using exact f(z={z_kaiser:.3f}) = {f_z_mgcamb:.6f} from MGCAMB (fsigma8={fsigma8_z:.6f}, sigma8={sigma8_z:.6f})')
    except Exception as e:
      print(f'[mg_hsc_y3.mg_hsc_y3_kaiser_RSD] Could not compute f(z): {e}')
    
    # 距離を設定
    ci.init_distances(z=self.z_interp_1D, chi=chi)
    
    # ▲▲▲ ここまで ▲▲▲
    
    # Kaiser単独ではレンズバイアスのみ必要（3x2pt用）
    self.set_lens_related(**params_values)
    # source_relatedは不要（Kaiser単独ではIAやshear calibrationは使わない）
    
    # Kaiser RSD用のバイアスを設定（3x2ptのバイアスとは別）
    # Nzビン分のバイアスパラメータを取得
    Nz = getattr(self, 'Nz', 10)  # デフォルトは10
    kaiser_bias = [
      params_values.get(f"hsc_B_KAISER_RSD_{i+1}", 1.0) 
      for i in range(Nz)
    ]
    ci.set_nuisance_kaiser_RSD_bias(B_KAISER_RSD=kaiser_bias)
    
    # Kaiser単独用のデータベクトル計算関数を呼び出す
    datavector = np.array(ci.compute_data_vector_kaiser_RSD_masked())

    # print処理
    if self.print_datavector:
      np.savetxt(self.print_datavector_file, np.column_stack((np.arange(len(datavector)), datavector)), fmt=['%d', '%1.8e'])
                
    return datavector