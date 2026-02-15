from cobaya.likelihoods.mg_hsc_y3._cosmolike_prototype_base import _cosmolike_prototype_base
import numpy as np
import os
import cosmolike_mg_hsc_y3_interface as ci
from getdist import IniFile
from scipy.interpolate import interp1d

class mg_hsc_y3_velocileptors(_cosmolike_prototype_base):
  """
  velocileptorsによるパワースペクトル多極モーメント（P0, P2, P4）を計算する独立したlikelihood
  
  2x2ptや3x2ptとは独立して、velocileptorsだけを計算できる。
  YAMLの設定によって、velocileptors単独、または他の統計量と同時に計算可能。
  """
  
  # Cobayaが認識できるようにクラス変数として宣言（値は設定しない）
  velocileptors_k_bins_file: str
  
  def initialize(self):
    # velocileptors専用の初期化（C++側の初期化は最小限）
    # .datasetファイルから情報を読み込む
    dataset_file_path = os.path.join(self.path, self.data_file)
    ini = IniFile(dataset_file_path)
    
    # 基底クラスの属性を設定（C++側の初期化はスキップ）
    self.probe = "velocileptors"
    
    # velocileptors設定（.datasetファイルから読み込み、デフォルト値付き）
    self.use_velocileptors = getattr(self, 'use_velocileptors', False)
    if ini.hasKey('velocileptors_k_min'):
      self.velocileptors_k_min = ini.float('velocileptors_k_min')
    else:
      self.velocileptors_k_min = getattr(self, 'velocileptors_k_min', 5e-3)
    
    if ini.hasKey('velocileptors_k_max'):
      self.velocileptors_k_max = ini.float('velocileptors_k_max')
    else:
      self.velocileptors_k_max = getattr(self, 'velocileptors_k_max', 0.3)
    
    if ini.hasKey('velocileptors_nk'):
      self.velocileptors_nk = ini.int('velocileptors_nk')
    else:
      self.velocileptors_nk = getattr(self, 'velocileptors_nk', 50)
    
    # k_bins_fileが指定されていれば、公式k-binsを使用
    if ini.hasKey('velocileptors_k_bins_file'):
      self.velocileptors_k_bins_file = ini.string('velocileptors_k_bins_file')
    else:
      self.velocileptors_k_bins_file = getattr(self, 'velocileptors_k_bins_file', None)
    
    if ini.hasKey('velocileptors_z_eff'):
      self.velocileptors_z_eff = ini.float('velocileptors_z_eff')
    else:
      self.velocileptors_z_eff = getattr(self, 'velocileptors_z_eff', 0.5)
    
    # velocileptors_pathは.datasetファイルからは読み込まない（YAMLから指定）
    self.velocileptors_path = getattr(self, 'velocileptors_path', 
        './external_modules/code/velocileptors')
    
    # velocileptors_bias_defaultsも.datasetファイルからは読み込まない（YAMLから指定）
    self.velocileptors_bias_defaults = getattr(self, 'velocileptors_bias_defaults', {})
    
    # 基底クラスが必要とする属性を設定（set_cosmo_related用）
    self.accuracyboost = getattr(self, 'accuracyboost', 1.0)
    self.samplingboost = getattr(self, 'samplingboost', 1.0)
    self.integration_accuracy = getattr(self, 'integration_accuracy', 0)
    self.kmax_boltzmann = getattr(self, 'kmax_boltzmann', 5.0)
    self.non_linear_emul = getattr(self, 'non_linear_emul', 2)
    
    # z_interp_2Dとk_interp_2Dを設定（set_cosmo_related用）
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
    
    # z_interp_1Dを設定（set_cosmo_related用）
    self.z_interp_1D = np.linspace(0, 2.0, 1000)
    self.z_interp_1D = np.concatenate((self.z_interp_1D,
      np.linspace(2.0, 10.1, 200)), axis=0)
    self.z_interp_1D[0] = 0
    
    # キャッシュ用の配列
    self.do_cache_lnPL = np.zeros(self.len_pkz_interp_2D)
    self.do_cache_lnPNL = np.zeros(self.len_pkz_interp_2D)
    self.do_cache_chi = np.zeros(len(self.z_interp_1D))
    self.do_cache_cosmo = np.zeros(2)
    
    # C++側の初期化は最小限（set_cosmo_relatedで必要）
    ci.initial_setup()
    ci.init_accuracy_boost(self.accuracyboost, self.samplingboost, self.integration_accuracy)
    # 注意: init_probesは呼ばない（velocileptorsはC++側でサポートされていない）
    
    # データファイル、共分散行列、マスクファイルのパスを取得（オプション）
    # データファイルがない場合は理論値のみ計算
    self.has_data = ini.hasKey('velocileptors_data_file')
    if self.has_data:
      self.data_file_name = ini.relativeFileName('velocileptors_data_file')
      self.cov_file = ini.relativeFileName('velocileptors_cov_file')
      self.mask_file = ini.relativeFileName('velocileptors_mask_file') if ini.hasKey('velocileptors_mask_file') else None
      
      # データと共分散をロード
      self.data_vector = np.loadtxt(self.data_file_name)
      if self.data_vector.ndim == 2:
        # 2列以上ある場合は最後の列をデータとして採用
        self.data_vector = self.data_vector[:, -1]
      
      self.cov_matrix = np.loadtxt(self.cov_file)
      if self.cov_matrix.ndim == 1:
        # 対角共分散の場合
        self.cov_matrix = np.diag(self.cov_matrix)
      
      # マスクがある場合は適用
      if self.mask_file is not None:
        self.mask = np.loadtxt(self.mask_file, dtype=bool)
        self.data_vector = self.data_vector[self.mask]
        self.cov_matrix = self.cov_matrix[np.ix_(self.mask, self.mask)]
      else:
        self.mask = None
      
      # 逆共分散行列を計算
      self.inv_cov = np.linalg.inv(self.cov_matrix)
    else:
      # データファイルがない場合は理論値のみ計算
      self.data_vector = None
      self.cov_matrix = None
      self.inv_cov = None
      self.mask = None
      self.log.info("No data file specified. Will compute theory only.")
    
    # 出力設定
    self.print_velocileptors_pk = getattr(self, 'print_velocileptors_pk', True)
    self.velocileptors_pk_output_file = getattr(self, 'velocileptors_pk_output_file',
        './projects/mg_hsc_y3/data/velocileptors_pk_multipoles.txt')

  def get_requirements(self):
    """
    velocileptors専用のrequirements
    DESIスタイルのパラメータ化に対応: (1+b1)*sigma8, b2*sigma8^2, bs*sigma8^2
    
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
    
    # DESIスタイルのパラメータを要求
    req['velocileptors_b1plus1_sigma8'] = None  # (1+b1)*sigma8
    req['velocileptors_b2_sigma8_sq'] = None    # b2*sigma8^2
    req['velocileptors_bs_sigma8_sq'] = None    # bs*sigma8^2
    req['velocileptors_b3'] = None
    req['velocileptors_alpha0'] = None
    req['velocileptors_alpha2'] = None
    req['velocileptors_alpha4'] = None
    req['velocileptors_SN0'] = None
    req['velocileptors_SN2'] = None
    
    # sigma8が必要（DESIスタイルのパラメータ変換用）
    req['sigma8'] = None
    
    # fsigma8とsigma8_zは親クラスで要求されているので、ここでは追加不要
    
    return req

  def logp(self, **params_values):
    """
    Log-likelihoodを計算
    
    データファイルがある場合: chi2を計算
    データファイルがない場合: 0を返す（理論値のみ計算）
    """
    if not self.has_data:
      # データがない場合は理論値のみ計算（likelihoodは0）
      self.internal_get_datavector(**params_values)
      return 0.0
    
    datavector = self.internal_get_datavector(**params_values)
    return self.compute_logp(datavector)

  def get_datavector(self, **params_values):
    """
    データベクトル（理論値）を取得
    """
    datavector = self.internal_get_datavector(**params_values)
    return np.array(datavector)

  def internal_get_datavector(self, **params_values):
    """
    velocileptorsによるパワースペクトル多極モーメントを計算
    DESIスタイルのパラメータ化に対応
    
    Returns:
    --------
    datavector: array
        P0, P2, P4を結合した配列（kの順に並ぶ）
        [P0(k1), P2(k1), P4(k1), P0(k2), P2(k2), P4(k2), ...]
    """
    # velocileptors計算に必要な部分のみ実行（C++側の初期化は不要）
    # MGパラメータを設定（必要に応じて）
    try:
      mg_params = self.provider.get_MG_params()
      z_table = mg_params["z"]
      mu_table = mg_params["mu"]
      Sigma_table = mg_params["Sigma"]
      ci.set_mg_params(z_table, mu_table, Sigma_table)
    except:
      # MGパラメータが取得できない場合はスキップ（標準宇宙論の場合）
      pass
    
    # velocileptorsで直接計算（set_cosmo_relatedは呼ばない）
    z_eff = self.velocileptors_z_eff
    
    # z_effでのσ8を取得（compute_velocileptors_pkと同じ方法）
    # DESIスタイルのパラメータ変換に必要
    sigma8_array = self.provider.get_sigma8_z(self.z_interp_2D)
    sigma8 = interp1d(self.z_interp_2D, sigma8_array, 
                      kind='linear', fill_value='extrapolate', 
                      assume_sorted=True)(z_eff)
    
    # DESIスタイルのパラメータから実際のb1, b2, bsを計算
    bias_params = {}
    
    # (1+b1)*sigma8 から b1 を計算
    b1plus1_sigma8 = self.provider.get_param('velocileptors_b1plus1_sigma8')
    bias_params['b1'] = (b1plus1_sigma8 / sigma8) - 1.0
    
    # b2*sigma8^2 から b2 を計算
    b2_sigma8_sq = self.provider.get_param('velocileptors_b2_sigma8_sq')
    bias_params['b2'] = b2_sigma8_sq / (sigma8 ** 2)
    
    # bs*sigma8^2 から bs を計算
    bs_sigma8_sq = self.provider.get_param('velocileptors_bs_sigma8_sq')
    bias_params['bs'] = bs_sigma8_sq / (sigma8 ** 2)
    
    # その他のパラメータは直接取得
    for param_name in ['b3', 'alpha0', 'alpha2', 'alpha4', 'SN0', 'SN2']:
      param_key = f'velocileptors_{param_name}'
      try:
        bias_params[param_name] = self.provider.get_param(param_key)
      except:
        defaults = {'b3': 0.0, 'alpha0': 0.0, 'alpha2': 0.0, 'alpha4': 0.0, 'SN0': 0.0, 'SN2': 0.0}
        bias_params[param_name] = defaults.get(param_name, 0.0)
    
    kv, p0, p2, p4 = self.compute_velocileptors_pk(z_eff, bias_params)
    
    if kv is None:
      raise RuntimeError("velocileptors calculation failed. Check use_velocileptors=True and velocileptors parameters.")
    
    # データベクトルとして結合: [P0(k1), P2(k1), P4(k1), P0(k2), P2(k2), P4(k2), ...]
    datavector = np.zeros(len(kv) * 3)
    datavector[0::3] = p0  # P0
    datavector[1::3] = p2  # P2
    datavector[2::3] = p4  # P4
    
    # 出力処理
    if self.print_velocileptors_pk:
      output_data = np.column_stack([kv, p0, p2, p4])
      np.savetxt(self.velocileptors_pk_output_file, output_data,
                fmt='%.8e',
                header='k[h/Mpc]  P0  P2  P4',
                comments='# ')
      self.log.info(f'Velocileptors P(k) multipoles saved to {self.velocileptors_pk_output_file}')
    
    return datavector

  def compute_logp(self, datavector):
    """
    Log-likelihoodを計算（データがある場合のみ）
    
    Parameters:
    -----------
    datavector: array
        理論値のデータベクトル
    
    Returns:
    --------
    logp: float
        Log-likelihood (-0.5 * chi2)
    """
    if not self.has_data:
      return 0.0
    
    # 理論ベクトルにマスクを適用（データと同じサイズにする）
    if self.mask is not None:
      datavector = datavector[self.mask]
    
    # chi2 = (theory - data)^T * inv_cov * (theory - data)
    theory = datavector
    diff = theory - self.data_vector
    chi2 = np.dot(diff, np.dot(self.inv_cov, diff))
    
    return -0.5 * chi2
