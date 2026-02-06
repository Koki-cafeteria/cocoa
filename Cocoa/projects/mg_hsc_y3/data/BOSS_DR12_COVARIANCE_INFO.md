# BOSS DR12 共分散行列データの取得方法

## データリリースの場所

BOSS DR12のパワースペクトルデータと共分散行列は、SDSS Data Release 12の公式サイトから取得できます：

- **SDSS DR12 ホームページ**: https://www.sdss.org/dr12/
- **データアクセスポートal**: https://www.sdss.org/dr12/data_access/

## 共分散行列ファイルの場所

BOSS DR12のパワースペクトルデータは通常、以下のような構造になっています：

```
BOSS_DR12_PowerSpectrum/
├── Power_Spectrum_cmass_ngc_v5.txt      # NGCサンプルのパワースペクトル
├── Power_Spectrum_cmass_sgc_v5.txt      # SGCサンプルのパワースペクトル
├── Covariance_Matrix_cmass_ngc_v5.txt   # NGCサンプルの共分散行列
├── Covariance_Matrix_cmass_sgc_v5.txt   # SGCサンプルの共分散行列
└── ...
```

## ファイル名の例

- `Covariance_Matrix_cmass_ngc_v5.txt` - CMASS NGCサンプル用
- `Covariance_Matrix_cmass_sgc_v5.txt` - CMASS SGCサンプル用
- `Covariance_Matrix_cmass_combined_v5.txt` - CMASS結合サンプル用

## 共分散行列の形式

BOSS DR12の共分散行列は通常、以下のいずれかの形式です：

1. **モノポール+クアドルポール形式** (2*nk x 2*nk)
   - 順序: [P0(k1), ..., P0(nk), P2(k1), ..., P2(nk)]
   - この場合は `--format mono_quad` を使用

2. **モノポール+クアドルポール+ヘキサデカポール形式** (3*nk x 3*nk)
   - 順序: [P0(k1), ..., P0(nk), P2(k1), ..., P2(nk), P4(k1), ..., P4(nk)]
   - この場合は `--format full` を使用

## ダウンロードと変換の手順

1. **共分散行列ファイルをダウンロード**
   ```bash
   # 例: wgetやcurlを使用
   wget https://[BOSS_DR12_URL]/Covariance_Matrix_cmass_ngc_v5.txt
   ```

2. **ファイルを変換**
   ```bash
   cd /home/tanida/cocoa_v41/cocoa/Cocoa/projects/mg_hsc_y3/data
   
   # モノポール+クアドルポール形式の場合
   python convert_boss_dr12_covariance.py \
       Covariance_Matrix_cmass_ngc_v5.txt \
       boss_dr12_cmass_covariance.txt \
       --nk 40 \
       --format mono_quad
   
   # または、モノポール+クアドルポール+ヘキサデカポール形式の場合
   python convert_boss_dr12_covariance.py \
       Covariance_Matrix_cmass_ngc_v5.txt \
       boss_dr12_cmass_covariance.txt \
       --nk 40 \
       --format full
   ```

3. **`.dataset`ファイルの確認**
   - `boss_dr12_cmass_kaiser_rsd.dataset`で`cov_file`が正しく設定されているか確認

## 参考リンク

- SDSS DR12: https://www.sdss.org/dr12/
- BOSS DR12 Power Spectrum論文: 
  - Beutler et al. (2017) "The clustering of galaxies in the completed SDSS-III Baryon Oscillation Spectroscopic Survey: Baryon Acoustic Oscillations in the Fourier space"
  - データファイルの正確な場所は論文の補足資料に記載されている場合があります

## トラブルシューティング

### 共分散行列のサイズが合わない場合

ファイルの最初の数行を確認して、実際の形式を判断してください：

```bash
head -20 Covariance_Matrix_cmass_ngc_v5.txt
```

- コメント行（#で始まる）の数を確認
- データ行の数を確認
- 1行あたりの列数を確認

### 形式が不明な場合

`convert_boss_dr12_covariance.py`スクリプトは自動的に形式を検出しようとしますが、手動で確認する場合は：

```python
import numpy as np
data = np.loadtxt('Covariance_Matrix_cmass_ngc_v5.txt', comments='#')
print(f"Matrix shape: {data.shape}")
print(f"Expected for mono_quad: (80, 80)")
print(f"Expected for full: (120, 120)")
```

