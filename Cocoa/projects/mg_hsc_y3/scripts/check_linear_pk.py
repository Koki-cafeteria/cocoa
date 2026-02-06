#!/usr/bin/env python3
"""
KaiserとVelociléptorsで使用している線形P(k)の値を確認
"""
import numpy as np
import sys
sys.path.insert(0, '/home/tanida/cocoa_v41/cocoa/Cocoa')

# CAMBを初期化して線形P(k)を取得
from cobaya.model import get_model

# Kaiser RSD YAMLをロード
kaiser_info = "/home/tanida/cocoa_v41/cocoa/Cocoa/projects/mg_hsc_y3/EXAMPLE_EVALUATE_kaiser_RSD.yaml"

print("=" * 80)
print("線形パワースペクトル Pₘ(k) の値を確認")
print("=" * 80)

print("\n📊 CAMBからP(k)を取得しています...")
model = get_model(kaiser_info)

# パラメータを取得
params = {
    'H0': 67.089,
    'omegab': 0.049434,
    'omegam': 0.3156,
    'As_1e9': 2.2065,
    'ns': 0.9645,
}

# CAMBを実行
point = model.logposterior(params)

# Providerから線形P(k)を取得
provider = model.provider
h = 67.089 / 100.0
z_test = 0.5893

# テスト用のk値（BOSS DR12の最初の数点）
k_test = np.array([0.007424, 0.010, 0.020, 0.030, 0.050])

print(f"\n📊 z = {z_test}, h = {h:.4f}")
print(f"\n{'k [h/Mpc]':<15} {'k [1/Mpc]':<15} {'Pₘ(k) [Mpc³]':<20} {'Pₘ(k) [(Mpc/h)³]':<20}")
print("-" * 80)

PKL = provider.get_Pk_interpolator(("delta_tot", "delta_tot"), nonlinear=False)

for k_h in k_test:
    k_mpc = k_h / h  # h/Mpc → 1/Mpc
    pk_mpc3 = PKL.P(z_test, k_mpc)  # Mpc³
    pk_mpch3 = pk_mpc3 * (h**3)  # (Mpc/h)³
    print(f"{k_h:<15.6f} {k_mpc:<15.6f} {pk_mpc3:<20.2f} {pk_mpch3:<20.2f}")

print("\n" + "=" * 80)
print("Kaiser RSD の期待値")
print("=" * 80)

b = 1.234
f = 0.790178
beta = f / b
kaiser_coeff_P0 = 1 + (2/3)*beta + (1/5)*beta**2
overall_coeff = b**2 * kaiser_coeff_P0

print(f"\nb = {b}, f = {f:.6f}, β = {beta:.4f}")
print(f"Kaiser係数 (P₀) = {overall_coeff:.4f}")
print(f"\n{'k [h/Mpc]':<15} {'Pₘ [(Mpc/h)³]':<20} {'P₀予測 [(Mpc/h)³]':<25}")
print("-" * 80)

for k_h in k_test:
    k_mpc = k_h / h
    pk_mpc3 = PKL.P(z_test, k_mpc)
    pk_mpch3 = pk_mpc3 * (h**3)
    p0_expected = overall_coeff * pk_mpch3
    print(f"{k_h:<15.6f} {pk_mpch3:<20.2f} {p0_expected:<25.2f}")

print("\n✅ この値とKaiser RSDの実際の出力を比較してください！")












