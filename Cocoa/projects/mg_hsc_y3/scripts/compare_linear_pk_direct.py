#!/usr/bin/env python3
"""
Kaiser RSDとVelocileptorsが使用している線形P(k)を直接比較
"""
import numpy as np

print("=" * 80)
print("線形P(k)の直接比較")
print("=" * 80)

# Kaiser RSDの結果から線形P(k)を逆算
kaiser_file = './projects/mg_hsc_y3/data/mg_hsc_y3_theory_kaiser_rsd_boss_dr12.modelvector'
kaiser_data = np.loadtxt(kaiser_file)
values = kaiser_data[:, 1]

Nk = 40
Nz = 1
total_size = Nk * Nz

P0_kaiser = values[:total_size]
P2_kaiser = values[total_size:2*total_size]

# パラメータ
b1 = 1.2338
f = 0.790178
beta = f / (1 + b1)

# Kaiser係数
P0_coeff = (1 + 2/3*beta + 1/5*beta**2) * (1+b1)**2
P2_coeff = (4/3*beta + 4/7*beta**2) * (1+b1)**2

# 線形P(k)を逆算
Plin_from_P0 = P0_kaiser / P0_coeff
Plin_from_P2 = P2_kaiser / P2_coeff

print(f"\nパラメータ:")
print(f"   b1 = {b1}")
print(f"   f = {f}")
print(f"   β = f/(1+b1) = {beta:.4f}")
print(f"   (1+b1)² = {(1+b1)**2:.4f}")

print(f"\nKaiser係数:")
print(f"   P0係数 = {P0_coeff:.4f}")
print(f"   P2係数 = {P2_coeff:.4f}")

print(f"\n📊 Kaiser RSDから逆算した線形P(k):")
print(f"   Plin[0] (from P0) = {Plin_from_P0[0]:.2f}")
print(f"   Plin[0] (from P2) = {Plin_from_P2[0]:.2f}")
print(f"   比率 (P2/P0) = {Plin_from_P2[0] / Plin_from_P0[0]:.4f}")

# 最初の10点を表示
print(f"\n最初の10点:")
print(f"{'Index':<8} {'Plin(P0)':<15} {'Plin(P2)':<15} {'比率':<10}")
print("-" * 60)
for i in range(min(10, len(Plin_from_P0))):
    ratio = Plin_from_P2[i] / Plin_from_P0[i] if Plin_from_P0[i] != 0 else 0
    print(f"{i:<8} {Plin_from_P0[i]:<15.2f} {Plin_from_P2[i]:<15.2f} {ratio:<10.4f}")

# Velocileptors
velo_data = np.loadtxt('./projects/mg_hsc_y3/data/velocileptors_pk_multipoles.txt')
k_velo = velo_data[:, 0]
P0_velo = velo_data[:, 1]
P2_velo = velo_data[:, 2]

Plin_velo_from_P0 = P0_velo / P0_coeff
Plin_velo_from_P2 = P2_velo / P2_coeff

print(f"\n📊 Velocileptorsから逆算した線形P(k):")
print(f"   Plin[0] (from P0) = {Plin_velo_from_P0[0]:.2f}")
print(f"   Plin[0] (from P2) = {Plin_velo_from_P2[0]:.2f}")
print(f"   比率 (P2/P0) = {Plin_velo_from_P2[0] / Plin_velo_from_P0[0]:.4f}")

print("\n" + "=" * 80)
print("🎯 結論")
print("=" * 80)

plin_kaiser = Plin_from_P0[0]
plin_velo = Plin_velo_from_P0[0]
ratio = plin_velo / plin_kaiser
diff_pct = abs(plin_velo - plin_kaiser) / plin_kaiser * 100

print(f"""
Kaiser RSDが使っている線形P(k):  {plin_kaiser:.2f}
Velocileptorsが使っている線形P(k): {plin_velo:.2f}
比率: {ratio:.4f}
差: {diff_pct:.1f}%

**線形P(k)自体が{diff_pct:.1f}%違います！**

これは以下のいずれかが原因です：
1. Kaiser RSDとVelocileptorsで異なる線形P(k)を使っている
2. Kaiser RSDのmatter_pk(k,z)の取得方法に問題がある
3. Velocileptorsの線形P(k)の取得方法に問題がある
4. k-binsの位置が異なる（補間誤差）

次のステップ：
- Kaiser RSDのmatter_pk()関数が正しいP(k)を返しているか確認
- VelocileptorsのPKL.P(z,k)と同じ値か確認
- k-binsが一致しているか確認
""")












