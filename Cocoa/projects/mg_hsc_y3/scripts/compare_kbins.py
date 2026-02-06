#!/usr/bin/env python3
import numpy as np
import os

dataset = "/home/tanida/cocoa_v41/cocoa/Cocoa/projects/mg_hsc_y3/data/boss_dr12_cmass_kaiser_rsd.dataset"
boss_ps = "/home/tanida/cocoa_v41/cocoa/Cocoa/projects/mg_hsc_y3/data/Power_Spectrum_cmass_ngc_v5.txt"

def read_dataset_params(path):
    Nk = None; k_min = None; k_max = None
    with open(path, "r") as f:
        for line in f:
            s = line.split("#",1)[0].strip()
            if "=" in s:
                k, v = [t.strip() for t in s.split("=",1)]
                if k == "Nk":
                    Nk = int(v)
                elif k == "k_min":
                    k_min = float(v)
                elif k == "k_max":
                    k_max = float(v)
    if None in (Nk, k_min, k_max):
        raise RuntimeError("Nk/k_min/k_max が .dataset から取得できませんでした")
    return Nk, k_min, k_max

def read_boss_kcols(path, Nwant=40):
    k_center = []
    k_eff = []
    with open(path, "r") as f:
        for line in f:
            if line.startswith("#") or len(line.strip())==0:
                continue
            parts = line.split()
            if len(parts) >= 2:
                k_center.append(float(parts[0]))
                k_eff.append(float(parts[1]))
                if len(k_center) >= Nwant:
                    break
    return np.asarray(k_center), np.asarray(k_eff)

Nk, kmin, kmax = read_dataset_params(dataset)
k_cocoa = np.linspace(kmin, kmax, Nk)  # いまのCOCOAは線形ビン
k_center, k_eff = read_boss_kcols(boss_ps, Nwant=Nk)

print(f"N = {Nk}")
print(f"COCOA k: [{k_cocoa[0]:.6f}, ..., {k_cocoa[-1]:.6f}]")
print(f"BOSS  k_center: [{k_center[0]:.6f}, ..., {k_center[-1]:.6f}]")
print(f"BOSS  k_eff   : [{k_eff[0]:.6f}, ..., {k_eff[-1]:.6f}]")

diff_center = k_cocoa - k_center
diff_keff   = k_cocoa - k_eff

print("\nmax|COCOA - k_center| =", np.max(np.abs(diff_center)))
print("max|COCOA - k_eff   | =", np.max(np.abs(diff_keff)))

# 必要なら行ごとの一覧を出力
print("\nidx   k_cocoa         k_center        k_eff           dc=cocoa-center     de=cocoa-keff")
for i in range(Nk):
    print(f"{i:2d}  {k_cocoa[i]:.9f}  {k_center[i]:.9f}  {k_eff[i]:.9f}  {diff_center[i]: .3e}  {diff_keff[i]: .3e}")
