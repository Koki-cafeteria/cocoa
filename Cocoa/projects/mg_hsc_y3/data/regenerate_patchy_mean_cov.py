#!/usr/bin/env python3
import glob
import os

import numpy as np


DATA_GLOB = "beutler_data/ps1D_patchy_NGC_z1_COMPnbar_TSC_V6C_*_renorm.dat"
K_BINS_FILE = "boss_dr12_cmass_k_bins_0.02_0.20.txt"
OUT_MEAN = "boss_dr12_NGC_z1_patchy_mean_vector_0p02_0p20.txt"
OUT_COV = "boss_dr12_NGC_z1_patchy_cov_0p02_0p20.txt"


def _is_float(token: str) -> bool:
    try:
        float(token)
        return True
    except ValueError:
        return False


def load_patchy_file(path: str):
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if len(parts) != 13:
                continue
            if not _is_float(parts[0]):
                continue
            rows.append([float(x) for x in parts])
    if not rows:
        raise RuntimeError(f"No numeric rows found in {path}")
    arr = np.array(rows, dtype=float)

    # Header definition (0-based columns):
    # 0:k_central, 1:k_mean, 2:P0, 3:sigP0, 4:P1_im, 5:sigP1,
    # 6:P2, 7:sigP2, 8:P3_im, 9:sigP3, 10:P4, 11:sigP4, 12:N_modes
    k = arr[:, 1]
    p0 = arr[:, 2]
    p2 = arr[:, 6]
    p4 = arr[:, 10]
    return k, p0, p2, p4


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    os.chdir(here)

    files = sorted(glob.glob(DATA_GLOB))
    if not files:
        raise RuntimeError(f"No files matched: {DATA_GLOB}")

    k_target = np.loadtxt(K_BINS_FILE, dtype=float)
    if k_target.ndim != 1:
        k_target = k_target.reshape(-1)

    all_vectors = []
    for path in files:
        k, p0, p2, p4 = load_patchy_file(path)
        order = np.argsort(k)
        k = k[order]
        p0 = p0[order]
        p2 = p2[order]
        p4 = p4[order]

        p0_i = np.interp(k_target, k, p0)
        p2_i = np.interp(k_target, k, p2)
        p4_i = np.interp(k_target, k, p4)

        vec = np.empty(3 * len(k_target), dtype=float)
        vec[0::3] = p0_i
        vec[1::3] = p2_i
        vec[2::3] = p4_i
        all_vectors.append(vec)

    mat = np.asarray(all_vectors, dtype=float)
    mean = np.mean(mat, axis=0)
    cov = np.cov(mat, rowvar=False)

    np.savetxt(OUT_MEAN, mean, fmt="%.18e")
    np.savetxt(OUT_COV, cov, fmt="%.18e")

    print(f"Loaded mocks: {mat.shape[0]}")
    print(f"Data-vector size: {mat.shape[1]}")
    print(f"Wrote mean: {OUT_MEAN}")
    print(f"Wrote cov : {OUT_COV}")
    print("First triplet [P0, P2, P4] at k1:")
    print(mean[0], mean[1], mean[2])


if __name__ == "__main__":
    main()
