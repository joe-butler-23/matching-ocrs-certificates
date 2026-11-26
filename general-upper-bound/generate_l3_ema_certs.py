#!/usr/bin/env python3
"""
Generator for Round-3 Level-3 EMA Cubic-Potential General Matching OCRS Certificates.

Generates three exact-integer certificates inside this bundle directory:
  1. cert_l3_ema_paley_t1_lean.json  (T2=1, Paley 2-coset, 8 free dual-weight numerators, q=100_000_007 == 3 mod 4)
  2. cert_l3_ema_general_t20.json    (T2=20, 21 cosets, 46 free dual-weight numerators, q=10_000_000_000_411 == 1 mod 21)
  3. cert_l3_ema_general_t40.json    (T2=40, 41 cosets, 86 free dual-weight numerators, q=10_000_000_000_313 == 1 mod 41)

Each certificate also specifies structural parameters (T2, q, P, K0, N_L, L, beta_den, eps1_den, eps2_den, D_int).
Zero `assert` statements — all checks use explicit `require()` raising `ValueError`.
"""
from __future__ import annotations

import json
import math
import os
import time
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog, minimize_scalar
from scipy.sparse import csc_matrix


def require(cond: bool, msg: str) -> None:
    if not cond:
        raise ValueError(f"GENERATOR CHECK FAILED: {msg}")


def eval_batch_ema_cubic(S_arr: np.ndarray, wa, wb, wA, wB, aA, aB, wdir, kappa):
    T2 = len(wA)
    dt = 0.5 / T2
    p = -np.expm1(-S_arr * dt)
    lam = kappa * (S_arr * S_arr)
    V0 = np.zeros_like(S_arr)
    V1 = aB * S_arr + 0.5 * lam
    V2 = aA * S_arr + 0.5 * lam
    V3 = wdir + (aA + aB) * S_arr + lam
    for t in range(T2 - 1, -1, -1):
        V3 = V3 + p * np.maximum(0.0, wB[t] + V2 - V3)
        V1 = V1 + p * np.maximum(0.0, wB[t] + V0 - V1)
        V3 = V3 + p * np.maximum(0.0, wA[t] + V1 - V3)
        V2 = V2 + p * np.maximum(0.0, wA[t] + V0 - V2)
    Vmid = V3 + 0.5 * np.maximum(0.0, wb + V2 - V3)
    Vpair = Vmid + 0.5 * np.maximum(0.0, wa + V1 - Vmid)
    return Vpair - lam * S_arr


def eval_scalar_ema_cubic(S: float, wa, wb, wA, wB, aA, aB, wdir, kappa) -> float:
    T2 = len(wA)
    dt = 0.5 / T2
    p = -math.expm1(-S * dt)
    lam = kappa * S * S
    V0 = 0.0
    V1 = aB * S + 0.5 * lam
    V2 = aA * S + 0.5 * lam
    V3 = wdir + (aA + aB) * S + lam
    for t in range(T2 - 1, -1, -1):
        g = wB[t] + V2 - V3
        if g > 0:
            V3 += p * g
        g = wB[t] + V0 - V1
        if g > 0:
            V1 += p * g
        g = wA[t] + V1 - V3
        if g > 0:
            V3 += p * g
        g = wA[t] + V0 - V2
        if g > 0:
            V2 += p * g
    gb = wb + V2 - V3
    Vmid = V3 + (0.5 * gb if gb > 0 else 0.0)
    ga = wa + V1 - Vmid
    Vpair = Vmid + (0.5 * ga if ga > 0 else 0.0)
    return Vpair - lam * S


def solve_adaptive_ema_cubic(T2: int, max_rounds: int = 8):
    dt = 0.5 / T2
    S_pts = sorted(
        set(
            list(np.linspace(0.0, 0.50, 20))
            + list(np.linspace(0.50, 0.72, 75))
            + list(np.linspace(0.72, 1.0, 20))
        )
    )
    S_dense = np.linspace(0.0, 1.0, 20001)
    h = 1.0 / 20000.0
    best_sup = 1.0
    best_sol = None

    for _rnd in range(max_rounds):
        vid = 1
        i_wa, i_wb = vid, vid + 1
        vid += 2
        i_wA = list(range(vid, vid + T2))
        vid += T2
        i_wB = list(range(vid, vid + T2))
        vid += T2
        i_aA, i_aB, i_wdir, i_kap = vid, vid + 1, vid + 2, vid + 3
        vid += 4

        re, ce, ve, be, ru, cu, vu, bu = [], [], [], [], [], [], [], []
        r_eq = r_ub = 0

        def eq(terms, rhs):
            nonlocal r_eq
            for c, v in terms:
                re.append(r_eq); ce.append(c); ve.append(v)
            be.append(rhs); r_eq += 1

        def ub(terms, rhs):
            nonlocal r_ub
            for c, v in terms:
                ru.append(r_ub); cu.append(c); vu.append(v)
            bu.append(rhs); r_ub += 1

        eq(
            [(i_wa, 0.5), (i_wb, 0.5), (i_aA, 1.0), (i_aB, 1.0), (i_wdir, 1.0)]
            + [(i_wA[t], dt) for t in range(T2)]
            + [(i_wB[t], dt) for t in range(T2)],
            1.0,
        )

        for s in S_pts:
            s2 = s * s
            s3 = s2 * s
            p = -math.expm1(-s * dt)
            V1, V2, V3 = vid, vid + 1, vid + 2
            vid += 3
            eq([(V1, 1), (i_aB, -s), (i_kap, -0.5 * s2)], 0)
            eq([(V2, 1), (i_aA, -s), (i_kap, -0.5 * s2)], 0)
            eq([(V3, 1), (i_wdir, -1), (i_aA, -s), (i_aB, -s), (i_kap, -s2)], 0)
            for t in range(T2 - 1, -1, -1):
                mB3, mB1, V3m, V1n = vid, vid + 1, vid + 2, vid + 3
                vid += 4
                ub([(i_wB[t], 1), (V2, 1), (V3, -1), (mB3, -1)], 0)
                eq([(V3m, 1), (V3, -1), (mB3, -p)], 0)
                ub([(i_wB[t], 1), (V1, -1), (mB1, -1)], 0)
                eq([(V1n, 1), (V1, -1), (mB1, -p)], 0)
                mA3, mA2, V3n, V2n = vid, vid + 1, vid + 2, vid + 3
                vid += 4
                ub([(i_wA[t], 1), (V1n, 1), (V3m, -1), (mA3, -1)], 0)
                eq([(V3n, 1), (V3m, -1), (mA3, -p)], 0)
                ub([(i_wA[t], 1), (V2, -1), (mA2, -1)], 0)
                eq([(V2n, 1), (V2, -1), (mA2, -p)], 0)
                V1, V2, V3 = V1n, V2n, V3n
            mb, Vmid, ma, Vpair = vid, vid + 1, vid + 2, vid + 3
            vid += 4
            ub([(i_wb, 1), (V2, 1), (V3, -1), (mb, -1)], 0)
            eq([(Vmid, 1), (V3, -1), (mb, -0.5)], 0)
            ub([(i_wa, 1), (V1, 1), (Vmid, -1), (ma, -1)], 0)
            eq([(Vpair, 1), (Vmid, -1), (ma, -0.5)], 0)
            ub([(Vpair, 1), (i_kap, -s3), (0, -1)], 0)

        c_obj = np.zeros(vid)
        c_obj[0] = 1.0
        bounds = [(0.0, None)] * vid
        bounds[0] = (None, None)
        res = linprog(
            c_obj,
            A_ub=csc_matrix((vu, (ru, cu)), shape=(r_ub, vid)),
            b_ub=np.array(bu),
            A_eq=csc_matrix((ve, (re, ce)), shape=(r_eq, vid)),
            b_eq=np.array(be),
            bounds=bounds,
            method="highs",
        )
        require(res.success, f"Adaptive LP failed: {res.message}")
        x = res.x
        wa = float(x[i_wa])
        wb = float(x[i_wb])
        wA = np.array([float(x[j]) for j in i_wA])
        wB = np.array([float(x[j]) for j in i_wB])
        aA = float(x[i_aA])
        aB = float(x[i_aB])
        wdir = float(x[i_wdir])
        kap = float(x[i_kap])
        lp_val = float(res.fun)

        vals = eval_batch_ema_cubic(S_dense, wa, wb, wA, wB, aA, aB, wdir, kap)
        loc_idx = [int(np.argmax(vals))]
        for k in range(1, len(vals) - 1):
            if vals[k] >= vals[k - 1] and vals[k] >= vals[k + 1] and vals[k] > lp_val + 1e-8:
                loc_idx.append(k)
        peaks = []
        seen = set()
        for idx in sorted(set(loc_idx), key=lambda i: -vals[i])[:12]:
            a_b = max(0.0, S_dense[idx] - 2 * h)
            b_b = min(1.0, S_dense[idx] + 2 * h)
            opt = minimize_scalar(
                lambda s: -eval_scalar_ema_cubic(s, wa, wb, wA, wB, aA, aB, wdir, kap),
                bounds=(a_b, b_b),
                method="bounded",
                options={"xatol": 1e-12},
            )
            sp = float(opt.x)
            vp = float(-opt.fun)
            key = round(sp, 6)
            if key not in seen:
                seen.add(key)
                peaks.append((vp, sp))
        peaks.sort(key=lambda t: -t[0])
        true_sup = peaks[0][0]
        if true_sup < best_sup:
            best_sup = true_sup
            best_sol = (wa, wb, wA.copy(), wB.copy(), aA, aB, wdir, kap)
        if true_sup - lp_val < 1e-7:
            break
        for _, sp in peaks[:8]:
            S_pts.append(sp)
        S_pts = sorted(set(S_pts))
    return best_sup, best_sol


def quantize_cert(
    T2: int,
    q: int,
    P: int,
    K0: int,
    N_L: int,
    L: int,
    beta_den: int,
    eps1_den: int,
    eps2_den: int,
    D_scale: int,
    claimed_ub_str: str,
    out_path: str,
) -> None:
    t0 = time.time()
    sup_cont, (wa, wb, wA, wB, aA, aB, wdir, kap) = solve_adaptive_ema_cubic(T2)
    D_base = D_scale
    wa_int = max(0, int(round(wa * D_base)))
    wb_int = max(0, int(round(wb * D_base)))
    if (wa_int + wb_int) % 2 != 0:
        wb_int += 1
    wA_ints = [max(0, int(round(float(v) * D_base))) for v in wA]
    wB_ints = [max(0, int(round(float(v) * D_base))) for v in wB]
    rem = (sum(wA_ints) + sum(wB_ints)) % (2 * T2)
    if rem != 0:
        wA_ints[0] += 2 * T2 - rem
    aA_int = max(0, int(round(aA * D_base)))
    aB_int = max(0, int(round(aB * D_base)))
    wdir_int = max(0, int(round(wdir * D_base)))
    D_int = (wa_int + wb_int) // 2 + (sum(wA_ints) + sum(wB_ints)) // (2 * T2) + aA_int + aB_int + wdir_int
    kappa_int = max(0, int(round(kap * D_int)))

    norm = (
        Fraction(wa_int + wb_int, 2 * D_int)
        + Fraction(sum(wA_ints) + sum(wB_ints), 2 * T2 * D_int)
        + Fraction(aA_int + aB_int + wdir_int, D_int)
    )
    require(norm == Fraction(1), f"Exact rational normalization = {norm} != 1")

    Sd = np.linspace(0.0, 1.0, 40001)
    vals = eval_batch_ema_cubic(
        Sd,
        wa_int / D_int,
        wb_int / D_int,
        np.array(wA_ints) / D_int,
        np.array(wB_ints) / D_int,
        aA_int / D_int,
        aB_int / D_int,
        wdir_int / D_int,
        kappa_int / D_int,
    )
    q_sup = float(np.max(vals))

    cert = {
        "schema": "ocrs_l3_ema_cubic_v1",
        "T2": T2,
        "q": q,
        "P": P,
        "K0": K0,
        "N_L": N_L,
        "L": L,
        "beta_den": beta_den,
        "eps1_den": eps1_den,
        "eps2_den": eps2_den,
        "D_int": D_int,
        "wa_int": wa_int,
        "wb_int": wb_int,
        "wA_ints": wA_ints,
        "wB_ints": wB_ints,
        "aA_int": aA_int,
        "aB_int": aB_int,
        "wdir_int": wdir_int,
        "kappa_int": kappa_int,
        "claimed_ub": claimed_ub_str,
    }
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(cert, f, indent=2)
    print(f"Wrote {out_path} ({time.time()-t0:.1f}s): cont_sup={sup_cont:.7f}, quant_sup={q_sup:.7f}, "
          f"D_int={D_int}, kappa={kappa_int/D_int:.6f}, claimed_ub={claimed_ub_str}")


def main() -> None:
    bundle_dir = os.path.dirname(os.path.abspath(__file__))
    quantize_cert(
        T2=1, q=100_000_007, P=3, K0=10_000, N_L=10_000, L=10**11,
        beta_den=10**6, eps1_den=10**9, eps2_den=10**9,
        D_scale=100_000, claimed_ub_str="0.37975",
        out_path=os.path.join(bundle_dir, "cert_l3_ema_paley_t1_lean.json"),
    )
    quantize_cert(
        T2=20, q=10_000_000_000_411, P=3, K0=10_000, N_L=10_000, L=10**12,
        beta_den=10**7, eps1_den=10**10, eps2_den=10**10,
        D_scale=1_000_000, claimed_ub_str="0.378685",
        out_path=os.path.join(bundle_dir, "cert_l3_ema_general_t20.json"),
    )
    quantize_cert(
        T2=40, q=10_000_000_000_313, P=3, K0=10_000, N_L=10_000, L=10**12,
        beta_den=10**7, eps1_den=10**10, eps2_den=10**10,
        D_scale=1_000_000, claimed_ub_str="0.378685",
        out_path=os.path.join(bundle_dir, "cert_l3_ema_general_t40.json"),
    )


if __name__ == "__main__":
    main()
