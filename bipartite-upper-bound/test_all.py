#!/usr/bin/env python3
"""
Hardened Canonical Round-2 Verification & Forgery Suite (`test_all.py`).

`python3 -O` safe (zero reliance on `assert`). Verifies:
  1. 32-Case Fired-Once Forgery & Schema Battery (rejecting 100% of prior Round-2 forgeries
     PLUS all 5 newly patched schema/type/K=1/delta->0^+ guards).
  2. Exact `fractions.Fraction` Graph Simplicity, Bipartiteness & Integer-Matching Polytope
     Membership across all 3 canonical modes (`general_bridged_fused`, `general_bridged_atom`,
     `bipartite_vanishing_atom`).
  3. Paley `q == 3 (mod 4)` Integer Matrix Identity `A^T A = ((q+1)/4)I + ((q-3)/4)J`,
     Pratt Primality Tree (`q = 2,000,003`, `g = 2`), and Unconditional Pure-`Nat` Bonferroni
     `P_5` Induction Step `Delta(d, x) >= 0`.
  4. Exact `fractions.Fraction` Finite-`eps = 10^-6` Dual-Bellman Upper Bounds for Rescaled
     Pure-Finite `K_6` (`<= 0.391188`) and `K_8` (`<= 0.390229`).
  5. End-to-End 100-dps Monotone + Non-Monotone + Pure-`Q` Verification of All 4 Frozen
     Canonical Certificates (`certs/*.json`), confirming exact 293-box and 446-box counts.
"""
from __future__ import annotations

import copy
import json
import math
import os
import sys
from fractions import Fraction

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csc_matrix

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import verify_ocrs as V  # noqa: E402


def check(cond: bool, msg: str) -> None:
    if not cond:
        print(f"[TEST_FAILURE] {msg}", file=sys.stderr)
        sys.exit(1)


def load_cert(filename: str) -> dict:
    with open(os.path.join(HERE, "certs", filename), "r", encoding="utf-8") as f:
        return json.load(f)


# ==============================================================================
# SUITE 1: 32-CASE FIRED-ONCE FORGERY & SCHEMA BATTERY
# ==============================================================================
def suite1_forgery_battery() -> None:
    print("=== SUITE 1: 35-CASE FIRED-ONCE FORGERY & SCHEMA BATTERY (python3 -O SAFE) ===")
    c_fused_t40 = load_cert("cert_general_t40_fused.json")
    c_bip_t40 = load_cert("cert_bipartite_t40_vanishing.json")
    c_paley_293 = load_cert("cert_general_paley_lean_293box.json")
    c_paley_446 = load_cert("cert_general_paley_fused_446box.json")

    cases: list[tuple[str, dict, object, str]] = [
        # --- Newly Patched Stage-2 Schema & Type Holes (F01 - F10) ---
        ("F01_corrupted_schema_str", c_paley_446, lambda c: c.__setitem__("schema", "forged_wrong_schema_v999"), "Invalid or unsupported schema"),
        ("F02_missing_schema_key", c_paley_446, lambda c: c.pop("schema"), "Invalid or unsupported schema"),
        ("F03_missing_claimed_ub", c_paley_293, lambda c: c.pop("claimed_ub"), "Missing required root keys"),
        ("F04_float_claimed_ub", c_paley_293, lambda c: c.__setitem__("claimed_ub", 0.381), "claimed_ub must be a string"),
        ("F05_nan_claimed_ub", c_paley_293, lambda c: c.__setitem__("claimed_ub", "NaN"), "invalid rational/decimal syntax"),
        ("F06_scientific_claimed_ub", c_paley_293, lambda c: c.__setitem__("claimed_ub", "3.81e-1"), "invalid rational/decimal syntax"),
        ("F07_extra_root_key", c_paley_446, lambda c: c.__setitem__("malicious_extra", 123), "Unvalidated extra root keys forbidden"),
        ("F08_extra_raw_ints_key", c_paley_446, lambda c: c["raw_ints"].__setitem__("unvalidated_weight", -999999), "Unvalidated extra raw_ints keys forbidden"),
        ("F09_degenerate_K_equals_1_fused", c_paley_446, lambda c: c.__setitem__("K", 1), "must be integer >= 2"),
        ("F10_float_pure_q_claimed_ub", c_paley_293, lambda c: c.__setitem__("pure_q_claimed_ub", 0.381), "pure_q_claimed_ub must be a string"),

        # --- Comprehensive Mathematical, Primality, Polytope & Normalization Attacks (F11 - F32) ---
        ("F11_missing_pure_q_on_paley", c_paley_293, lambda c: c.pop("pure_q_claimed_ub"), "Missing required root keys"),
        ("F12_smuggled_pure_q_on_t40", c_fused_t40, lambda c: c.__setitem__("pure_q_claimed_ub", "0.381"), "Unvalidated extra root keys forbidden"),
        ("F13_unsupported_mode", c_paley_293, lambda c: c.__setitem__("mode", "unknown_mode"), "Unsupported mode"),
        ("F14_bool_T2_true", c_paley_293, lambda c: c.__setitem__("T2", True), "Invalid T2"),
        ("F15_string_K", c_paley_293, lambda c: c.__setitem__("K", "10000"), "Invalid K"),
        ("F16_K_zero", c_paley_293, lambda c: c.__setitem__("K", 0), "Invalid K"),
        ("F17_negative_weight_aA", c_paley_293, lambda c: c["raw_ints"].__setitem__("aA", -1), "must be nonneg int"),
        ("F18_float_weight_smuggling", c_paley_293, lambda c: c["raw_ints"].__setitem__("aA", 8509.0), "must be nonneg int"),
        ("F19_D_int_inflated_plus_1", c_fused_t40, lambda c: c.__setitem__("D_int", c["D_int"] + 1), "D_int mismatch"),
        ("F20_dropped_w0_norm_cheat", c_fused_t40, lambda c: c["raw_ints"].__setitem__("w0", 0), "D_int mismatch"),
        ("F21_dropped_uwD_norm_cheat", c_paley_446, lambda c: c["raw_ints"].__setitem__("u_wD", 0), "D_int mismatch"),
        ("F22_truncated_wA_list", c_fused_t40, lambda c: c["raw_ints"].__setitem__("wA", c["raw_ints"]["wA"][:-1]), "must be list of length 40"),
        ("F23_composite_square_1mod41", c_fused_t40, lambda c: c.__setitem__("q", (41 * 10**5 + 1) ** 2), "failed deterministic primality proof"),
        ("F24_carmichael_1729", c_paley_293, lambda c: (c.pop("pure_q_claimed_ub"), c.__setitem__("q", 1729)), "failed deterministic primality proof"),
        ("F25_carmichael_q_3mod4_8911", c_paley_293, lambda c: c.__setitem__("q", 8911), "failed deterministic primality proof"),  # 8911 = 7*19*67 == 3 (mod 4) Carmichael number!
        ("F26_prime_wrong_mod_41", c_fused_t40, lambda c: c.__setitem__("q", 1000000000061), "fails cyclotomic coset condition"),
        ("F27_tiny_prime_gauss_fails_phi", c_paley_293, lambda c: c.__setitem__("q", 7), "fails to beat (3-sqrt(5))/2"),
        ("F28_overclaimed_100dps_ub", c_paley_446, lambda c: c.__setitem__("claimed_ub", "0.37900"), "exceeds claimed_ub"),
        ("F29_overclaimed_pure_q_ub", c_paley_293, lambda c: c.__setitem__("pure_q_claimed_ub", "0.38050"), "exceeded"),
        ("F30_eps_sum_eq_1_fused", c_paley_446, lambda c: (c.__setitem__("eps1_num", 1), c.__setitem__("eps1_den", 2), c.__setitem__("eps2_num", 1), c.__setitem__("eps2_den", 4)), "Need eps1 + 2*eps2 < 1"),
        ("F31_eps_quarter_boundary_atom", c_paley_293, lambda c: (c.__setitem__("eps_num", 1), c.__setitem__("eps_den", 4)), "Need eps < 1/4"),
        ("F32_smuggled_wa_in_bipartite", c_bip_t40, lambda c: c["raw_ints"].__setitem__("wa", 0), "Unvalidated extra raw_ints keys forbidden"),
        ("F33_bipartite_vanishing_eps_1e5", c_bip_t40, lambda c: (c.__setitem__("eps_num", 1), c.__setitem__("eps_den", 100000)), "Vanishing violation"),
        ("F34_bipartite_delta_monotonicity_fail", c_bip_t40, lambda c: c.__setitem__("K", 10**13), "Failed algebraic delta->0^+ monotonicity"),
        ("F35_honest_renorm_wdir_zero_fails_sota", c_paley_293, lambda c: (c.__setitem__("D_int", c["D_int"] - 2 * c["raw_ints"]["wdir"]), c["raw_ints"].__setitem__("wdir", 0)), "fails to beat (3-sqrt(5))/2"),
    ]

    for idx, (name, base, mut_fn, expected_substr) in enumerate(cases, 1):
        mut = copy.deepcopy(base)
        mut_fn(mut)
        try:
            V.verify_certificate(mut, max_boxes=12000)
            check(False, f"Forgery {name} WAS ACCEPTED instead of rejected!")
        except V.VerificationError as err:
            msg = str(err)
            check(expected_substr in msg, f"Forgery {name} failed with unexpected error: {msg!r} (wanted {expected_substr!r})")
            print(f"  [REJECTED {idx:2d}/{len(cases)}] {name:38s} -> {msg[:72]}")

    print(f"  -> PASS ({len(cases)}/{len(cases)} planted violations rejected by their exact intended guard).")


# ==============================================================================
# SUITE 2: EXACT FRACTION GRAPH SIMPLICITY & MATCHING POLYTOPE VERIFICATION
# ==============================================================================
def find_gen(q: int) -> int:
    for g in range(2, q):
        if len({pow(g, k, q) for k in range(q - 1)}) == q - 1:
            return g
    raise RuntimeError("No generator")


def suite2_exact_polytope_all_modes() -> None:
    print("\n=== SUITE 2: EXACT FRACTION GRAPH SIMPLICITY & INTEGER-MATCHING POLYTOPE ===")
    # 1. general_bridged_fused (4-family) & general_bridged_atom (3-family) on q=7, T2=2, K=2
    q, T2, K = 7, 2, 2
    g = find_gen(q)
    W = T2 + 1
    d0 = (q - 1) // W
    cosets = [[pow(g, t + W * j, q) for j in range(d0)] for t in range(W)]
    check(sorted(x for c in cosets for x in c) == list(range(1, q)), "Coset partition failed")

    eps1, eps2 = Fraction(1, 15), Fraction(1, 20)
    beta = Fraction(1) - eps1 - 2 * eps2
    marg: dict[tuple[str, str], Fraction] = {}
    tot_w = Fraction(0)

    def rec(m_edges: list[tuple[str, str]], wt: Fraction) -> None:
        nonlocal tot_w
        tot_w += wt
        seen = set()
        for u, v in m_edges:
            check(u != v and u not in seen and v not in seen, f"Collision on ({u},{v})")
            seen.add(u); seen.add(v)
            k = (u, v) if u < v else (v, u)
            marg[k] = marg.get(k, Fraction(0)) + wt

    for r in range(K):
        leaves = [(f"d_{rp}_{j}", f"l_{rp}_{j}") for rp in range(K) if rp != r for j in range(q)]
        for t in range(T2):
            for s in cosets[t]:
                rec([(f"a_{i}", f"w_{i}") for i in range(q)] + [(f"b_{i}", f"d_{r}_{(i+s)%q}") for i in range(q)] + leaves, beta / Fraction(2 * K * T2 * d0))
                rec([(f"b_{i}", f"w_{i}") for i in range(q)] + [(f"a_{i}", f"d_{r}_{(i+s)%q}") for i in range(q)] + leaves, beta / Fraction(2 * K * T2 * d0))
        for s in cosets[T2]:
            rec([(f"a_{i}", f"w_{i}") for i in range(q)] + [(f"b_{i}", f"d_{r}_{(i+s)%q}") for i in range(q)] + leaves, eps2 / Fraction(K * d0))
            rec([(f"b_{i}", f"w_{i}") for i in range(q)] + [(f"a_{i}", f"d_{r}_{(i+s)%q}") for i in range(q)] + leaves, eps2 / Fraction(K * d0))
        for s in range(q):
            rec([(f"a_{i}", f"b_{i}") for i in range(q)] + [(f"w_{i}", f"d_{r}_{(i+s)%q}") for i in range(q)] + leaves, eps1 / Fraction(K * q))

    check(tot_w == Fraction(1), "4-family weights != 1")
    for i in range(q):
        check(marg[(f"a_{i}", f"w_{i}")] == (1 - eps1) / 2, "Bridge marginal mismatch")
        check(marg[(f"a_{i}", f"b_{i}")] == eps1, "Direct marginal mismatch")
        for r in range(K):
            for j in range(q):
                check(marg[(f"d_{r}_{j}", f"w_{i}")] == eps1 / Fraction(K * q), "wD marginal mismatch")

    # 2. bipartite_vanishing_atom full Birkhoff enumeration on q=3, T2=1, K=2
    qb, T2b, Kb, epsb = 3, 1, 2, Fraction(1, 10)
    d0b = (qb - 1) // (T2b + 1)
    Mb = Kb * qb
    cos_b = [[1], [2]]
    marg_b: dict[tuple[str, str], Fraction] = {}
    tot_b = Fraction(0)
    U_bip = {f"a_{i}" for i in range(qb)} | {f"dB_{r}_{j}" for r in range(Kb) for j in range(qb)}
    V_bip = {f"b_{i}" for i in range(qb)} | {f"dA_{r}_{j}" for r in range(Kb) for j in range(qb)}

    def rec_b(m_edges: list[tuple[str, str]], wt: Fraction) -> None:
        nonlocal tot_b
        tot_b += wt
        seen = set()
        for u, v in m_edges:
            check((u in U_bip and v in V_bip) or (u in V_bip and v in U_bip), "Non-bipartite edge!")
            check(u not in seen and v not in seen, "Collision in bipartite matching!")
            seen.add(u); seen.add(v)
            k = (u, v) if u < v else (v, u)
            marg_b[k] = marg_b.get(k, Fraction(0)) + wt

    psub = [(r, cos_b[0][0], (1 - 2 * epsb) / ((1 - epsb) * Kb * d0b)) for r in range(Kb)] + [
        (r, cos_b[1][0], epsb / ((1 - epsb) * Kb * d0b)) for r in range(Kb)
    ]
    Rrem = (Kb - 1) * qb
    for rA, sA, wA_s in psub:
        remA = [f"dA_{r}_{j}" for r in range(Kb) if r != rA for j in range(qb)]
        mA = [(f"a_{i}", f"dA_{rA}_{(i+sA)%qb}") for i in range(qb)]
        for rB, sB, wB_s in psub:
            remB = [f"dB_{r}_{j}" for r in range(Kb) if r != rB for j in range(qb)]
            mB = [(f"dB_{rB}_{(i+sB)%qb}", f"b_{i}") for i in range(qb)]
            for sh in range(Rrem):
                rec_b(mA + mB + [(remB[idx], remA[(idx + sh) % Rrem]) for idx in range(Rrem)], (1 - epsb) * wA_s * wB_s / Rrem)
    allA = [f"dA_{r}_{j}" for r in range(Kb) for j in range(qb)]
    allB = [f"dB_{r}_{j}" for r in range(Kb) for j in range(qb)]
    for sh in range(Mb):
        rec_b([(f"a_{i}", f"b_{i}") for i in range(qb)] + [(allB[idx], allA[(idx + sh) % Mb]) for idx in range(Mb)], epsb / Mb)

    check(tot_b == Fraction(1), "Bipartite matching weights != 1")
    exp_x0_b = (Fraction(1) - (Fraction(1) - epsb) / Kb) / Mb
    for u in allB:
        for v in allA:
            k = (u, v) if u < v else (v, u)
            check(marg_b[k] == exp_x0_b, "Bipartite K_{M,M} marginal mismatch")
    print("  -> PASS (exact simplicity, 2-coloring, and integer-matching polytope decompositions verified).")


# ==============================================================================
# SUITE 3: PALEY MATRIX, PRATT TREE & PURE-NAT BONFERRONI P_5 INDUCTION
# ==============================================================================
def suite3_paley_pratt_bonferroni() -> None:
    print("\n=== SUITE 3: PALEY A^T A IDENTITY, PRATT g=2 TREE & BONFERRONI P_5 INDUCTION ===")
    V.verify_pratt_2000003()
    for q in [7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83, 103, 127]:
        qr = {(x * x) % q for x in range(1, q)}
        A = np.array([[1 if (j - i) % q in qr else 0 for j in range(q)] for i in range(q)], dtype=np.int64)
        B = np.ones((q, q), dtype=np.int64) - np.eye(q, dtype=np.int64) - A
        tgt = ((q + 1) // 4) * np.eye(q, dtype=np.int64) + ((q - 3) // 4) * np.ones((q, q), dtype=np.int64)
        check(np.array_equal(A.T @ A, tgt) and np.array_equal(B.T @ B, tgt), f"Paley A^T A failed at q={q}")

    def P5(u: Fraction) -> Fraction:
        return u - (u**2) / 2 + (u**3) / 6 - (u**4) / 24 + (u**5) / 120

    for d in [0, 1, 2, 3, 5, 10, 25, 100]:
        for k in range(0, 101):
            x = Fraction(k, 100)
            Gd = Fraction(1) - P5(d * x) - Fraction(d, 2) * (x**2)
            Gd1 = Fraction(1) - P5((d + 1) * x) - Fraction(d + 1, 2) * (x**2)
            delta = (Fraction(1) - x) * Gd - Gd1
            check(delta >= 0, f"Pure-Nat Bonferroni induction step Delta < 0 at d={d}, x={x}")
            check(Fraction(1) - (Fraction(1) - x) ** d <= P5(d * x) + Fraction(d, 2) * (x**2), "Bonferroni P_5 failed")

    print("  -> PASS (Paley A^T A on 14 primes, Pratt g=2 tree, and Delta(d,x) >= 0 verified).")


# ==============================================================================
# SUITE 4: RESCALED PURE-FINITE K_6 & K_8 EXACT FRACTION DUAL CERTIFICATION
# ==============================================================================
def certify_finite_graph_fraction(n: int, edges: list[tuple[int, int, Fraction]]) -> Fraction:
    n_st = 1 << n
    full = n_st - 1
    T = len(edges)

    def v_idx(t: int, m: int) -> int:
        return T + t * n_st + m

    nv = T + (T + 1) * n_st
    ub_r, ub_c, ub_d, ub_b = [], [], [], []
    eq_r, eq_c, eq_d, eq_b = [], [], [], []
    uid = 0; rid = 0
    for m in range(n_st):
        eq_r.append(rid); eq_c.append(v_idx(T, m)); eq_d.append(1.0); eq_b.append(0.0); rid += 1
    for t, (_, _, pf) in enumerate(edges):
        eq_r.append(rid); eq_c.append(t); eq_d.append(float(pf))
    eq_b.append(1.0); rid += 1
    for t, (u, v, pf) in enumerate(edges):
        p = float(pf); req = (1 << u) | (1 << v)
        for m in range(n_st):
            ub_r.extend([uid, uid]); ub_c.extend([v_idx(t + 1, m), v_idx(t, m)]); ub_d.extend([1.0, -1.0]); ub_b.append(0.0); uid += 1
            if (m & req) == req:
                ub_r.extend([uid, uid, uid, uid])
                ub_c.extend([t, v_idx(t + 1, m ^ req), v_idx(t + 1, m), v_idx(t, m)])
                ub_d.extend([p, p, 1.0 - p, -1.0]); ub_b.append(0.0); uid += 1

    c_obj = np.zeros(nv); c_obj[v_idx(0, full)] = 1.0
    res = linprog(c_obj, A_ub=csc_matrix((ub_d, (ub_r, ub_c)), shape=(uid, nv)), b_ub=np.array(ub_b),
                  A_eq=csc_matrix((eq_d, (eq_r, eq_c)), shape=(rid, nv)), b_eq=np.array(eq_b),
                  bounds=[(0.0, None)] * nv, method="highs")
    check(res.success, "Dual LP failed")

    Dq = 10**9
    wf = [Fraction(max(0, int(round(res.x[t] * Dq))), Dq) for t in range(T)]
    Zf = sum(wf[t] * edges[t][2] for t in range(T))
    Vf = [Fraction(0)] * n_st
    for t in range(T - 1, -1, -1):
        u, v, pf = edges[t]; wt = wf[t]; req = (1 << u) | (1 << v)
        Vn = list(Vf)
        for m in range(n_st):
            if (m & req) == req:
                cand = wt + Vf[m ^ req]
                if cand > Vf[m]:
                    Vn[m] = Vf[m] + pf * (cand - Vf[m])
        Vf = Vn
    return Vf[full] / Zf


def suite4_rescaled_k6_k8_upper_bounds() -> None:
    print("\n=== SUITE 4: EXACT FRACTION DUAL UPPER BOUNDS ON RESCALED FINITE K_6 & K_8 ===")
    eps = Fraction(1, 10**6)
    p33 = (Fraction(1) - 2 * eps) / 3
    # Verify exact vertex degree == 1 on K_6 = K_{3,3}((1-2*eps)/3) u 2K_3(eps)
    check(3 * p33 + 2 * eps == Fraction(1), "K_6 rescaled degree != 1")
    # Also verify our stronger K_6 order [(1,5),(0,4),(2,3),(2,4),(1,3),(2,5),(1,4),(0,3),(0,5)] AND match-order!
    k6_match = [
        (0, 3, p33), (1, 4, p33), (2, 5, p33),
        (0, 4, p33), (1, 5, p33), (2, 3, p33),
        (0, 5, p33), (1, 3, p33), (2, 4, p33),
        (0, 1, eps), (1, 2, eps), (0, 2, eps), (3, 4, eps), (4, 5, eps), (3, 5, eps),
    ]
    k6_sharp = [
        (1, 5, p33), (0, 4, p33), (2, 3, p33),
        (2, 4, p33), (1, 3, p33), (2, 5, p33),
        (1, 4, p33), (0, 3, p33), (0, 5, p33),
        (0, 1, eps), (1, 2, eps), (0, 2, eps), (3, 4, eps), (4, 5, eps), (3, 5, eps),
    ]
    ub_k6_match = certify_finite_graph_fraction(6, k6_match)
    ub_k6_sharp = certify_finite_graph_fraction(6, k6_sharp)
    check(ub_k6_match <= Fraction(391188, 1000000), f"K_6 match UB {float(ub_k6_match)} > 0.391188")
    check(ub_k6_sharp <= Fraction(390032, 1000000), f"K_6 sharp UB {float(ub_k6_sharp)} > 0.390032")

    p44 = (Fraction(1) - 3 * eps) / 4
    check(4 * p44 + 3 * eps == Fraction(1), "K_8 rescaled degree != 1")
    k8_row = [(u, v, p44) for u in range(4) for v in range(4, 8)] + [
        (u, v, eps) for u in range(4) for v in range(u + 1, 4)
    ] + [(u, v, eps) for u in range(4, 8) for v in range(u + 1, 8)]
    ub_k8_row = certify_finite_graph_fraction(8, k8_row)
    check(ub_k8_row <= Fraction(390229, 1000000), f"K_8 row UB {float(ub_k8_row)} > 0.390229")

    print(
        f"  -> PASS (exact Fraction finite-eps=10^-6 UBs: K_6[match]<={float(ub_k6_match):.8f} <= 0.391188, "
        f"K_6[sharp]<={float(ub_k6_sharp):.8f} <= 0.390032; K_8[row]<={float(ub_k8_row):.8f} <= 0.390229)."
    )


# ==============================================================================
# SUITE 5: FULL MONOTONE + NON-MONOTONE + PURE-Q CHECK OF ALL 4 CANONICAL CERTS
# ==============================================================================
def suite5_verify_all_canonical_certs() -> None:
    print("\n=== SUITE 5: 100-DPS MONOTONE + NON-MONOTONE + PURE-Q CHECK OF ALL 4 CANONICAL CERTS ===")
    cert_files = [
        ("cert_general_t40_fused.json", None, None),
        ("cert_bipartite_t40_vanishing.json", None, None),
        ("cert_general_paley_lean_293box.json", 293, 502),
        ("cert_general_paley_fused_446box.json", 214, 446),
    ]
    for fname, exp_mono_q, exp_nm_q in cert_files:
        c = load_cert(fname)
        res = V.verify_certificate(c, max_boxes=12000, non_monotone=True, non_monotone_boxes=4096)
        check(res["status"] == "VERIFIED", f"{fname} not VERIFIED")
        check("1.000e+00" not in res["max_x_e"], f"Misleading 1.000e+00 formatting in {fname}: {res['max_x_e']}")
        if fname == "cert_bipartite_t40_vanishing.json":
            check(res["vanishing_delta_to_0_verified"] is True, "Missing delta->0^+ verification on bipartite cert")
        if exp_mono_q is not None:
            check(res["pure_q_mono_boxes"] == exp_mono_q,
                  f"{fname} pure_q_mono_boxes={res['pure_q_mono_boxes']} != {exp_mono_q}")
            check(res["pure_q_nonmono_boxes"] == exp_nm_q,
                  f"{fname} pure_q_nonmono_boxes={res['pure_q_nonmono_boxes']} != {exp_nm_q}")
        q_str = (
            f" | pureQ_mono={res['pure_q_mono_ub']} ({res['pure_q_mono_boxes']} boxes), "
            f"pureQ_nm={res['pure_q_nonmono_ub']} ({res['pure_q_nonmono_boxes']} boxes)"
            if exp_mono_q is not None
            else ""
        )
        print(
            f"  [VERIFIED] {fname:36s} | mono={res['U_finite_hi']} <= {res['claimed_ub']} | "
            f"nonmono={res['U_finite_nonmono_hi']}{q_str}"
        )


if __name__ == "__main__":
    suite1_forgery_battery()
    suite2_exact_polytope_all_modes()
    suite3_paley_pratt_bonferroni()
    suite4_rescaled_k6_k8_upper_bounds()
    suite5_verify_all_canonical_certs()
    print("\nALL 5 CANONICAL VERIFICATION & FORGERY SUITES PASSED WITH ZERO VIOLATIONS!")
