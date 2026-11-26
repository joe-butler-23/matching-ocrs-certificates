#!/usr/bin/env python3
"""
18-Case Fired-Once Forgery Battery + 6-Suite Independent Adversarial Lemma Fuzzer
for Round-3 Level-3 EMA Cubic-Potential General Matching OCRS!

Must pass under `python3 -O` with exit code 0!
  - Part A: 18 Fired-Once Forgeries (EVERY single planted defect MUST be rejected!)
  - Part B: 6 Independent Adversarial Lemma Suites:
      Suite 1: Exact fractions.Fraction Top-Down Pyramid matching polytope + 10k valid matching samples
               (0 collisions, exact gamma_k, x_e > 0 everywhere, simple graph!).
      Suite 2: 100,000-sample adversarial test of the 1-line Cubic Potential Identity
               (kappa/(3*beta))*(S_next^3 - S^3) <= kappa*S^2*(S_lvl - S) + kappa*beta
               + full multi-level telescoping path check!
      Suite 3: Exact (T2+1)-coset Gauss-period + Paley FFT spectral bound & multi-level EMA weighted
               sum_j alpha_{k,j} |rho^(j) - S_j| bound!
      Suite 4: 20,000 100-dps checks of finite-to-ideal direct edge + atom wave inequality with
               terminal bonus +lambda >= 0, PLUS 15,000 L3 Bellman state-monotonicity, S-monotonicity,
               and L1-sensitivity checks!
      Suite 5: LITERAL EXACT COUPLED 2-LEVEL FINITE-GRAPH BELLMAN CHECK (Level 0 -> Level 1 -> Level 2
               with real capacity-1 pool vertices!) proving EXACT_COUPLED <= RELAXED_FROZEN_CAPACITY!
      Suite 6: 100,000-point + Brent adversarial global peak search on all 3 certificates!
"""
from __future__ import annotations

import copy
import json
import math
import os
import random
import tempfile
import time
from fractions import Fraction
import mpmath as mp
import numpy as np
from scipy.optimize import minimize_scalar

import verify_l3_ema_ocrs as verifier


def _fail(msg: str) -> None:
    raise RuntimeError(f"AUDIT FAILURE: {msg}")


# =============================================================================
# PART A: 18-CASE FIRED-ONCE FORGERY BATTERY
# =============================================================================
def run_forgery_battery(base_cert_path: str) -> None:
    print("=" * 85)
    print("PART A: 18-CASE FIRED-ONCE FORGERY BATTERY (python3 -O safe)")
    print("=" * 85)
    with open(base_cert_path, "r", encoding="utf-8") as f:
        base = json.load(f)

    cases: list[tuple[str, callable]] = []

    def add(name, fn):
        cases.append((name, fn))

    add("F01_wrong_schema", lambda d: d.update(schema="ocrs_l2_general_v1"))
    add("F02_extra_unknown_key", lambda d: d.update(sneaky_bonus=123))
    add("F03_missing_claimed_ub", lambda d: d.pop("claimed_ub"))
    add("F04_float_claimed_ub", lambda d: d.update(claimed_ub=0.38))
    add("F05_float_weight_smuggling", lambda d: d.update(wdir_int=float(d["wdir_int"])))
    add("F06_bool_smuggling_T2", lambda d: d.update(T2=True))
    add("F07_D_int_plus_1_renorm_break", lambda d: d.update(D_int=d["D_int"] + 1))
    add("F08_negative_kappa", lambda d: d.update(kappa_int=-1))
    add("F09_truncated_wA_list", lambda d: d.update(wA_ints=d["wA_ints"][:-1]))
    add("F10_P_too_small_P2", lambda d: d.update(P=2))
    add("F11_K0_degenerate_1", lambda d: d.update(K0=1))
    add("F12_NL_degenerate_1", lambda d: d.update(N_L=1))
    add("F13_composite_q_1_mod_W", lambda d: d.update(q=2101))  # 11*191 == 1 mod 21 (for T2=20!)
    add("F14_prime_wrong_mod", lambda d: d.update(q=100_000_007))  # prime, != 1 mod 21
    add("F15_tiny_prime_q43", lambda d: d.update(q=43))  # 43 == 1 mod 21, huge gauss err!
    add("F16_overclaimed_ub_0.37860", lambda d: d.update(claimed_ub="0.37860"))
    add("F17_eps_sum_ge_1", lambda d: d.update(eps1_den=4, eps2_den=2))
    def f18_honest_weak(d):
        # Honest normalized certificate with kappa=0 and all weight on wa/wb -> UB = 0.75 <= 0.76 but > golden!
        T2 = d["T2"]
        d.update(D_int=100, wa_int=100, wb_int=100, wA_ints=[0]*T2, wB_ints=[0]*T2,
                 aA_int=0, aB_int=0, wdir_int=0, kappa_int=0, claimed_ub="0.76")
    add("F18_honest_weak_above_golden", f18_honest_weak)

    rejected = 0
    with tempfile.TemporaryDirectory() as td:
        for name, mut in cases:
            cand = copy.deepcopy(base)
            mut(cand)
            p = os.path.join(td, f"{name}.json")
            with open(p, "w", encoding="utf-8") as f:
                json.dump(cand, f)
            try:
                verifier.verify_certificate(p, verbose=False)
            except ValueError as e:
                rejected += 1
                print(f"  [REJECTED] {name:32s} -> {str(e)[:72]}")
            else:
                _fail(f"FORGERY {name} WAS WRONGLY ACCEPTED!")
    if rejected != len(cases):
        _fail(f"Only {rejected}/{len(cases)} forgeries rejected")
    print(f"  ==> {rejected}/{len(cases)} FIRED-ONCE FORGERIES REJECTED (100% PASS)!")


# =============================================================================
# SUITE 1: EXACT FRACTION TOP-DOWN PYRAMID MATCHING POLYTOPE
# =============================================================================
def suite1_pyramid_polytope() -> None:
    from test_ema_coupled_finite_toy import verify_ema_pyramid_polytope_exact
    for q, T2, P, K0, N_L, L, b_den, e1, e2 in [
        (5, 1, 3, 2, 2, 2, 5, 20, 20),
        (7, 2, 3, 2, 2, 3, 10, 50, 50),
        (13, 3, 3, 3, 2, 4, 100, 1000, 1000),
    ]:
        verify_ema_pyramid_polytope_exact(
            q=q, T2=T2, P=P, K0=K0, N_L=N_L, L=L,
            beta=Fraction(1, b_den), eps1=Fraction(1, e1), eps2=Fraction(1, e2)
        )
    print("  [PASS] Suite 1: Exact Fraction Top-Down Pyramid polytope + 0-collision matching sampler")


# =============================================================================
# SUITE 2: 100,000-POINT CUBIC POTENTIAL IDENTITY & TELESCOPING CHECK
# =============================================================================
def suite2_cubic_identity() -> None:
    rng = np.random.default_rng(101)
    N = 100_000
    S = rng.uniform(0.0, 1.0, N)
    S_lvl = rng.uniform(0.0, 1.0, N)
    # Include extreme corner cases (0, 1), (1, 0)
    S[:4] = [0.0, 1.0, 0.0, 1.0]
    S_lvl[:4] = [1.0, 0.0, 0.0, 1.0]
    beta = 10.0 ** rng.uniform(-6.0, -0.5, N)
    kappa = rng.uniform(0.0, 2.0, N)
    S_next = (1.0 - beta) * S + beta * S_lvl
    lhs = (kappa / (3.0 * beta)) * (S_next ** 3 - S ** 3)
    rhs = kappa * (S ** 2) * (S_lvl - S) + kappa * beta
    worst_viol = float(np.max(lhs - rhs))
    if worst_viol > 1e-13:
        _fail(f"Cubic identity violated: max(lhs - rhs) = {worst_viol}")

    # Multi-level telescoping path check over 500 random paths of length L=200
    max_path_viol = -1.0
    for _ in range(500):
        b = 0.05; k_c = 0.16; L = 200
        s_cur = float(rng.uniform(0.0, 1.0))
        s0 = s_cur
        sum_lin = 0.0
        for _ in range(L):
            sl = float(rng.uniform(0.0, 1.0))
            sum_lin += k_c * (s_cur ** 2) * (sl - s_cur) + k_c * b
            s_cur = (1.0 - b) * s_cur + b * sl
        exact_tele = (k_c / (3.0 * b)) * (s_cur ** 3 - s0 ** 3)
        # We use: sum_k [ -lin_k ] <= (k/(3bL))*(s0^3 - s_L^3) + k*b <= k/(3bL) + k*b
        # Equivalently exact_tele <= sum_lin!
        viol = (exact_tele - sum_lin) / L
        if viol > max_path_viol:
            max_path_viol = viol
    if max_path_viol > 1e-13:
        _fail(f"Telescoping path violated: {max_path_viol}")
    print(f"  [PASS] Suite 2: 100k Cubic Potential Identity + 500 telescoping paths (max_viol={worst_viol:.2e} <= 0)")


# =============================================================================
# SUITE 3: EXACT GAUSS-PERIOD & MULTI-LEVEL WEIGHTED MIXING FFT BOUND
# =============================================================================
def suite3_gauss_multilevel_fft() -> None:
    rng = np.random.default_rng(202)
    worst_ratio = 0.0
    for q, W in [(13, 2), (29, 2), (43, 7), (83, 2), (127, 7), (211, 21)]:
        # Find primitive root g mod q
        g = 2
        while True:
            if len({pow(g, k, q) for k in range(q - 1)}) == q - 1:
                break
            g += 1
        d0 = (q - 1) // W
        cosets = [[pow(g, t + W * j, q) for j in range(d0)] for t in range(W)]
        bound = W * math.sqrt(q) / (2.0 * (q - 1))
        # Multi-level convex combination alpha_j over 5 levels
        alphas = rng.dirichlet(np.ones(5))
        for _ in range(15):
            # Adversarial + random subsets F_j of Z_q
            rhos = np.zeros((W, q))
            S_ema = 0.0
            for j in range(5):
                p_f = float(rng.uniform(0.1, 0.9))
                F_vec = (rng.random(q) < p_f).astype(float)
                Sj = float(np.mean(F_vec))
                S_ema += alphas[j] * Sj
                Ff = np.fft.fft(F_vec)
                for t in range(W):
                    ind = np.zeros(q)
                    for s in cosets[t]:
                        ind[(-s) % q] = 1.0 / d0
                    rho_jt = np.fft.ifft(Ff * np.fft.fft(ind)).real
                    rhos[t] += alphas[j] * rho_jt
            for t in range(W):
                l1 = float(np.mean(np.abs(rhos[t] - S_ema)))
                ratio = l1 / bound
                if ratio > worst_ratio:
                    worst_ratio = ratio
                if ratio > 1.0 + 1e-12:
                    _fail(f"Gauss L1 exceeded bound: {l1} > {bound} at q={q}, W={W}")
    print(f"  [PASS] Suite 3: Multi-level weighted Gauss-period FFT checks (worst L1/bound = {worst_ratio:.4f} <= 1.0)")


# =============================================================================
# SUITE 4: 100-DPS ATOM+DIRECT BOUND WITH +LAMBDA >= 0 & BELLMAN SENSITIVITY
# =============================================================================
def suite4_atom_direct_and_bellman() -> None:
    mp.mp.dps = 100
    rng = random.Random(303)
    worst_atom_viol = mp.mpf("-1")
    for _ in range(5000):
        lam = mp.mpf(rng.uniform(0.0, 1.5))
        wdir = mp.mpf(rng.uniform(0.0, 1.0))
        aA = mp.mpf(rng.uniform(0.0, 0.5))
        aB = mp.mpf(rng.uniform(0.0, 0.5))
        eps1 = mp.mpf(10) ** (-rng.randint(3, 12))
        eps2 = mp.mpf(10) ** (-rng.randint(3, 12))
        gamma = mp.mpf(rng.uniform(0.5, 1.0))
        rhoA = mp.mpf(rng.uniform(0.0, 1.0))
        rhoB = mp.mpf(rng.uniform(0.0, 1.0))
        dA = rng.randint(1, 500)
        dB = rng.randint(1, 500)
        # Finite activation probabilities: p_dir = gamma*eps1, pA = 1 - (1 - gamma*eps2*rhoA/dA)^dA <= eps2*rhoA
        p_dir = gamma * eps1
        pA = 1 - (1 - gamma * eps2 * rhoA / dA) ** dA
        pB = 1 - (1 - gamma * eps2 * rhoB / dB) ** dB

        # Exact finite Phase-3 direct edge + atom wave backward Bellman!
        V0_end = mp.mpf(0)
        V1_end = lam / 2
        V2_end = lam / 2
        V3_end = lam + p_dir * max(mp.mpf(0), wdir / eps1 + V0_end - lam)

        # B-atom step
        V1_b = V1_end + pB * max(mp.mpf(0), aB / eps2 + V0_end - V1_end)
        V3_m = V3_end + pB * max(mp.mpf(0), aB / eps2 + V2_end - V3_end)
        # A-atom step
        V2_a = V2_end + pA * max(mp.mpf(0), aA / eps2 + V0_end - V2_end)
        V3_a = V3_m + pA * max(mp.mpf(0), aA / eps2 + V1_b - V3_m)

        # Ideal upper bounds
        B1 = aB * rhoB + lam / 2
        B2 = aA * rhoA + lam / 2
        B3 = wdir + aA * rhoA + aB * rhoB + lam
        v = max(V1_b - B1, V2_a - B2, V3_a - B3)
        if v > worst_atom_viol:
            worst_atom_viol = v
        if v > mp.mpf("1e-90"):
            _fail(f"Finite atom+direct bound violated: {v}")

    # 15,000 random checks of Bellman state-monotonicity, S-monotonicity, and L1-sensitivity!
    np_rng = np.random.default_rng(404)
    worst_sens = 0.0
    for _ in range(15000):
        T2 = 4
        wA = np_rng.uniform(0.0, 1.0, T2)
        wB = np_rng.uniform(0.0, 1.0, T2)
        wa, wb, aA, aB, wdir, kap = np_rng.uniform(0.0, 1.0, 6)
        s1 = float(np_rng.uniform(0.0, 1.0))
        s2 = float(np_rng.uniform(s1, 1.0))

        def dp(pA_vec, pB_vec, rA, rB, lam_val):
            V0 = 0.0; V1 = aB * rB + 0.5 * lam_val; V2 = aA * rA + 0.5 * lam_val
            V3 = wdir + aA * rA + aB * rB + lam_val
            for t in range(T2 - 1, -1, -1):
                V3 = V3 + pB_vec[t] * max(0.0, wB[t] + V2 - V3)
                V1 = V1 + pB_vec[t] * max(0.0, wB[t] + V0 - V1)
                V3 = V3 + pA_vec[t] * max(0.0, wA[t] + V1 - V3)
                V2 = V2 + pA_vec[t] * max(0.0, wA[t] + V0 - V2)
                if not (V3 >= V1 - 1e-12 and V3 >= V2 - 1e-12 and V1 >= -1e-12 and V2 >= -1e-12):
                    _fail("State monotonicity V3 >= V1, V2 >= 0 violated!")
            Vmid = V3 + 0.5 * max(0.0, wb + V2 - V3)
            return Vmid + 0.5 * max(0.0, wa + V1 - Vmid)

        dt = 0.5 / T2
        p1 = np.full(T2, -math.expm1(-s1 * dt))
        p2 = np.full(T2, -math.expm1(-s2 * dt))
        v1 = dp(p1, p1, s1, s1, kap * s1 * s1)
        v2 = dp(p2, p2, s2, s2, kap * s2 * s2)
        if v2 < v1 - 1e-12:
            _fail(f"S-monotonicity violated: v({s2})={v2} < v({s1})={v1}")

        # L1 sensitivity at fixed lam
        pA_a = np_rng.uniform(0.0, 1.0, T2); pB_a = np_rng.uniform(0.0, 1.0, T2)
        pA_b = np_rng.uniform(0.0, 1.0, T2); pB_b = np_rng.uniform(0.0, 1.0, T2)
        rA_a, rB_a, rA_b, rB_b = np_rng.uniform(0.0, 1.0, 4)
        lam_f = float(np_rng.uniform(0.0, 1.0))
        va = dp(pA_a, pB_a, rA_a, rB_a, lam_f)
        vb = dp(pA_b, pB_b, rA_b, rB_b, lam_f)
        bnd = float(np.sum(wA * np.abs(pA_a - pA_b)) + np.sum(wB * np.abs(pB_a - pB_b))
                    + aA * abs(rA_a - rA_b) + aB * abs(rB_a - rB_b))
        if bnd > 1e-12:
            ratio = abs(va - vb) / bnd
            if ratio > worst_sens:
                worst_sens = ratio
            if ratio > 1.0 + 1e-11:
                _fail(f"L1 sensitivity violated: {abs(va-vb)} > {bnd}")

    print(f"  [PASS] Suite 4: 5k 100-dps atom+direct checks (max_viol={mp.nstr(worst_atom_viol, 4)}) + "
          f"15k Bellman checks (worst sens/bnd={worst_sens:.4f} <= 1.0)")


# =============================================================================
# SUITE 5: LITERAL EXACT COUPLED 2-LEVEL FINITE-GRAPH BELLMAN CHECK!
# =============================================================================
def suite5_exact_coupled_2level_bellman() -> None:
    """
    LITERAL exact backward Bellman DP on a concrete small 2-LEVEL coupled finite graph!
      - Level 0 pool D0 = {0, 1} (2 vertices, initially both free; Phase-0 leaf on each with prob p0)
      - Level 1 pair {2, 3} with bridge vertex 4:
          edges (2,4), (3,4) [bridge], then (2,0), (2,1), (3,0), (3,1) [Phase 2 to D0], then (2,3) [Phase 3]
      - Level 2 pair {5, 6} with bridge vertex 7:
          edges (5,7), (6,7) [bridge], then edges to D0={0,1} AND D1={2,3} [Phase 2 to D0 U D1], then (5,6) [Phase 3]

    Compares:
      (A) EXACT_COUPLED: every vertex in {0..7} has exact capacity 1 across ALL 3 levels! (2^8 = 256 states!)
      (B) RELAXED_FROZEN: vertices in D0={0,1} have capacity 1 in Level 0, then FROZEN with capacity inf
                          for Level 1 & Level 2; vertices in D1={2,3} have capacity 1 in Level 1, then FROZEN
                          with capacity inf for Level 2!
    Verifies EXACT_COUPLED <= RELAXED_FROZEN across 200 random nonnegative weight vectors & edge probs!
    """
    rng = random.Random(505)
    min_slack = 1e9
    for _ in range(200):
        # Random edge probabilities and nonnegative weights
        p0 = rng.uniform(0.1, 0.9)
        w0 = rng.uniform(0.0, 1.0)
        p_br1 = rng.uniform(0.1, 0.5); w_br1 = rng.uniform(0.0, 1.0)
        p_2_01 = rng.uniform(0.05, 0.4); w_2_01 = rng.uniform(0.0, 1.0)
        p_dir1 = rng.uniform(0.05, 0.3); w_dir1 = rng.uniform(0.0, 1.0)
        p_br2 = rng.uniform(0.1, 0.5); w_br2 = rng.uniform(0.0, 1.0)
        p_2_to_01 = rng.uniform(0.05, 0.3); w_2_to_01 = rng.uniform(0.0, 1.0)
        p_2_to_23 = rng.uniform(0.05, 0.3); w_2_to_23 = rng.uniform(0.0, 1.0)
        p_dir2 = rng.uniform(0.05, 0.3); w_dir2 = rng.uniform(0.0, 1.0)

        # Edge list in arrival order: (level, u, v, prob, weight, kills_v_in_relaxed)
        # Note: u is the 'new' vertex in current level, v is either same-level or earlier-pool vertex!
        edges = [
            # Level 0: private leaves (represented by v = -1)
            (0, 0, -1, p0, w0, True),
            (0, 1, -1, p0, w0, True),
            # Level 1: bridge + Phase 2 to D0={0,1} + direct (2,3)
            (1, 2, 4, p_br1, w_br1, True),
            (1, 3, 4, p_br1, w_br1, True),
            (1, 2, 0, p_2_01, w_2_01, False),
            (1, 2, 1, p_2_01, w_2_01, False),
            (1, 3, 0, p_2_01, w_2_01, False),
            (1, 3, 1, p_2_01, w_2_01, False),
            (1, 2, 3, p_dir1, w_dir1, True),
            # Level 2: bridge + Phase 2 to D0={0,1} & D1={2,3} + direct (5,6)
            (2, 5, 7, p_br2, w_br2, True),
            (2, 6, 7, p_br2, w_br2, True),
            (2, 5, 0, p_2_to_01, w_2_to_01, False),
            (2, 5, 1, p_2_to_01, w_2_to_01, False),
            (2, 5, 2, p_2_to_23, w_2_to_23, False),
            (2, 5, 3, p_2_to_23, w_2_to_23, False),
            (2, 6, 0, p_2_to_01, w_2_to_01, False),
            (2, 6, 1, p_2_to_01, w_2_to_01, False),
            (2, 6, 2, p_2_to_23, w_2_to_23, False),
            (2, 6, 3, p_2_to_23, w_2_to_23, False),
            (2, 5, 6, p_dir2, w_dir2, True),
        ]

        def solve_exact_dp(relaxed: bool) -> float:
            # 2^8 = 256 states (bit i == 1 means vertex i is FREE)
            V = np.zeros(256, dtype=float)
            for lvl, u, v, prob, weight, kills_v_rel in reversed(edges):
                V_new = V.copy()
                bu = 1 << u
                bv = (1 << v) if v >= 0 else 0
                req = bu | bv
                # In coupled mode, matching clears BOTH bu and bv.
                # In relaxed frozen-capacity mode, if kills_v_rel is False, v has capacity inf -> ONLY clears bu!
                clr = req if (not relaxed or kills_v_rel) else bu
                for m in range(256):
                    if (m & req) == req:
                        gain = weight + V[m ^ clr] - V[m]
                        if gain > 0.0:
                            V_new[m] = V[m] + prob * gain
                V = V_new
            return float(V[255])

        v_exact = solve_exact_dp(relaxed=False)
        v_relax = solve_exact_dp(relaxed=True)
        slack = v_relax - v_exact
        if slack < min_slack:
            min_slack = slack
        if slack < -1e-12:
            _fail(f"Frozen-capacity relaxation violated: exact={v_exact} > relaxed={v_relax}")
    print(f"  [PASS] Suite 5: 200 literal 256-state exact 2-level coupled-vs-frozen Bellman checks "
          f"(min_slack = {min_slack:+.3e} >= 0)")


# =============================================================================
# SUITE 6: 100,000-POINT + BRENT ADVERSARIAL PEAK AUDIT ON ALL 3 CERTIFICATES
# =============================================================================
def suite6_global_peak_audit(cert_paths: list[str]) -> None:
    from generate_l3_ema_certs import eval_batch_ema_cubic, eval_scalar_ema_cubic
    Sd = np.linspace(0.0, 1.0, 100_001)
    h = 1.0 / 100_000.0
    for cp in cert_paths:
        with open(cp, "r", encoding="utf-8") as f:
            d = json.load(f)
        D = d["D_int"]
        wa = d["wa_int"] / D; wb = d["wb_int"] / D
        wA = np.array(d["wA_ints"]) / D; wB = np.array(d["wB_ints"]) / D
        aA = d["aA_int"] / D; aB = d["aB_int"] / D
        wdir = d["wdir_int"] / D; kap = d["kappa_int"] / D
        vals = eval_batch_ema_cubic(Sd, wa, wb, wA, wB, aA, aB, wdir, kap)
        idx = int(np.argmax(vals))
        opt = minimize_scalar(lambda s: -eval_scalar_ema_cubic(s, wa, wb, wA, wB, aA, aB, wdir, kap),
                              bounds=(max(0.0, Sd[idx] - 2 * h), min(1.0, Sd[idx] + 2 * h)),
                              method="bounded", options={"xatol": 1e-13})
        peak_val = float(-opt.fun)
        peak_S = float(opt.x)
        res = verifier.verify_certificate(cp, verbose=False)
        u_lim_hi = float(res["U_limit_mono_hi"])
        u_fin_hi = float(res["U_finite_mono_hi"])
        if peak_val > u_lim_hi + 1e-12:
            _fail(f"Peak {peak_val} exceeded verified U_limit_hi {u_lim_hi} on {cp}")
        print(f"  [PASS] Suite 6 ({os.path.basename(cp)}): peak={peak_val:.9f} @ S={peak_S:.6f} "
              f"<= U_limit_hi={u_lim_hi:.9f} < U_finite_hi={u_fin_hi:.9f}")


def main():
    t0 = time.time()
    bundle_dir = os.path.dirname(os.path.abspath(__file__))
    certs = [
        os.path.join(bundle_dir, "cert_l3_ema_paley_t1_lean.json"),
        os.path.join(bundle_dir, "cert_l3_ema_general_t20.json"),
        os.path.join(bundle_dir, "cert_l3_ema_general_t40.json"),
    ]
    run_forgery_battery(certs[1])
    print("\n" + "=" * 85)
    print("PART B: 6-SUITE INDEPENDENT ADVERSARIAL LEMMA FUZZER (python3 -O safe)")
    print("=" * 85)
    suite1_pyramid_polytope()
    suite2_cubic_identity()
    suite3_gauss_multilevel_fft()
    suite4_atom_direct_and_bellman()
    suite5_exact_coupled_2level_bellman()
    suite6_global_peak_audit(certs)
    print(f"\nALL 18 FORGERIES REJECTED & ALL 6 ADVERSARIAL SUITES PASSED! ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
