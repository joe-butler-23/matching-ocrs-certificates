#!/usr/bin/env python3
"""
STANDALONE ZERO-SHARED-CODE 100-DPS OUTWARD-ROUNDED VERIFIER
for Round-3 Level-3 EMA Cubic-Potential General Matching OCRS Certificates!

Hard Zero-Approximation Verification Standard (python3 -O safe, ZERO asserts for checks):
  1. Strict schema ("ocrs_l3_ema_cubic_v1"), exact allowed key set (no missing or extra keys),
     mandatory exact-string "claimed_ub", strict int (no bool/float) checks everywhere.
  2. Deterministic primality of q via exact trial division up to floor(sqrt(q)),
     exact congruence q == 1 (mod T2 + 1), strict P >= 3, K0 >= 2, N_L >= 2, L >= 2,
     strict beta = 1/beta_den in (0, 1/2], eps1 + 2*eps2 < 1, exact Fraction normalization == 1.
  3. 100-dps outward-rounded mpmath.iv evaluation of:
       F(S) = V_Bellman(S; w, +kappa * S^2) - kappa * S^3
     using BOTH:
       (a) Monotone exact rational box splitting on a gap-free partition of [0, 1]:
             F([a, b]) in [V(a_iv).a - kappa*b_iv^3.b,  V(b_iv).b - kappa*a_iv^3.a]
       (b) Independent RAW NON-MONOTONE adaptive branch-and-bound mpmath.iv pass
           evaluating V_Bellman([a, b]; w, +kappa*[a, b]^2) - kappa*[a, b]^3 directly
           with ZERO monotonicity assumption!
  4. Strictly outward-rounded finite-graph error terms:
       err_gauss     = (W_p2 + aA + aB) * (T2 + 1) * sqrt(q) / (2 * (q - 1))
       err_binom     = W_p2 * (T2 + 1) / (8 * T2 * N_L * (q - 1))
       err_cubic     = kappa * beta
       err_telescope = kappa / (3 * beta * L)
       Z_exact       = (1 - beta / (P - 1)) * (1 - eps1*(wa+wb)/2 - (eps1 + 2*eps2)*W_p2)
     Proves BOTH U_finite_hi_mono <= claimed_ub < (3 - sqrt(5))/2
     AND        U_finite_hi_raw  <= claimed_ub < (3 - sqrt(5))/2 !
"""
from __future__ import annotations

import json
import math
import sys
import time
from fractions import Fraction
import mpmath as mp


EXPECTED_SCHEMA = "ocrs_l3_ema_cubic_v1"
REQUIRED_KEYS = frozenset({
    "schema", "T2", "q", "P", "K0", "N_L", "L",
    "beta_den", "eps1_den", "eps2_den",
    "D_int", "wa_int", "wb_int", "wA_ints", "wB_ints",
    "aA_int", "aB_int", "wdir_int", "kappa_int", "claimed_ub",
})


def _fail(msg: str) -> None:
    raise ValueError(f"CERTIFICATE REJECTED: {msg}")


def _req_int(obj: dict, key: str, min_val: int = 0) -> int:
    v = obj.get(key)
    if isinstance(v, bool) or not isinstance(v, int):
        _fail(f"'{key}' must be an exact int (got {type(v).__name__}: {v!r})")
    if v < min_val:
        _fail(f"'{key}' = {v} < minimum {min_val}")
    return v


def _req_int_list(obj: dict, key: str, expected_len: int) -> list[int]:
    lst = obj.get(key)
    if not isinstance(lst, list) or len(lst) != expected_len:
        _fail(f"'{key}' must be a list of length {expected_len}")
    out = []
    for idx, v in enumerate(lst):
        if isinstance(v, bool) or not isinstance(v, int) or v < 0:
            _fail(f"'{key}[{idx}]' must be a non-negative exact int (got {v!r})")
        out.append(v)
    return out


def is_prime_deterministic(n: int) -> bool:
    """100% deterministic trial-division primality proof up to isqrt(n)."""
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    lim = math.isqrt(n)
    d = 5
    while d <= lim:
        if n % d == 0 or n % (d + 2) == 0:
            return False
        d += 6
    return True


def frac_to_iv(fr: Fraction) -> mp.iv.mpf:
    """Strictly outward-rounded 100-dps interval enclosure of exact Fraction."""
    return mp.iv.mpf(fr.numerator) / mp.iv.mpf(fr.denominator)


def iv_max0(x: mp.iv.mpf) -> mp.iv.mpf:
    """Exact outward-rounded enclosure of max(0, x) for interval x."""
    xa = mp.mpf(x.a)
    xb = mp.mpf(x.b)
    lo = xa if xa > 0 else mp.mpf(0)
    hi = xb if xb > 0 else mp.mpf(0)
    return mp.iv.mpf([lo, hi])


def eval_V_bellman_iv(
    S_iv: mp.iv.mpf,
    dt_iv: mp.iv.mpf,
    wa_iv: mp.iv.mpf,
    wb_iv: mp.iv.mpf,
    wA_ivs: list[mp.iv.mpf],
    wB_ivs: list[mp.iv.mpf],
    aA_iv: mp.iv.mpf,
    aB_iv: mp.iv.mpf,
    wdir_iv: mp.iv.mpf,
    kappa_iv: mp.iv.mpf,
) -> mp.iv.mpf:
    """
    Outward-rounded 100-dps interval evaluation of V_Bellman(S; w, +kappa*S^2).
    Valid for BOTH point intervals S_iv = [s, s] AND wide intervals S_iv = [a, b]!
    """
    half = mp.iv.mpf(1) / mp.iv.mpf(2)
    one = mp.iv.mpf(1)
    p_iv = one - mp.iv.exp(-S_iv * dt_iv)
    lam_iv = kappa_iv * S_iv * S_iv
    half_lam = half * lam_iv

    V0 = mp.iv.mpf(0)
    V1 = aB_iv * S_iv + half_lam
    V2 = aA_iv * S_iv + half_lam
    V3 = wdir_iv + (aA_iv + aB_iv) * S_iv + lam_iv

    T2 = len(wA_ivs)
    for t in range(T2 - 1, -1, -1):
        V3 = V3 + p_iv * iv_max0(wB_ivs[t] + V2 - V3)
        V1 = V1 + p_iv * iv_max0(wB_ivs[t] + V0 - V1)
        V3 = V3 + p_iv * iv_max0(wA_ivs[t] + V1 - V3)
        V2 = V2 + p_iv * iv_max0(wA_ivs[t] + V0 - V2)

    Vmid = V3 + half * iv_max0(wb_iv + V2 - V3)
    Vpair = Vmid + half * iv_max0(wa_iv + V1 - Vmid)
    return Vpair


def eval_F_monotone_box(
    a_fr: Fraction, b_fr: Fraction, params: tuple
) -> tuple[mp.mpf, mp.mpf]:
    """
    Monotone exact splitting on [a, b]:
      V_Bellman(S; w, +kappa*S^2) is provably non-decreasing in S (since w >= 0, kappa >= 0),
      -kappa * S^3 is provably non-increasing in S.
    Therefore F([a, b]) is enclosed in [V(a).a - kappa*b^3.b,  V(b).b - kappa*a^3.a].
    """
    dt_iv, wa_iv, wb_iv, wA_ivs, wB_ivs, aA_iv, aB_iv, wdir_iv, kappa_iv = params
    a_iv = frac_to_iv(a_fr)
    b_iv = frac_to_iv(b_fr)
    Va = eval_V_bellman_iv(a_iv, dt_iv, wa_iv, wb_iv, wA_ivs, wB_ivs, aA_iv, aB_iv, wdir_iv, kappa_iv)
    Vb = eval_V_bellman_iv(b_iv, dt_iv, wa_iv, wb_iv, wA_ivs, wB_ivs, aA_iv, aB_iv, wdir_iv, kappa_iv)
    a3 = kappa_iv * a_iv * a_iv * a_iv
    b3 = kappa_iv * b_iv * b_iv * b_iv
    lo_iv = mp.iv.mpf(Va.a) - mp.iv.mpf(b3.b)
    hi_iv = mp.iv.mpf(Vb.b) - mp.iv.mpf(a3.a)
    return mp.mpf(lo_iv.a), mp.mpf(hi_iv.b)


def eval_F_nonmonotone_box(
    a_fr: Fraction, b_fr: Fraction, params: tuple
) -> tuple[mp.mpf, mp.mpf]:
    """
    100% assumption-free RAW interval evaluation of F([a, b]) = V([a, b]) - kappa*[a, b]^3!
    Uses ZERO monotonicity assumption — pure outward-rounded mpmath.iv arithmetic!
    """
    dt_iv, wa_iv, wb_iv, wA_ivs, wB_ivs, aA_iv, aB_iv, wdir_iv, kappa_iv = params
    a_iv = frac_to_iv(a_fr)
    b_iv = frac_to_iv(b_fr)
    S_box = mp.iv.mpf([mp.mpf(a_iv.a), mp.mpf(b_iv.b)])
    V_box = eval_V_bellman_iv(S_box, dt_iv, wa_iv, wb_iv, wA_ivs, wB_ivs, aA_iv, aB_iv, wdir_iv, kappa_iv)
    F_box = V_box - kappa_iv * S_box * S_box * S_box
    return mp.mpf(F_box.a), mp.mpf(F_box.b)


def verify_certificate(cert_path: str, verbose: bool = True) -> dict:
    t0 = time.time()
    mp.mp.dps = 100
    mp.iv.dps = 100

    with open(cert_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        _fail("JSON root must be an object")

    # 1. Exact schema & key-set check
    if data.get("schema") != EXPECTED_SCHEMA:
        _fail(f"Invalid schema: expected {EXPECTED_SCHEMA!r}, got {data.get('schema')!r}")
    keys = set(data.keys())
    if keys != REQUIRED_KEYS:
        missing = sorted(REQUIRED_KEYS - keys)
        extra = sorted(keys - REQUIRED_KEYS)
        _fail(f"Exact key-set mismatch: missing={missing}, extra={extra}")

    # Mandatory exact string claimed_ub
    claimed_str = data["claimed_ub"]
    if not isinstance(claimed_str, str) or not claimed_str.strip():
        _fail(f"'claimed_ub' must be a non-empty string (got {type(claimed_str).__name__})")
    try:
        claimed_fr = Fraction(claimed_str)
    except Exception as exc:
        _fail(f"'claimed_ub' is not a valid exact rational string: {claimed_str!r} ({exc})")
    if claimed_fr <= 0 or claimed_fr >= 1:
        _fail(f"'claimed_ub' = {claimed_fr} out of (0, 1)")

    # 2. Strict integer parameter validation
    T2 = _req_int(data, "T2", min_val=1)
    q = _req_int(data, "q", min_val=3)
    P = _req_int(data, "P", min_val=3)
    K0 = _req_int(data, "K0", min_val=2)
    N_L = _req_int(data, "N_L", min_val=2)
    L = _req_int(data, "L", min_val=2)
    beta_den = _req_int(data, "beta_den", min_val=2)
    eps1_den = _req_int(data, "eps1_den", min_val=4)
    eps2_den = _req_int(data, "eps2_den", min_val=4)

    D_int = _req_int(data, "D_int", min_val=1)
    wa_int = _req_int(data, "wa_int", min_val=0)
    wb_int = _req_int(data, "wb_int", min_val=0)
    wA_ints = _req_int_list(data, "wA_ints", T2)
    wB_ints = _req_int_list(data, "wB_ints", T2)
    aA_int = _req_int(data, "aA_int", min_val=0)
    aB_int = _req_int(data, "aB_int", min_val=0)
    wdir_int = _req_int(data, "wdir_int", min_val=0)
    kappa_int = _req_int(data, "kappa_int", min_val=0)

    # Exact primality & congruence
    if (q - 1) % (T2 + 1) != 0:
        _fail(f"q - 1 = {q - 1} is NOT divisible by T2 + 1 = {T2 + 1}")
    if not is_prime_deterministic(q):
        _fail(f"q = {q} is COMPOSITE (failed deterministic trial division)")

    # Exact rational polytope feasibility checks
    beta_fr = Fraction(1, beta_den)
    eps1_fr = Fraction(1, eps1_den)
    eps2_fr = Fraction(1, eps2_den)
    if eps1_fr + 2 * eps2_fr >= Fraction(1):
        _fail(f"eps1 + 2*eps2 = {eps1_fr + 2*eps2_fr} >= 1 (invalid matching polytope)")

    wa_fr = Fraction(wa_int, D_int)
    wb_fr = Fraction(wb_int, D_int)
    wA_frs = [Fraction(v, D_int) for v in wA_ints]
    wB_frs = [Fraction(v, D_int) for v in wB_ints]
    aA_fr = Fraction(aA_int, D_int)
    aB_fr = Fraction(aB_int, D_int)
    wdir_fr = Fraction(wdir_int, D_int)
    kappa_fr = Fraction(kappa_int, D_int)

    W_p2_fr = sum(wA_frs + wB_frs, Fraction(0)) / Fraction(2 * T2)
    norm_fr = (wa_fr + wb_fr) / Fraction(2) + W_p2_fr + aA_fr + aB_fr + wdir_fr
    if norm_fr != Fraction(1):
        _fail(f"Exact rational normalization = {norm_fr} != 1")

    # 3. Finite-graph error terms & denominator (100-dps outward-rounded!)
    q_iv = mp.iv.mpf(q)
    qm1_iv = mp.iv.mpf(q - 1)
    T2_iv = mp.iv.mpf(T2)
    W_cos_iv = mp.iv.mpf(T2 + 1)
    W_p2_iv = frac_to_iv(W_p2_fr)
    aA_iv = frac_to_iv(aA_fr)
    aB_iv = frac_to_iv(aB_fr)
    wa_iv = frac_to_iv(wa_fr)
    wb_iv = frac_to_iv(wb_fr)
    wA_ivs = [frac_to_iv(v) for v in wA_frs]
    wB_ivs = [frac_to_iv(v) for v in wB_frs]
    wdir_iv = frac_to_iv(wdir_fr)
    kappa_iv = frac_to_iv(kappa_fr)
    beta_iv = frac_to_iv(beta_fr)
    eps1_iv = frac_to_iv(eps1_fr)
    eps2_iv = frac_to_iv(eps2_fr)
    dt_iv = frac_to_iv(Fraction(1, 2 * T2))

    # Gauss-period expander error: (W_p2 + aA + aB) * (T2 + 1) * sqrt(q) / (2*(q - 1))
    err_gauss_iv = (W_p2_iv + aA_iv + aB_iv) * W_cos_iv * mp.iv.sqrt(q_iv) / (mp.iv.mpf(2) * qm1_iv)
    # Binomial-vs-Poisson error: W_p2 * (T2 + 1) / (8 * T2 * N_L * (q - 1))
    err_binom_iv = W_p2_iv * W_cos_iv / (mp.iv.mpf(8) * T2_iv * mp.iv.mpf(N_L) * qm1_iv)
    # Exact 1-line cubic polynomial remainder: kappa * beta
    err_cubic_iv = kappa_iv * beta_iv
    # Exact telescoping boundary error: kappa / (3 * beta * L)
    err_tele_iv = kappa_iv / (mp.iv.mpf(3) * beta_iv * mp.iv.mpf(L))
    total_err_iv = err_gauss_iv + err_binom_iv + err_cubic_iv + err_tele_iv

    # Exact denominator lower bound: (1 - beta/(P-1)) * (1 - eps1*(wa+wb)/2 - (eps1+2*eps2)*W_p2)
    gamma_lb_fr = Fraction(1) - beta_fr / Fraction(P - 1)
    Z_lvl_fr = Fraction(1) - eps1_fr * (wa_fr + wb_fr) / Fraction(2) - (eps1_fr + 2 * eps2_fr) * W_p2_fr
    Z_fr = gamma_lb_fr * Z_lvl_fr
    if Z_fr <= 0:
        _fail(f"Denominator Z = {Z_fr} <= 0")
    Z_iv = frac_to_iv(Z_fr)

    claimed_iv = frac_to_iv(claimed_fr)
    golden_iv = (mp.iv.mpf(3) - mp.iv.sqrt(mp.iv.mpf(5))) / mp.iv.mpf(2)

    # Target limit threshold: we need (U_limit_hi + total_err_iv.b) / Z_iv.a <= claimed_iv.a
    # => U_limit_hi <= (claimed_iv * Z_iv - total_err_iv).a  (strictly conservative lower bound!)
    u_lim_target_iv = claimed_iv * Z_iv - total_err_iv
    u_lim_target = mp.mpf(u_lim_target_iv.a)
    if u_lim_target <= 0:
        _fail("Finite error terms exceed claimed_ub * Z")

    params = (dt_iv, wa_iv, wb_iv, wA_ivs, wB_ivs, aA_iv, aB_iv, wdir_iv, kappa_iv)

    # 4A. PASS 1: Monotone gap-free adaptive rational box verification over [0, 1]
    # Start with 256 uniform exact rational boxes, adaptively refine highest-upper-bound boxes!
    init_n = 256
    boxes_mono = [(Fraction(k, init_n), Fraction(k + 1, init_n)) for k in range(init_n)]
    mono_evals = []
    for a, b in boxes_mono:
        lo, hi = eval_F_monotone_box(a, b, params)
        mono_evals.append([hi, lo, a, b])

    # Refine worst boxes up to 6,000 total boxes (adaptive branch-and-bound on peak!)
    max_boxes_mono = 6000
    while len(mono_evals) < max_boxes_mono:
        mono_evals.sort(key=lambda r: r[0], reverse=True)
        worst_hi, worst_lo, a, b = mono_evals[0]
        # Stop early once global width < 1e-7 AND well below target
        max_lo = max(r[1] for r in mono_evals)
        if (worst_hi - max_lo) < mp.mpf("2e-7") and worst_hi <= u_lim_target:
            break
        # Refine top 32 boxes at once for speed
        batch = min(32, (max_boxes_mono - len(mono_evals)) // 1)
        next_new = []
        for idx in range(batch):
            if idx >= len(mono_evals):
                break
            hi_k, lo_k, ak, bk = mono_evals[idx]
            if hi_k <= max_lo and hi_k <= u_lim_target:
                next_new.append(mono_evals[idx])
                continue
            mid = (ak + bk) / Fraction(2)
            l1, h1 = eval_F_monotone_box(ak, mid, params)
            l2, h2 = eval_F_monotone_box(mid, bk, params)
            next_new.append([h1, l1, ak, mid])
            next_new.append([h2, l2, mid, bk])
        mono_evals = next_new + mono_evals[batch:]

    # Verify exact gap-free partition of [0, 1] for Pass 1
    mono_sorted = sorted(mono_evals, key=lambda r: r[2])
    if mono_sorted[0][2] != Fraction(0) or mono_sorted[-1][3] != Fraction(1):
        _fail("Monotone partition does not cover [0, 1]")
    for k in range(len(mono_sorted) - 1):
        if mono_sorted[k][3] != mono_sorted[k + 1][2]:
            _fail(f"Gap/overlap in monotone partition at box {k}")
    U_lim_mono_hi = max(r[0] for r in mono_sorted)
    U_lim_mono_lo = max(r[1] for r in mono_sorted)
    U_fin_mono_iv = (mp.iv.mpf([U_lim_mono_lo, U_lim_mono_hi]) + total_err_iv) / Z_iv

    # Fast-fail monotone check BEFORE running Pass 2 raw non-monotone pass
    u_fin_mono_hi = mp.mpf(U_fin_mono_iv.b)
    claimed_lo = mp.mpf(claimed_iv.a)
    golden_lo = mp.mpf(golden_iv.a)
    if not (u_fin_mono_hi <= claimed_lo):
        _fail(f"Monotone U_finite_hi = {mp.nstr(u_fin_mono_hi, 15)} > claimed_ub = {claimed_str}")
    if not (u_fin_mono_hi < golden_lo):
        _fail(f"Does NOT beat golden ratio (3-sqrt(5))/2: U_fin_hi = {mp.nstr(u_fin_mono_hi, 15)}")

    # 4B. PASS 2: 100% Assumption-Free RAW NON-MONOTONE Interval Cross-Check over [0, 1]!
    # Adaptive branch-and-bound: split ONLY boxes with raw_hi > u_lim_target!
    init_raw = 128
    raw_stack = [(Fraction(k, init_raw), Fraction(k + 1, init_raw)) for k in range(init_raw)]
    raw_done = []
    max_boxes_raw = 25000
    raw_evals_count = 0
    while raw_stack:
        a, b = raw_stack.pop()
        lo, hi = eval_F_nonmonotone_box(a, b, params)
        raw_evals_count += 1
        if hi <= u_lim_target:
            raw_done.append((hi, lo, a, b))
        else:
            if len(raw_done) + len(raw_stack) >= max_boxes_raw or (b - a) < Fraction(1, 10**12):
                raw_done.append((hi, lo, a, b))
                continue
            mid = (a + b) / Fraction(2)
            raw_stack.append((mid, b))
            raw_stack.append((a, mid))

    # Verify exact gap-free partition of [0, 1] for Pass 2
    raw_done.sort(key=lambda r: r[2])
    if raw_done[0][2] != Fraction(0) or raw_done[-1][3] != Fraction(1):
        _fail("Non-monotone partition does not cover [0, 1]")
    for k in range(len(raw_done) - 1):
        if raw_done[k][3] != raw_done[k + 1][2]:
            _fail(f"Gap/overlap in non-monotone partition at box {k}")
    U_lim_raw_hi = max(r[0] for r in raw_done)
    U_lim_raw_lo = max(r[1] for r in raw_done)
    U_fin_raw_iv = (mp.iv.mpf([U_lim_raw_lo, U_lim_raw_hi]) + total_err_iv) / Z_iv

    # 5. Final Hard Checks (BOTH monotone AND raw non-monotone MUST pass!)
    u_fin_raw_hi = mp.mpf(U_fin_raw_iv.b)
    if not (u_fin_raw_hi <= claimed_lo):
        _fail(f"Raw non-monotone U_finite_hi = {mp.nstr(u_fin_raw_hi, 15)} > claimed_ub = {claimed_str}")
    if not (u_fin_raw_hi < golden_lo):
        _fail(f"Raw non-monotone does NOT beat golden ratio (3-sqrt(5))/2: U_fin_raw_hi = {mp.nstr(u_fin_raw_hi, 15)}")

    margin_golden = mp.mpf((golden_iv - U_fin_mono_iv).a)
    elapsed = time.time() - t0
    res = {
        "file": cert_path,
        "T2": T2,
        "q": q,
        "mono_boxes": len(mono_sorted),
        "raw_boxes": len(raw_done),
        "U_limit_mono_hi": mp.nstr(U_lim_mono_hi, 15),
        "U_finite_mono_hi": mp.nstr(u_fin_mono_hi, 15),
        "U_finite_raw_hi": mp.nstr(u_fin_raw_hi, 15),
        "claimed_ub": claimed_str,
        "err_gauss": mp.nstr(mp.mpf(err_gauss_iv.b), 6),
        "err_binom": mp.nstr(mp.mpf(err_binom_iv.b), 6),
        "err_cubic": mp.nstr(mp.mpf(err_cubic_iv.b), 6),
        "err_tele": mp.nstr(mp.mpf(err_tele_iv.b), 6),
        "Z_lo": mp.nstr(mp.mpf(Z_iv.a), 12),
        "margin_below_golden": mp.nstr(margin_golden, 12),
        "elapsed_s": round(elapsed, 2),
    }
    if verbose:
        print(f"[VERIFIED] {cert_path}")
        print(f"  T2={T2}, W_cosets={T2+1}, q={q} (prime == 1 mod {T2+1}), P={P}, K0={K0}, N_L={N_L}, L={L}")
        print(f"  Monotone Pass ({len(mono_sorted)} boxes):     U_limit_hi = {res['U_limit_mono_hi']} | U_finite_hi = {res['U_finite_mono_hi']}")
        print(f"  Raw Non-Mono Pass ({len(raw_done)} boxes): U_finite_hi = {res['U_finite_raw_hi']} <= claimed_ub = {claimed_str}")
        print(f"  Errors: gauss={res['err_gauss']}, binom={res['err_binom']}, cubic={res['err_cubic']}, tele={res['err_tele']} | Z_lo={res['Z_lo']}")
        print(f"  Margin below (3-sqrt(5))/2: +{res['margin_below_golden']}  ({elapsed:.2f}s)")
    return res


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 -O verify_l3_ema_ocrs.py <cert1.json> [cert2.json ...]")
        sys.exit(2)
    for path in sys.argv[1:]:
        verify_certificate(path, verbose=True)


if __name__ == "__main__":
    main()
