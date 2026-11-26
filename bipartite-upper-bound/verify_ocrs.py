#!/usr/bin/env python3
"""
Unified, Hardened Round-2 Matching OCRS Certificate Verifier (`verify_ocrs.py`).

Zero external LP/solver imports (`fractions`, `heapq`, `math`, `re`, `mpmath` only).
`python3 -O` safe: every check raises `VerificationError` via `require(...)`, never `assert`.

Supports all three canonical Round-2 certificate modes under schema `"ocrs.r2_canonical.v1"`:
  1. `general_bridged_fused`  : 4-family simple graph (pendant Phase-0 leaves + Phase-1 bridges
                                + T2 bulk Gauss-period waves + explicit t=1^- atom wave
                                + Phase-3A vanishing (w_i, D) K_{q,M} probe + Phase-3B direct).
                                Note: By Bellman translation invariance (V + C*1 on the state-11
                                branch where w_i survives Phase 1), Phase-3A (w_i, D) arrival and
                                post-Phase-1 arrival produce the exact same decoupled value
                                V3_p1 = V3 + u_wD * S, while Phase 2 only further depletes F in
                                the coupled finite graph.
  2. `general_bridged_atom`   : 3-family simple graph (omits (w_i, D) and Family 4; minimal
                                8-integer 293-box Lean-4 Paley architecture at T2=1, q=2,000,003).
  3. `bipartite_vanishing_atom`: 2-family Birkhoff bipartite simple graph (Phase-0 K_{M,M} + T2
                                bulk waves + explicit t=1^- atom wave + Phase-3 direct, verifying
                                both max_e x_e <= 10^-12 AND algebraic delta -> 0^+ vanishing
                                monotonicity for K(delta) = max(K, ceil(1/delta))).

Hardened Schema & Type Guards (incorporating all 5 Stage-2 audit fixes):
  - Exact mandatory `"schema": "ocrs.r2_canonical.v1"` check.
  - Mandatory `isinstance(claimed_ub, str)` (and `pure_q_claimed_ub` on Paley T2=1 certs) with
    strict decimal/fraction syntax `^\\d+(\\.\\d+|/\\d+)?$` (no float coercion or silent skip).
  - Exact key-set equality on both `cert_dict` and `raw_ints` (zero extra or missing keys).
  - Strict `K >= 2` across all modes so every edge satisfies `0 < x_e < 1`.
  - Exact rational `max_x_e` reporting (`num/den` and `1 - max_x_e`, never rounding `0.999999`
    to `1.000e+00`) plus first-principles 100-dps enclosures of Ma-MacRury-Nuti's beta & gamma.
"""
from __future__ import annotations

import argparse
import hashlib
import heapq
import json
import math
import re
import sys
from fractions import Fraction

import mpmath as mp

CANONICAL_SCHEMA = "ocrs.r2_canonical.v1"
STRICT_RAT_STR_RE = re.compile(r"^\d+(?:\.\d+|/\d+)?$")


class VerificationError(ValueError):
    """Raised whenever any certificate schema, algebraic, or interval check fails."""


def require(cond: bool, msg: str) -> None:
    if not cond:
        raise VerificationError(msg)


def is_strict_int(val: object) -> bool:
    return isinstance(val, int) and not isinstance(val, bool)


def parse_strict_rat_str(val: object, field_name: str) -> Fraction:
    require(isinstance(val, str), f"{field_name} must be a string, got {type(val).__name__}: {val!r}")
    s = str(val)
    require(bool(STRICT_RAT_STR_RE.match(s)), f"{field_name} has invalid rational/decimal syntax: {s!r}")
    if "/" in s:
        num_s, den_s = s.split("/", 1)
        den_i = int(den_s)
        require(den_i > 0, f"{field_name} denominator must be > 0")
        return Fraction(int(num_s), den_i)
    return Fraction(s)


def deterministic_is_prime(n: int) -> bool:
    """Deterministic 6k +/- 1 trial-division primality proof up to floor(sqrt(n))."""
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    limit = math.isqrt(n)
    k = 5
    while k <= limit:
        if n % k == 0 or n % (k + 2) == 0:
            return False
        k += 6
    return True


def verify_pratt_2000003() -> None:
    """Explicit zero-search Lucas/Pratt primality proof for q = 2,000,003 with g = 2."""
    tree = [
        (2_000_003, [2, 101, 9901], (2**1) * (101**1) * (9901**1)),
        (9901, [2, 3, 5, 11], (2**2) * (3**2) * (5**2) * 11),
        (101, [2, 5], (2**2) * (5**2)),
    ]
    for p, prime_divs, prod_check in tree:
        require(p - 1 == prod_check, f"Pratt factorization mismatch at {p}")
        require(pow(2, p - 1, p) == 1, f"Fermat check failed at {p}")
        for qd in prime_divs:
            require(pow(2, (p - 1) // qd, p) != 1, f"Order check failed at {p} for divisor {qd}")


def frac_to_iv(f: Fraction):
    return mp.iv.mpf(f.numerator) / mp.iv.mpf(f.denominator)


def iv_max(x, y):
    return mp.iv.mpf([max(x.a, y.a), max(x.b, y.b)])


def parse_and_validate_cert(cert_dict: object) -> dict:
    require(isinstance(cert_dict, dict), "Certificate root must be a JSON object")

    # Guard 1: Mandatory exact schema check
    schema = cert_dict.get("schema")
    require(isinstance(schema, str) and schema == CANONICAL_SCHEMA,
            f"Invalid or unsupported schema: {schema!r} (expected {CANONICAL_SCHEMA!r})")

    mode = cert_dict.get("mode")
    require(
        isinstance(mode, str) and mode in ("general_bridged_fused", "general_bridged_atom", "bipartite_vanishing_atom"),
        f"Unsupported mode: {mode!r}",
    )

    T2 = cert_dict.get("T2")
    require(is_strict_int(T2) and T2 >= 1, f"Invalid T2: {T2!r}")
    q = cert_dict.get("q")
    require(is_strict_int(q) and q > T2 + 1, f"Invalid prime q: {q!r}")

    is_paley_lean_eligible = (
        mode in ("general_bridged_fused", "general_bridged_atom") and T2 == 1 and (q % 4 == 3)
    )

    # Guard 3a: Strict exact root key-set check (no missing keys, no extra keys)
    base_keys = {"schema", "mode", "T2", "q", "K", "D_int", "raw_ints", "claimed_ub"}
    if mode == "general_bridged_fused":
        eps_keys = {"eps1_num", "eps1_den", "eps2_num", "eps2_den"}
    else:
        eps_keys = {"eps_num", "eps_den"}
    expected_root_keys = base_keys | eps_keys | ({"pure_q_claimed_ub"} if is_paley_lean_eligible else set())

    actual_root_keys = set(cert_dict.keys())
    missing_root = sorted(expected_root_keys - actual_root_keys)
    extra_root = sorted(actual_root_keys - expected_root_keys)
    require(not missing_root, f"Missing required root keys: {missing_root}")
    require(not extra_root, f"Unvalidated extra root keys forbidden: {extra_root}")

    # Guard 2: Mandatory strict string `claimed_ub` (and `pure_q_claimed_ub` when Paley T2=1)
    claimed_ub = parse_strict_rat_str(cert_dict["claimed_ub"], "claimed_ub")
    require(Fraction(0) < claimed_ub < Fraction(1), f"claimed_ub out of (0, 1): {claimed_ub}")
    pure_q_claimed_ub = (
        parse_strict_rat_str(cert_dict["pure_q_claimed_ub"], "pure_q_claimed_ub")
        if is_paley_lean_eligible
        else None
    )

    # Guard 4: Strict K >= 2 (`x_e > 0` on all edges) + deterministic primality & coset divisibility
    K = cert_dict["K"]
    require(is_strict_int(K) and K >= 2, f"Invalid K (must be integer >= 2 so x_e > 0): {K!r}")
    D_int = cert_dict["D_int"]
    require(is_strict_int(D_int) and D_int > 0, f"Invalid D_int: {D_int!r}")

    W_cosets = T2 + 1
    require((q - 1) % W_cosets == 0, f"q={q} fails cyclotomic coset condition q == 1 (mod {W_cosets})")
    require(deterministic_is_prime(q), f"q={q} failed deterministic primality proof")
    if q == 2_000_003:
        verify_pratt_2000003()

    # Guard 3b: Strict exact `raw_ints` key-set check
    raw = cert_dict["raw_ints"]
    require(isinstance(raw, dict), "raw_ints must be a JSON object")
    if mode == "general_bridged_fused":
        expected_raw_keys = {"w0", "wa", "wb", "wA", "wB", "aA", "aB", "u_wD", "wdir"}
    elif mode == "general_bridged_atom":
        expected_raw_keys = {"w0", "wa", "wb", "wA", "wB", "aA", "aB", "wdir"}
    else:
        expected_raw_keys = {"w0", "wA", "wB", "aA", "aB", "wdir"}

    actual_raw_keys = set(raw.keys())
    missing_raw = sorted(expected_raw_keys - actual_raw_keys)
    extra_raw = sorted(actual_raw_keys - expected_raw_keys)
    require(not missing_raw, f"Missing required raw_ints keys: {missing_raw}")
    require(not extra_raw, f"Unvalidated extra raw_ints keys forbidden: {extra_raw}")

    def get_int(k: str) -> int:
        v = raw[k]
        require(is_strict_int(v) and v >= 0, f"raw_ints[{k!r}] must be nonneg int, got {v!r}")
        return v

    def get_int_list(k: str) -> list[int]:
        arr = raw[k]
        require(isinstance(arr, list) and len(arr) == T2, f"raw_ints[{k!r}] must be list of length {T2}")
        out = []
        for i, v in enumerate(arr):
            require(is_strict_int(v) and v >= 0, f"raw_ints[{k!r}][{i}] must be nonneg int, got {v!r}")
            out.append(v)
        return out

    w0_i = get_int("w0")
    wA_i = get_int_list("wA")
    wB_i = get_int_list("wB")
    aA_i = get_int("aA")
    aB_i = get_int("aB")
    wdir_i = get_int("wdir")
    use_bridge = mode in ("general_bridged_fused", "general_bridged_atom")
    has_uwD = (mode == "general_bridged_fused")
    wa_i = get_int("wa") if use_bridge else 0
    wb_i = get_int("wb") if use_bridge else 0
    uwD_i = get_int("u_wD") if has_uwD else 0

    if use_bridge:
        sc = 2 * T2
        L_frac = Fraction(1, 2)
        expected_D = (
            sc * (w0_i + aA_i + aB_i + uwD_i + wdir_i)
            + T2 * (wa_i + wb_i)
            + sum(wA_i) + sum(wB_i)
        )
    else:
        sc = T2
        L_frac = Fraction(1, 1)
        expected_D = sc * (w0_i + aA_i + aB_i + wdir_i) + sum(wA_i) + sum(wB_i)

    require(D_int == expected_D, f"D_int mismatch: got {D_int}, expected {expected_D}")

    dt2 = L_frac / Fraction(T2)
    w0 = Fraction(sc * w0_i, D_int)
    wa = Fraction(sc * wa_i, D_int)
    wb = Fraction(sc * wb_i, D_int)
    wA = [Fraction(sc * v, D_int) for v in wA_i]
    wB = [Fraction(sc * v, D_int) for v in wB_i]
    aA = Fraction(sc * aA_i, D_int)
    aB = Fraction(sc * aB_i, D_int)
    u_wD = Fraction(sc * uwD_i, D_int)
    wdir = Fraction(sc * wdir_i, D_int)

    W_bulk = (sum(wA) + sum(wB)) * dt2
    norm_sum = w0 + (Fraction(1, 2) * (wa + wb) if use_bridge else Fraction(0)) + W_bulk + aA + aB + u_wD + wdir
    require(norm_sum == Fraction(1, 1), f"Exact rational normalization failed: {norm_sum} != 1")

    d0 = (q - 1) // W_cosets
    M = K * q
    vanishing_family_verified = False

    if mode == "general_bridged_fused":
        for nm in ("eps1_num", "eps1_den", "eps2_num", "eps2_den"):
            require(is_strict_int(cert_dict[nm]) and cert_dict[nm] > 0, f"Invalid {nm}: {cert_dict[nm]!r}")
        eps1 = Fraction(cert_dict["eps1_num"], cert_dict["eps1_den"])
        eps2 = Fraction(cert_dict["eps2_num"], cert_dict["eps2_den"])
        require(eps1 + 2 * eps2 < Fraction(1, 1), f"Need eps1 + 2*eps2 < 1, got {eps1 + 2 * eps2}")
        beta = Fraction(1, 1) - eps1 - 2 * eps2
        x0 = Fraction(1, 1) - Fraction(1, K)
        x_br = (Fraction(1, 1) - eps1) / 2
        x_bulk = beta / Fraction(2 * K * T2 * d0)
        x_atom = eps2 / Fraction(K * d0)
        x_wD = eps1 / Fraction(M)
        x_dir = eps1
        edge_probs = [x0, x_br, x_bulk, x_atom, x_wD, x_dir]
        require(all(Fraction(0) < x < Fraction(1) for x in edge_probs), "Edge marginal out of (0, 1)")
        # Verify exact 4-family matching polytope identities
        require(x_br + T2 * K * d0 * x_bulk + K * d0 * x_atom + x_dir == Fraction(1), "deg(a_i) != 1")
        require(2 * x_br + M * x_wD == Fraction(1), "deg(w_i) != 1")
        require(x0 + 2 * T2 * d0 * x_bulk + 2 * d0 * x_atom + q * x_wD == Fraction(1), "deg(d) != 1")
        require(2 * x_br + x_dir == Fraction(1), "3-blossom on (a_i, b_i, w_i) != 1")
        Z_exact = (
            Fraction(1, 1)
            - w0 / K
            - eps1 * (wa + wb) * Fraction(1, 2)
            - (eps1 + 2 * eps2) * W_bulk
        )
        eps_meta = {"eps1": str(eps1), "eps2": str(eps2)}
    elif mode == "general_bridged_atom":
        for nm in ("eps_num", "eps_den"):
            require(is_strict_int(cert_dict[nm]) and cert_dict[nm] > 0, f"Invalid {nm}: {cert_dict[nm]!r}")
        eps = Fraction(cert_dict["eps_num"], cert_dict["eps_den"])
        require(eps < Fraction(1, 4), f"Need eps < 1/4 in general_bridged_atom, got {eps}")
        x0 = Fraction(1, 1) - (Fraction(1, 1) - eps) / K
        x_br = (Fraction(1, 1) - eps) / 2
        x_bulk = (Fraction(1, 1) - 3 * eps) / Fraction(2 * T2 * K * d0)
        x_atom = eps / Fraction(K * d0)
        x_dir = eps
        edge_probs = [x0, x_br, x_bulk, x_atom, x_dir]
        require(all(Fraction(0) < x < Fraction(1) for x in edge_probs), "Edge marginal out of (0, 1)")
        require(x_br + T2 * K * d0 * x_bulk + K * d0 * x_atom + x_dir == Fraction(1), "deg(a_i) != 1")
        require(2 * x_br == Fraction(1) - eps, "deg(w_i) != 1 - eps")
        require(x0 + 2 * T2 * d0 * x_bulk + 2 * d0 * x_atom == Fraction(1), "deg(d) != 1")
        require(2 * x_br + x_dir == Fraction(1), "3-blossom != 1")
        Z_exact = (
            Fraction(1, 1)
            - (Fraction(1, 1) - eps) * w0 / K
            - eps * (wa + wb) * Fraction(1, 2)
            - 3 * eps * W_bulk
        )
        eps_meta = {"eps": str(eps)}
    else:
        for nm in ("eps_num", "eps_den"):
            require(is_strict_int(cert_dict[nm]) and cert_dict[nm] > 0, f"Invalid {nm}: {cert_dict[nm]!r}")
        eps = Fraction(cert_dict["eps_num"], cert_dict["eps_den"])
        require(eps < Fraction(1, 4), f"Need eps < 1/4, got {eps}")
        x0 = (Fraction(1, 1) - (Fraction(1, 1) - eps) / K) / Fraction(M)
        x_bulk = (Fraction(1, 1) - 2 * eps) / Fraction(T2 * K * d0)
        x_atom = eps / Fraction(K * d0)
        x_dir = eps
        edge_probs = [x0, x_bulk, x_atom, x_dir]
        require(all(Fraction(0) < x < Fraction(1) for x in edge_probs), "Edge marginal out of (0, 1)")
        require(T2 * K * d0 * x_bulk + K * d0 * x_atom + x_dir == Fraction(1), "deg(a_i) != 1")
        require(M * x0 + T2 * d0 * x_bulk + d0 * x_atom == Fraction(1), "deg(d) != 1")
        require((Fraction(1) - eps) * Fraction(K - 1, K) / M + eps / M == x0, "Birkhoff K_{M,M} marginal mismatch")
        max_xe_bip = max(edge_probs)
        require(max_xe_bip <= Fraction(1, 10**12),
                f"Vanishing violation: max_e x_e = {max_xe_bip} > 1/10^12")
        # Guard 5b: Algebraic delta -> 0^+ vanishing-family monotonicity proof:
        # Setting eps(delta) = delta and K(delta) = max(K, ceil(1/delta)) for any delta in (0, eps]:
        #   - 1/K(delta) <= delta <= eps <= (1 - eps)/K => (1 - delta)/K(delta) < delta <= (1 - eps)/K
        #   - Therefore Z(q, K(delta), delta) >= Z(q, K, eps) > 0 AND err_binom(K(delta)) <= err_binom(K),
        #     and every edge marginal in G_delta is <= delta -> 0^+.
        require(max_xe_bip == eps, "Expected Phase-3 eps edge to be the unique maximum marginal")
        require((Fraction(1, 1) - eps) / K >= eps,
                "Failed algebraic delta->0^+ monotonicity condition (1 - eps)/K >= eps")
        vanishing_family_verified = True
        Z_exact = Fraction(1, 1) - (Fraction(1, 1) - eps) * w0 / K - 2 * eps * W_bulk
        eps_meta = {"eps": str(eps)}

    require(Z_exact > Fraction(0, 1), f"Non-positive Z_exact: {Z_exact}")
    max_x_e = max(edge_probs)

    return {
        "mode": mode,
        "use_bridge": use_bridge,
        "has_uwD": has_uwD,
        "is_paley_lean_eligible": is_paley_lean_eligible,
        "T2": T2,
        "W_cosets": W_cosets,
        "L_frac": L_frac,
        "dt2": dt2,
        "q": q,
        "K": K,
        "w0": w0,
        "wa": wa,
        "wb": wb,
        "wA": wA,
        "wB": wB,
        "aA": aA,
        "aB": aB,
        "u_wD": u_wD,
        "wdir": wdir,
        "W_bulk": W_bulk,
        "Z_exact": Z_exact,
        "max_x_e": max_x_e,
        "vanishing_family_verified": vanishing_family_verified,
        "claimed_ub": claimed_ub,
        "pure_q_claimed_ub": pure_q_claimed_ub,
        "eps_meta": eps_meta,
    }


def verify_pure_rational_paley(p: dict) -> dict:
    """
    Zero-float, zero-exp, zero-Complex, zero-sqrt pure-`fractions.Fraction` verification
    for T2 = 1, q == 3 (mod 4) Paley certificates:
      - Uses Paley integer matrix identity bound err_paley <= (W_bulk + aA + aB) * ceil(sqrt(q+1)) / (2*(q-1))
      - Uses unconditional degree-5 Bonferroni polynomial P_5(S/2) >= 1 - exp(-S/2)
      - Computes BOTH monotone-split box count (e.g. 293 or 214 boxes) AND non-monotone box count (562 or 446 boxes).
    """
    q = p["q"]
    K = p["K"]
    w0, wa, wb = p["w0"], p["wa"], p["wb"]
    wA0, wB0 = p["wA"][0], p["wB"][0]
    aA, aB, u_wD, wdir = p["aA"], p["aB"], p["u_wD"], p["wdir"]
    W_bulk = p["W_bulk"]
    Z_q = p["Z_exact"]

    sqrt_ub = math.isqrt(q + 1)
    if sqrt_ub * sqrt_ub < q + 1:
        sqrt_ub += 1
    require(sqrt_ub * sqrt_ub >= q + 1, "Invalid integer ceiling for sqrt(q+1)")

    err_paley = (W_bulk + aA + aB) * Fraction(sqrt_ub, 2 * (q - 1))
    err_binom = W_bulk * Fraction(1, 2 * K * (q - 1))
    target_381 = Fraction(381, 1000)
    # Algebraic proof that 381/1000 < (3 - sqrt(5))/2: 3 - 2*(381/1000) = 1119/500 and (1119/500)^2 > 5
    require((Fraction(3) - 2 * target_381) ** 2 > Fraction(5), "Failed rational golden-ratio comparison")

    pure_q_claim = p["pure_q_claimed_ub"]
    require(pure_q_claim is not None and pure_q_claim <= target_381,
            f"pure_q_claimed_ub={pure_q_claim} must be <= 381/1000")

    def p5_pt(u: Fraction) -> Fraction:
        return u - (u**2) / 2 + (u**3) / 6 - (u**4) / 24 + (u**5) / 120

    def vp_pt(S: Fraction) -> Fraction:
        pr = p5_pt(S / 2)
        V1 = aB * S
        V2 = aA * S
        V3 = wdir + (aA + aB) * S
        V3 += pr * max(Fraction(0), wB0 + V2 - V3)
        V1 += pr * max(Fraction(0), wB0 - V1)
        V3 += pr * max(Fraction(0), wA0 + V1 - V3)
        V2 += pr * max(Fraction(0), wA0 - V2)
        V3_p1 = V3 + u_wD * S
        Vm = V3_p1 + Fraction(1, 2) * max(Fraction(0), wb + V2 - V3_p1)
        return Vm + Fraction(1, 2) * max(Fraction(0), wa + V1 - Vm)

    # 1. Monotone pure-Q adaptive bisection (exact Lean-4 `decide` box count)
    stk = [(Fraction(0), Fraction(1))]
    mono_boxes = []
    while stk:
        a, b = stk.pop()
        ub = (w0 * (Fraction(1) - a) + vp_pt(b) + err_paley + err_binom) / Z_q
        if ub < target_381 and ub <= pure_q_claim:
            mono_boxes.append((a, b, ub))
        else:
            mid = (a + b) / 2
            stk.extend([(mid, b), (a, mid)])
            require(len(mono_boxes) + len(stk) <= 5000, "Pure-Q monotone bisection exceeded 5000 boxes")
    mono_boxes.sort()
    require(mono_boxes[0][0] == 0 and mono_boxes[-1][1] == 1, "Pure-Q monotone partition boundary error")
    for idx in range(len(mono_boxes) - 1):
        require(mono_boxes[idx][1] == mono_boxes[idx + 1][0], "Pure-Q monotone partition gap")
    mono_ub = max(x[2] for x in mono_boxes)

    # 2. Non-monotone pure-Q adaptive bisection (zero Bellman S-monotonicity assumption;
    #    note P_5'(u) = (1 - u/2)^2 + (u^2/4)*(1 - 2u/3) + u^4/24 > 0 on all of R, so [P5(u0), P5(u1)] is exact!)
    def qadd(x, y): return (x[0] + y[0], x[1] + y[1])
    def qsub(x, y): return (x[0] - y[1], x[1] - y[0])
    def qmulp(c_pos: Fraction, x): return (c_pos * x[0], c_pos * x[1])
    def qmul(x, y):
        ps = [x[0] * y[0], x[0] * y[1], x[1] * y[0], x[1] * y[1]]
        return (min(ps), max(ps))
    def qmax0(x): return (max(Fraction(0), x[0]), max(Fraction(0), x[1]))
    def p5_iv(u): return (p5_pt(u[0]), p5_pt(u[1]))

    stk2 = [(Fraction(0), Fraction(1))]
    nm_boxes = []
    while stk2:
        a, b = stk2.pop()
        S_iv = (a, b)
        pr_iv = p5_iv((a / 2, b / 2))
        V1 = qmulp(aB, S_iv)
        V2 = qmulp(aA, S_iv)
        V3 = qadd((wdir, wdir), qmulp(aA + aB, S_iv))
        V3 = qadd(V3, qmul(pr_iv, qmax0(qsub(qadd((wB0, wB0), V2), V3))))
        V1 = qadd(V1, qmul(pr_iv, qmax0(qsub((wB0, wB0), V1))))
        V3 = qadd(V3, qmul(pr_iv, qmax0(qsub(qadd((wA0, wA0), V1), V3))))
        V2 = qadd(V2, qmul(pr_iv, qmax0(qsub((wA0, wA0), V2))))
        V3_p1 = qadd(V3, qmulp(u_wD, S_iv))
        Vm = qadd(V3_p1, qmulp(Fraction(1, 2), qmax0(qsub(qadd((wb, wb), V2), V3_p1))))
        Vp = qadd(Vm, qmulp(Fraction(1, 2), qmax0(qsub(qadd((wa, wa), V1), Vm))))
        F_iv = qadd(qmulp(w0, qsub((Fraction(1), Fraction(1)), S_iv)), Vp)
        ub_nm = (F_iv[1] + err_paley + err_binom) / Z_q
        if ub_nm < target_381:
            nm_boxes.append((a, b, ub_nm))
        else:
            mid = (a + b) / 2
            stk2.extend([(mid, b), (a, mid)])
            require(len(nm_boxes) + len(stk2) <= 5000, "Pure-Q non-monotone bisection exceeded 5000 boxes")
    nm_boxes.sort()
    require(nm_boxes[0][0] == 0 and nm_boxes[-1][1] == 1, "Pure-Q non-monotone boundary error")
    for idx in range(len(nm_boxes) - 1):
        require(nm_boxes[idx][1] == nm_boxes[idx + 1][0], "Pure-Q non-monotone partition gap")
    nm_ub = max(x[2] for x in nm_boxes)

    require(mono_ub <= pure_q_claim and mono_ub < target_381, "Pure-Q monotone bound exceeded claim")
    require(nm_ub < target_381, "Pure-Q non-monotone bound exceeded 381/1000")

    return {
        "pure_q_mono_boxes": len(mono_boxes),
        "pure_q_mono_ub": f"{float(mono_ub):.10f}",
        "pure_q_nonmono_boxes": len(nm_boxes),
        "pure_q_nonmono_ub": f"{float(nm_ub):.10f}",
    }


def verify_certificate(
    cert_dict: object,
    max_boxes: int = 12000,
    non_monotone: bool = False,
    non_monotone_boxes: int = 4096,
) -> dict:
    mp.mp.dps = 100
    mp.iv.dps = 100
    p = parse_and_validate_cert(cert_dict)

    T2 = p["T2"]
    w0_iv = frac_to_iv(p["w0"])
    wa_iv = frac_to_iv(p["wa"])
    wb_iv = frac_to_iv(p["wb"])
    wA_iv = [frac_to_iv(v) for v in p["wA"]]
    wB_iv = [frac_to_iv(v) for v in p["wB"]]
    aA_iv = frac_to_iv(p["aA"])
    aB_iv = frac_to_iv(p["aB"])
    uwD_iv = frac_to_iv(p["u_wD"])
    wdir_iv = frac_to_iv(p["wdir"])
    dt_iv = frac_to_iv(p["dt2"])
    half_iv = frac_to_iv(Fraction(1, 2))
    one_iv = mp.iv.mpf(1)
    zero_iv = mp.iv.mpf(0)

    def eval_vpair(S_iv):
        pr = one_iv - mp.iv.exp(-S_iv * dt_iv)
        qc = one_iv - pr
        V0 = zero_iv
        V1 = aB_iv * S_iv
        V2 = aA_iv * S_iv
        V3 = wdir_iv + (aA_iv + aB_iv) * S_iv
        for t in range(T2 - 1, -1, -1):
            V3 = qc * V3 + pr * iv_max(V3, wB_iv[t] + V2)
            V1 = qc * V1 + pr * iv_max(V1, wB_iv[t] + V0)
            V3 = qc * V3 + pr * iv_max(V3, wA_iv[t] + V1)
            V2 = qc * V2 + pr * iv_max(V2, wA_iv[t] + V0)
        if p["has_uwD"]:
            V3 = V3 + uwD_iv * S_iv
        if p["use_bridge"]:
            Vm = half_iv * V3 + half_iv * iv_max(V3, wb_iv + V2)
            return half_iv * Vm + half_iv * iv_max(Vm, wa_iv + V1)
        return V3

    def run_bnb(eval_fn, cap: int, tol: float):
        n0 = 64
        pq = []
        for k in range(n0):
            a = Fraction(k, n0)
            b = Fraction(k + 1, n0)
            iv = eval_fn(a, b)
            heapq.heappush(pq, (-float(iv.b), a, b, iv.a, iv.b))
        best_lo = max(it[3] for it in pq)
        used = n0
        while used < cap:
            neg_hi, a, b, lo, hi = heapq.heappop(pq)
            if hi - best_lo < mp.mpf(tol):
                heapq.heappush(pq, (neg_hi, a, b, lo, hi))
                break
            mid = (a + b) / 2
            for la, lb in ((a, mid), (mid, b)):
                iv = eval_fn(la, lb)
                used += 1
                if iv.a > best_lo:
                    best_lo = iv.a
                heapq.heappush(pq, (-float(iv.b), la, lb, iv.a, iv.b))
        ivs = sorted((it[1], it[2]) for it in pq)
        require(ivs[0][0] == Fraction(0) and ivs[-1][1] == Fraction(1), "Interval partition boundary error")
        for idx in range(len(ivs) - 1):
            require(ivs[idx][1] == ivs[idx + 1][0], "Interval partition gap detected")
        return best_lo, max(it[4] for it in pq), used

    def eval_mono_box(a: Fraction, b: Fraction):
        a_iv = frac_to_iv(a)
        b_iv = frac_to_iv(b)
        m_iv = frac_to_iv((a + b) / 2)
        hi_iv = w0_iv * (one_iv - a_iv) + eval_vpair(b_iv)
        mid_iv = w0_iv * (one_iv - m_iv) + eval_vpair(m_iv)
        return mp.iv.mpf([mid_iv.a, hi_iv.b])

    best_lo, U_limit_hi, boxes_used = run_bnb(eval_mono_box, cap=max_boxes, tol=5e-7)

    q_iv = mp.iv.mpf(p["q"])
    K_iv = mp.iv.mpf(p["K"])
    W_iv = mp.iv.mpf(p["W_cosets"])
    L_iv = frac_to_iv(p["L_frac"])
    Wb_iv = frac_to_iv(p["W_bulk"])
    Wtot_iv = frac_to_iv(p["W_bulk"] + p["aA"] + p["aB"])
    Z_iv = frac_to_iv(p["Z_exact"])

    err_gauss = Wtot_iv * (W_iv * mp.iv.sqrt(q_iv) / (mp.iv.mpf(2) * (q_iv - one_iv)))
    err_binom = Wb_iv * (L_iv * W_iv / (mp.iv.mpf(2 * T2) * K_iv * (q_iv - one_iv)))
    U_finite_iv = (mp.iv.mpf(U_limit_hi) + err_gauss + err_binom) / Z_iv
    U_finite_hi = U_finite_iv.b

    # First-principles 100-dps SOTA threshold checks!
    phi_iv = (mp.iv.mpf(3) - mp.iv.sqrt(mp.iv.mpf(5))) / mp.iv.mpf(2)
    # Prove beta > 3895380533 / 10^10 and gamma > 3929353184 / 10^10 directly in 100-dps interval arithmetic!
    beta_lo_f = Fraction(3895380533, 10**10)
    gamma_lo_f = Fraction(3929353184, 10**10)
    b_iv = frac_to_iv(beta_lo_f)
    g_iv = frac_to_iv(gamma_lo_f)
    require((one_iv - mp.iv.mpf(2) * b_iv - b_iv * mp.iv.exp(mp.iv.mpf(2) - one_iv / b_iv)).a > 0,
            "Failed 100-dps lower-bound proof for MMN beta")
    require((one_iv - mp.iv.mpf(2) * g_iv - g_iv * mp.iv.exp(g_iv - one_iv)).a > 0,
            "Failed 100-dps lower-bound proof for MMN gamma")

    if p["use_bridge"]:
        require(U_finite_hi < phi_iv.a,
                f"U_finite_hi={mp.nstr(mp.mpf(U_finite_hi), 15)} fails to beat (3-sqrt(5))/2")
        sota_margin_str = f"margin_below_golden_ratio={mp.nstr(mp.mpf(phi_iv.a) - mp.mpf(U_finite_hi), 10)}"
    else:
        require(U_finite_hi < b_iv.a and U_finite_hi < g_iv.a,
                f"U_finite_hi={mp.nstr(mp.mpf(U_finite_hi), 15)} fails to beat MMN beta/gamma")
        m_beta = mp.mpf(b_iv.a) - mp.mpf(U_finite_hi)
        m_gamma = mp.mpf(g_iv.a) - mp.mpf(U_finite_hi)
        sota_margin_str = (
            f"margin_below_beta(0.3895380533)={mp.nstr(m_beta, 10)}, "
            f"margin_below_gamma(0.3929353184)={mp.nstr(m_gamma, 10)}"
        )

    claimed_iv = frac_to_iv(p["claimed_ub"])
    require(U_finite_hi <= claimed_iv.a,
            f"Certified U_finite_hi={mp.nstr(mp.mpf(U_finite_hi), 16)} exceeds claimed_ub={p['claimed_ub']}")

    # Guard 5a: Exact rational max_x_e formatting (never misleadingly rounding 0.999999 to 1.000e+00)
    mx = p["max_x_e"]
    if mx <= Fraction(1, 1000):
        max_xe_display = f"{mx.numerator}/{mx.denominator} ({float(mx):.3e})"
    else:
        gap1 = Fraction(1, 1) - mx
        max_xe_display = f"{mx.numerator}/{mx.denominator} (1 - {float(gap1):.3e})"

    out = {
        "status": "VERIFIED",
        "schema": CANONICAL_SCHEMA,
        "mode": p["mode"],
        "T2": T2,
        "q": p["q"],
        "K": p["K"],
        **p["eps_meta"],
        "max_x_e": max_xe_display,
        "vanishing_delta_to_0_verified": p["vanishing_family_verified"],
        "boxes_mono": boxes_used,
        "U_limit_lo": mp.nstr(mp.mpf(best_lo), 15),
        "U_limit_hi": mp.nstr(mp.mpf(U_limit_hi), 15),
        "err_gauss": mp.nstr(mp.mpf(err_gauss.b), 10),
        "err_binom": mp.nstr(mp.mpf(err_binom.b), 10),
        "Z_exact": f"{float(p['Z_exact']):.12f}",
        "U_finite_hi": mp.nstr(mp.mpf(U_finite_hi), 16),
        "claimed_ub": cert_dict["claimed_ub"],
        "sota_check": sota_margin_str,
    }

    if p["is_paley_lean_eligible"]:
        err_paley_iv = Wtot_iv * (mp.iv.sqrt(q_iv + one_iv) / (mp.iv.mpf(2) * (q_iv - one_iv)))
        u_paley_iv = (mp.iv.mpf(U_limit_hi) + err_paley_iv + err_binom) / Z_iv
        out["U_finite_paley_tight_hi"] = mp.nstr(mp.mpf(u_paley_iv.b), 16)
        out.update(verify_pure_rational_paley(p))

    if non_monotone:
        def eval_nm_box(a: Fraction, b: Fraction):
            a_iv = frac_to_iv(a)
            b_iv = frac_to_iv(b)
            m_iv = frac_to_iv((a + b) / 2)
            S_box = mp.iv.mpf([a_iv.a, b_iv.b])
            f_box = w0_iv * (one_iv - S_box) + eval_vpair(S_box)
            f_mid = w0_iv * (one_iv - m_iv) + eval_vpair(m_iv)
            return mp.iv.mpf([f_mid.a, f_box.b])

        _, U_nm_lim_hi, nm_used = run_bnb(eval_nm_box, cap=non_monotone_boxes, tol=1e-6)
        U_nm_fin_hi = ((mp.iv.mpf(U_nm_lim_hi) + err_gauss + err_binom) / Z_iv).b
        if p["use_bridge"]:
            require(U_nm_fin_hi < phi_iv.a, "Non-monotone 100-dps bound failed to beat golden ratio")
        else:
            require(U_nm_fin_hi < b_iv.a and U_nm_fin_hi < g_iv.a, "Non-monotone 100-dps bound failed to beat MMN beta/gamma")
        out["boxes_nonmono"] = nm_used
        out["U_finite_nonmono_hi"] = mp.nstr(mp.mpf(U_nm_fin_hi), 16)

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Unified Hardened Round-2 Matching OCRS Certificate Verifier")
    ap.add_argument("certs", nargs="+", help="Paths to canonical JSON certificate files")
    ap.add_argument("--non-monotone", action="store_true", help="Also run 100-dps zero-monotonicity interval verification")
    ap.add_argument("--max-boxes", type=int, default=12000, help="Max monotone rational partition boxes (default: 12000)")
    ap.add_argument("--nonmono-boxes", type=int, default=4096, help="Max non-monotone partition boxes (default: 4096)")
    args = ap.parse_args()

    for path in args.certs:
        with open(path, "rb") as f:
            raw_bytes = f.read()
        sha = hashlib.sha256(raw_bytes).hexdigest()
        cert = json.loads(raw_bytes.decode("utf-8"))
        res = verify_certificate(
            cert,
            max_boxes=args.max_boxes,
            non_monotone=args.non_monotone,
            non_monotone_boxes=args.nonmono_boxes,
        )
        print(f"PASS [{path}] (sha256={sha[:16]}...)")
        for k, v in res.items():
            print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
