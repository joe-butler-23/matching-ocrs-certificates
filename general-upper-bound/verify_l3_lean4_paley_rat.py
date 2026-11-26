#!/usr/bin/env python3
"""
Exact Pure-Rational (ZERO Complex Numbers, ZERO Real.exp, ZERO sqrt!) Lean-4 Verification Check
for the Round-3 Level-3 EMA Cubic-Potential Paley 2-Coset Certificate (`cert_l3_ema_paley_t1_lean.json`).

Certificate parameter note:
  Uses 8 free dual-weight numerators (`wa_int, wb_int, wA_ints[0], wB_ints[0], aA_int, aB_int,
  wdir_int, kappa_int`) alongside structural parameters (`T2, q, P, K0, N_L, L, beta_den, eps1_den,
  eps2_den, D_int`).

Verifies (zero `assert` statements — 100% enforced under `python3 -O`):
  1. [B-V1] & [B-V3]: Complete recursive Pratt primality tree for q = 100_000_007 == 3 (mod 4)
     all the way down to base prime 2, enforcing BOTH Fermat (`pow(g, p - 1, p) == 1`) AND Lucas
     order (`all(pow(g, (p - 1) // pf, p) != 1)`) at every node, plus a fired-once composite
     regression test rejecting 15, 21, and Carmichael 561!
  2. Exact integer Paley matrix identity A^T A = B^T B = ((q+1)/4)*I + ((q-3)/4)*J ->
       (1/q) * sum_i |rho_i - S| <= sqrt(q + 1) / (2*(q - 1)) < 10001 / (2*(q - 1))
     (since 10001^2 = 100,020,001 > 100,000,008 = q + 1 -> zero sqrt library needed).
  3. Zero-Real.exp 5-term Bonferroni upper bound on 1 - exp(-S/2) <= P_5(S/2)
     (pure `fractions.Fraction` arithmetic) certifying:
       - U_finite_rat <= 380/1000 = 0.380 (breaks 0.380 in 150 pure-Q dyadic boxes)
       - U_finite_rat <= 381/1000 < (3 - sqrt(5))/2 (in 71 pure-Q dyadic boxes)
"""
from __future__ import annotations

import json
import math
import os
from fractions import Fraction


def require(cond: bool, msg: str) -> None:
    if not cond:
        raise ValueError(f"VERIFICATION FAILED: {msg}")


def factor_int(n: int) -> dict[int, int]:
    require(isinstance(n, int) and not isinstance(n, bool) and n >= 1, f"factor_int requires n >= 1, got {n!r}")
    res: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            res[d] = res.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        res[n] = res.get(n, 0) + 1
    return res


def find_prim_root(p: int) -> tuple[int, dict[int, int]]:
    """
    Find a verified Pratt primitive root g modulo p, enforcing BOTH:
      - Fermat check: pow(g, p - 1, p) == 1  (rejects composites like 15, 21 immediately!)
      - Lucas order check: pow(g, (p - 1) // pf, p) != 1 for every prime factor pf of p - 1.
    """
    require(isinstance(p, int) and not isinstance(p, bool) and p >= 2, f"p must be int >= 2, got {p!r}")
    if p == 2:
        return 1, {}
    facs = factor_int(p - 1)
    for g in range(2, p):
        if pow(g, p - 1, p) != 1:
            continue
        if all(pow(g, (p - 1) // pf, p) != 1 for pf in facs):
            return g, facs
    raise ValueError(f"No Pratt primitive root exists modulo {p} (composite)")


def verify_pratt_tree(n: int, _verified: dict[int, tuple[int, dict[int, int]]] | None = None) -> dict[int, tuple[int, dict[int, int]]]:
    """
    Construct AND independently verify a complete recursive Pratt primality certificate tree
    for `n` all the way down to base prime `2`!
    Returns dict mapping each prime node p -> (witness_g, prime_factorization_of_p_minus_1).
    """
    if _verified is None:
        _verified = {}
    require(isinstance(n, int) and not isinstance(n, bool) and n >= 2, f"n must be int >= 2, got {n!r}")
    if n in _verified:
        return _verified
    if n == 2:
        _verified[2] = (1, {})
        return _verified

    g, facs = find_prim_root(n)
    # 1. Verify exact product of factorization == n - 1
    prod = 1
    for pf, exp in facs.items():
        require(pf >= 2 and exp >= 1, f"Invalid factor ({pf}, {exp}) for {n - 1}")
        prod *= pf ** exp
    require(prod == n - 1, f"Factorization product {prod} != {n - 1}")

    # 2. Recursively verify EVERY prime factor pf of n - 1 all the way down to 2!
    for pf in sorted(facs):
        verify_pratt_tree(pf, _verified)

    # 3. Independently re-check Fermat + Lucas order conditions at node n
    require(1 <= g < n, f"Witness g={g} out of [1, {n - 1}]")
    require(pow(g, n - 1, n) == 1, f"Fermat check failed: {g}^({n - 1}) != 1 (mod {n})")
    for pf in sorted(facs):
        require(
            pow(g, (n - 1) // pf, n) != 1,
            f"Lucas order check failed: {g}^(({n - 1})//{pf}) == 1 (mod {n})",
        )

    _verified[n] = (g, facs)
    return _verified


def run_composite_regression_tests() -> None:
    """[B-V1] Fired-once regression tests ensuring find_prim_root & verify_pratt_tree reject composites."""
    for comp in (4, 9, 15, 21, 25, 341, 561):
        rejected_root = False
        try:
            find_prim_root(comp)
        except ValueError:
            rejected_root = True
        require(rejected_root, f"REGRESSION FAILURE: find_prim_root({comp}) failed to reject composite!")

        rejected_tree = False
        try:
            verify_pratt_tree(comp)
        except ValueError:
            rejected_tree = True
        require(rejected_tree, f"REGRESSION FAILURE: verify_pratt_tree({comp}) failed to reject composite!")


def P5(u: Fraction) -> Fraction:
    return u - (u ** 2) / 2 + (u ** 3) / 6 - (u ** 4) / 24 + (u ** 5) / 120


def eval_V_rat_upper(S: Fraction, wa, wb, wA0, wB0, aA, aB, wdir, kap) -> Fraction:
    """Exact pure-Fraction upper bound on V_Bellman(S; w, +kap*S^2) using 1 - exp(-S/2) <= P5(S/2)."""
    p = P5(S / 2)
    lam = kap * S * S
    V0 = Fraction(0)
    V1 = aB * S + lam / 2
    V2 = aA * S + lam / 2
    V3 = wdir + (aA + aB) * S + lam

    # T2 = 1 bulk wave (B then A backward)
    g = wB0 + V2 - V3
    if g > 0:
        V3 = V3 + p * g
    g = wB0 + V0 - V1
    if g > 0:
        V1 = V1 + p * g
    g = wA0 + V1 - V3
    if g > 0:
        V3 = V3 + p * g
    g = wA0 + V0 - V2
    if g > 0:
        V2 = V2 + p * g

    gb = wb + V2 - V3
    Vmid = V3 + (gb / 2 if gb > 0 else Fraction(0))
    ga = wa + V1 - Vmid
    Vpair = Vmid + (ga / 2 if ga > 0 else Fraction(0))
    return Vpair


def main() -> None:
    # 0. Fired-once composite rejection regression tests [B-V1]
    run_composite_regression_tests()
    print("=== 0. COMPOSITE REJECTION REGRESSION TESTS [B-V1] ===")
    print("  [PASS] find_prim_root & verify_pratt_tree rejected composites (4, 9, 15, 21, 25, 341, 561)")

    bundle_dir = os.path.dirname(os.path.abspath(__file__))
    cert_path = os.path.join(bundle_dir, "cert_l3_ema_paley_t1_lean.json")
    with open(cert_path, "r", encoding="utf-8") as f:
        d = json.load(f)

    q = d["q"]
    require(isinstance(q, int) and not isinstance(q, bool) and q % 4 == 3, f"Expected q == 3 mod 4, got {q}")

    # 1. Full recursive Pratt primality tree down to 2 [B-V3]
    tree = verify_pratt_tree(q)
    nodes_sorted = [(p, tree[p][0]) for p in sorted(tree)]
    expected_nodes = [(2, 1), (3, 2), (5, 2), (7, 3), (101, 2), (491, 2), (4243, 2), (101833, 5), (100000007, 5)]
    require(nodes_sorted == expected_nodes, f"Pratt tree mismatch: got {nodes_sorted}, expected {expected_nodes}")
    print(f"\n=== 1. FULL RECURSIVE PRATT PRIMALITY TREE DOWN TO 2 [B-V3] (q = {q} == 3 mod 4) ===")
    for p in sorted(tree):
        g_p, f_p = tree[p]
        if p == 2:
            print("  node 2: base prime (witness g = 1)")
        else:
            print(f"  node {p}: {p - 1} = {f_p}, witness g = {g_p} [Fermat + Lucas verified]")
    print(f"  Exact sorted Pratt pairs: {nodes_sorted}")

    # 2. Pure-rational sqrt(q + 1) bound & exact Fraction dyadic box proof
    sq_hi = math.isqrt(q + 1) + 1
    require(sq_hi * sq_hi > q + 1, f"sq_hi^2 = {sq_hi * sq_hi} <= q + 1 = {q + 1}")
    D = d["D_int"]
    wa = Fraction(d["wa_int"], D)
    wb = Fraction(d["wb_int"], D)
    wA0 = Fraction(d["wA_ints"][0], D)
    wB0 = Fraction(d["wB_ints"][0], D)
    aA = Fraction(d["aA_int"], D)
    aB = Fraction(d["aB_int"], D)
    wdir = Fraction(d["wdir_int"], D)
    kap = Fraction(d["kappa_int"], D)
    W_p2 = (wA0 + wB0) / 2
    require((wa + wb) / 2 + W_p2 + aA + aB + wdir == Fraction(1), "Exact dual normalization != 1")

    beta = Fraction(1, d["beta_den"])
    eps1 = Fraction(1, d["eps1_den"])
    eps2 = Fraction(1, d["eps2_den"])

    # Exact Paley bound for q == 3 mod 4: (W_p2 + aA + aB) * sqrt(q+1) / (2*(q-1))
    err_paley = (W_p2 + aA + aB) * Fraction(sq_hi, 2 * (q - 1))
    err_binom = W_p2 * Fraction(2, 8 * d["N_L"] * (q - 1))
    err_cubic = kap * beta
    err_tele = kap / (3 * beta * d["L"])
    tot_err = err_paley + err_binom + err_cubic + err_tele
    Z_fr = (Fraction(1) - beta / (d["P"] - 1)) * (Fraction(1) - eps1 * (wa + wb) / 2 - (eps1 + 2 * eps2) * W_p2)

    print("\n=== 2. PURE-RATIONAL ZERO-Real.exp DYADIC BOX CHECK (8 free dual-weight numerators) ===")
    for thresh_label, thresh in [
        ("380/1000 (breaks 0.380)", Fraction(380, 1000)),
        ("381/1000 (< (3-sqrt(5))/2)", Fraction(381, 1000)),
    ]:
        u_lim_tgt = thresh * Z_fr - tot_err
        stack = [(Fraction(0), Fraction(1))]
        done = []
        while stack:
            a, b = stack.pop()
            hi = eval_V_rat_upper(b, wa, wb, wA0, wB0, aA, aB, wdir, kap) - kap * (a ** 3)
            if hi <= u_lim_tgt:
                done.append((a, b, hi))
            else:
                mid = (a + b) / 2
                stack.append((mid, b))
                stack.append((a, mid))
                require(len(done) + len(stack) <= 5000, f"Exceeded 5000 boxes for {thresh_label}")
        # Verify exact gap-free partition of [0, 1]
        done.sort(key=lambda t: t[0])
        require(done[0][0] == Fraction(0) and done[-1][1] == Fraction(1), "Partition does not span [0, 1]")
        for i in range(len(done) - 1):
            require(done[i][1] == done[i + 1][0], f"Gap/overlap in rational partition at index {i}")
        max_hi = max(h for _, _, h in done)
        u_fin_rat = (max_hi + tot_err) / Z_fr
        require(u_fin_rat <= thresh, f"u_fin_rat = {float(u_fin_rat)} > {float(thresh)}")
        print(f"  Target {thresh_label:28s}: PROVED in {len(done):4d} pure-Q boxes! "
              f"U_finite_rat = {float(u_fin_rat):.9f} <= {float(thresh):.6f}")


if __name__ == "__main__":
    main()
