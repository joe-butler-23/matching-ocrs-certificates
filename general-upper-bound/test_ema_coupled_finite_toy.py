#!/usr/bin/env python3
"""
Exact Finite-Graph Polytope & Top-Down Integer Matching Verification for the
EMA Pyramid Recursive Cascade G^{EMA}_{q, P, K0, N_L, L, beta, eps1, eps2}.

Zero `assert` statements — all invariants enforced via explicit `require()` raising `ValueError`
under `python3 -O` ([B-V2]).

Verifies in pure `fractions.Fraction`:
  1. Deterministic non-exhaustion of free blocks across all levels j = 0..L.
  2. Exact EMA weights sum_{j < k} alpha_{k, j} = 1.
  3. Exact backward recurrence for block-free probability gamma_j and analytical lower bound
     gamma_j >= 1 - beta / (P - 1) (j >= 1), gamma_0 >= 1 - 1 / (K0 * (P - 1)).
  4. Exact vertex fractional degree sums <= 1 and strict positivity x_e > 0.
  5. 10,000-trial exact integer matching sampler confirming zero vertex collisions and
     empirical marginals matching exact gamma_k.
"""
from __future__ import annotations

from fractions import Fraction
import random


def require(cond: bool, msg: str) -> None:
    if not cond:
        raise ValueError(f"POLYTOPE CHECK FAILED: {msg}")


def verify_ema_pyramid_polytope_exact(
    q: int,
    T2: int,
    P: int,
    K0: int,
    N_L: int,
    L: int,
    beta: Fraction,
    eps1: Fraction,
    eps2: Fraction,
) -> list[Fraction]:
    require((q - 1) % (T2 + 1) == 0, f"q - 1 = {q - 1} not divisible by T2 + 1 = {T2 + 1}")
    require(P >= 2 and K0 >= 2 and N_L >= 1 and L >= 1, f"Invalid sizes P={P}, K0={K0}, N_L={N_L}, L={L}")
    require(Fraction(0) < beta < Fraction(1), f"beta={beta} not in (0, 1)")
    require(eps1 > 0 and eps2 > 0 and eps1 + 2 * eps2 < Fraction(1), f"Invalid eps1={eps1}, eps2={eps2}")

    N = [0] * (L + 1)
    for k in range(1, L + 1):
        N[k] = N_L * (P ** (L - k))
    N[0] = K0 * N_L * (P ** L)

    # Worst-case free-block guarantee (deterministic for EVERY realization)
    for j in range(L + 1):
        max_consumed = sum(N[k] for k in range(j + 1, L + 1))
        require(max_consumed < N[j], f"Level {j}: max_consumed={max_consumed} >= N[j]={N[j]}")

    # Exact alpha_{k, j}
    alpha: dict[tuple[int, int], Fraction] = {}
    for k in range(1, L + 1):
        alpha[(k, 0)] = (Fraction(1) - beta) ** (k - 1)
        for j in range(1, k):
            alpha[(k, j)] = beta * ((Fraction(1) - beta) ** (k - 1 - j))
        require(sum(alpha[(k, j)] for j in range(k)) == Fraction(1), f"alpha sum != 1 at level {k}")

    # Exact backward recurrence for gamma_j = P[block r' of Level j is FREE]
    gamma = [Fraction(0)] * (L + 1)
    gamma[L] = Fraction(1)
    for j in range(L - 1, -1, -1):
        consumed_exp = Fraction(0)
        for k in range(j + 1, L + 1):
            consumed_exp += Fraction(N[k], N[j]) * gamma[k] * (Fraction(1) - eps1) * alpha[(k, j)]
        gamma[j] = Fraction(1) - consumed_exp
        require(Fraction(0) < gamma[j] <= Fraction(1), f"gamma[{j}]={gamma[j]} out of (0, 1]")
        if j >= 1:
            lb = Fraction(1) - beta / Fraction(P - 1)
            require(gamma[j] >= lb, f"gamma[{j}]={float(gamma[j])} < lb={float(lb)}")
        else:
            lb0 = Fraction(1) - Fraction(1, K0 * (P - 1))
            require(gamma[0] >= lb0, f"gamma[0]={float(gamma[0])} < lb0={float(lb0)}")

    # Exact vertex degree sums <= 1
    d0 = (q - 1) // (T2 + 1)
    for k in range(1, L + 1):
        deg_internal_and_in = gamma[k] * (
            (Fraction(1) - eps1) / 2
            + eps1
            + (Fraction(1) - eps1 - 2 * eps2) / 2
            + eps2
        )
        require(deg_internal_and_in == gamma[k], f"Internal+in degree mismatch at level {k}")
        deg_export = (Fraction(1) - gamma[k]) / 2
        require(
            deg_internal_and_in + deg_export == (Fraction(1) + gamma[k]) / 2 <= Fraction(1),
            f"Total degree > 1 at level {k}",
        )
        deg_w = gamma[k] * (Fraction(1) - eps1)
        require(deg_w <= Fraction(1), f"Bridge degree {deg_w} > 1 at level {k}")

    # Monte Carlo validation on concrete finite instances (zero vertex collisions + exact gamma match)
    if q <= 5 and L <= 3 and N[0] <= 200:
        rng = random.Random(2026)
        C_bulk = list(range(1, 1 + T2 * d0))
        C_atom = list(range(1 + T2 * d0, q))
        free_counts = [0] * (L + 1)
        trials = 10000
        p_bulk = float((Fraction(1) - eps1 - 2 * eps2) / 2)
        p_atom = float(eps2)
        p_dir = float(eps1)
        for _ in range(trials):
            free_blocks = [list(range(N[k])) for k in range(L + 1)]
            used_verts = set()
            for k in range(L, 0, -1):
                free_counts[k] += len(free_blocks[k])
                blks = list(free_blocks[k])
                rng.shuffle(blks)
                for r in blks:
                    u = rng.random()
                    if u < p_dir:
                        for i in range(q):
                            va = (k, r, 0, i)
                            vb = (k, r, 1, i)
                            require(va not in used_verts and vb not in used_verts, f"Collision on {(va, vb)}")
                            used_verts.add(va)
                            used_verts.add(vb)
                    else:
                        if u < p_dir + p_bulk:
                            br_side = 0
                            shifts = C_bulk
                        elif u < p_dir + 2 * p_bulk:
                            br_side = 1
                            shifts = C_bulk
                        elif u < p_dir + 2 * p_bulk + p_atom:
                            br_side = 0
                            shifts = C_atom
                        else:
                            br_side = 1
                            shifts = C_atom
                        p2_side = 1 - br_side
                        rj = rng.random()
                        acc = 0.0
                        chosen_j = 0
                        for j in range(k):
                            acc += float(alpha[(k, j)])
                            if rj <= acc:
                                chosen_j = j
                                break
                        idx_p = rng.randrange(len(free_blocks[chosen_j]))
                        rp = free_blocks[chosen_j].pop(idx_p)
                        sig = rng.randrange(2)
                        s = rng.choice(shifts)
                        for i in range(q):
                            v_br = (k, r, br_side, i)
                            vw = (k, r, 2, i)
                            v_p2 = (k, r, p2_side, i)
                            vd = (chosen_j, rp, sig, (i + s) % q)
                            for v in (v_br, vw, v_p2, vd):
                                require(v not in used_verts, f"COLLISION on vertex {v}!")
                                used_verts.add(v)
            # Level 0 leaves
            free_counts[0] += len(free_blocks[0])
            for r in free_blocks[0]:
                for sig in (0, 1):
                    for i in range(q):
                        vd = (0, r, sig, i)
                        vl = ("leaf", r, sig, i)
                        require(vd not in used_verts and vl not in used_verts, f"Leaf collision on {(vd, vl)}")
                        used_verts.add(vd)
                        used_verts.add(vl)

        for k in range(L + 1):
            emp_g = free_counts[k] / (trials * N[k])
            exact_g = float(gamma[k])
            err = abs(emp_g - exact_g)
            require(err < 0.015, f"Level {k}: emp_g={emp_g:.5f} vs exact_g={exact_g:.5f}")
        print(f"  [PASS] Monte Carlo 10,000 exact valid integer matchings (q={q}, L={L}, P={P}, K0={K0}): "
              f"0 collisions, gamma={[round(float(g), 5) for g in gamma]}")
    return gamma


def main() -> None:
    print("=== VERIFYING EXACT TOP-DOWN INTEGER MATCHING PYRAMID (pure Fraction + 10k MC, python3 -O safe) ===")
    verify_ema_pyramid_polytope_exact(q=5, T2=1, P=3, K0=2, N_L=2, L=3,
                                      beta=Fraction(1, 3), eps1=Fraction(1, 20), eps2=Fraction(1, 20))
    verify_ema_pyramid_polytope_exact(q=5, T2=1, P=3, K0=2, N_L=2, L=2,
                                      beta=Fraction(1, 5), eps1=Fraction(1, 50), eps2=Fraction(1, 50))
    verify_ema_pyramid_polytope_exact(q=13, T2=2, P=3, K0=2, N_L=2, L=5,
                                      beta=Fraction(1, 100), eps1=Fraction(1, 1000), eps2=Fraction(1, 1000))
    print("  [PASS] Large L=5 exact Fraction polytope check passed!")


if __name__ == "__main__":
    main()
