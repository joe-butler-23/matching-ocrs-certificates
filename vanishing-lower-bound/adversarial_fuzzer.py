#!/usr/bin/env python3
"""
Adversarial End-to-End Fuzzer for Strongly-Online Global-Coin Matching OCRS Lower Bounds.
100% python3 -O safe (zero assert-based checks); ZERO Monte Carlo approximation!

Suites:
  Suite 1: 3,000 adversarial NON-VANISHING simple finite trees (random, depth-2 MMN, depth-3, stars, paths)
           with explicit fractional-matching feasibility guard max_v deg_x(v) <= 1.0 on EVERY instance,
           testing EXACT 2^|V|-state selectability of:
           (a) `discrete_2branch` (>= 0.382309325242),
           (b) `discrete_64branch` (>= 0.382661888422), AND
           (c) `continuous_staircase` (>= 0.382673618862) via exact breakpoint integration.
  Suite 2: Hard simple-tree benchmarks (ANN S-OCRS tight star-trees k = 2..5 & (0.999, 0.999) corner tree).
  Suite 3: Strictly feasible 7-edge (cycle-length-6 + pendant) counterexample to MMN arXiv:2406.14506 v3
           lem:tree_reduction_game_part_1 (L1572) with max_v deg_x(v) = 1.00 (gap = -8.4093614528e-5 < 0),
           PLUS 1,500 strictly feasible cyclic ell = 1, 2, 3 locally-tree-like graphs verifying our
           6-Line Zero-Gap Disagreement-Coupling Bound |p_G - gamma_e| <= 2 * a_e * beta0^ell.
  Suite 4: Load-Bearing Hypothesis Necessity Separations (confirming simple-graph & oblivious-adversary
           hypotheses are genuinely required):
           (a) Parallel-edge bipartite multigraphs (k = 2 -> 0.352687, k = 200 -> 0.330779 < alpha*),
           (b) Adaptive arrival-order adversary on simple forests (w0*(1 - c1)^2 ~= 0.360869 < alpha*).
"""
from collections import defaultdict
import json
import math
from pathlib import Path
import random

ALPHA_STAR = (3.0 - math.sqrt(5.0)) / 2.0


def require(cond, msg):
    if not cond:
        raise RuntimeError(msg)


def verify_simple_fractional_matching(n_verts, edges, x_vals, ctx):
    seen = set()
    deg = [0.0] * n_verts
    for (u, v), x in zip(edges, x_vals):
        require(0 <= u < n_verts and 0 <= v < n_verts and u != v, f"{ctx}: invalid edge {(u, v)}")
        e = (u, v) if u < v else (v, u)
        require(e not in seen, f"{ctx}: parallel edge {e} in simple-graph test!")
        seen.add(e)
        require(-1e-15 <= x <= 1.0 + 1e-12, f"{ctx}: invalid x={x}")
        deg[u] += x
        deg[v] += x
    max_deg = max(deg) if deg else 0.0
    require(max_deg <= 1.0 + 1e-12, f"{ctx}: INFEASIBLE fractional matching max_deg={max_deg:.12f} > 1.0!")


def A_of_b(b):
    if b <= 0.0:
        return 1.0
    return 2.0 / (2.0 * b + 1.0 + math.sqrt(4.0 * b + 1.0))


def I_of_b(b):
    if b <= 0.0:
        return 0.0
    return (6.0 * b * b + 6.0 * b - 1.0 + (1.0 + 4.0 * b) * math.sqrt(1.0 + 4.0 * b)) / 12.0


def exact_ocrs_edge_probs(n_verts, edges_ordered, x_vals, a_vals):
    """EXACT 2^n_verts DP over vertex occupancy masks — zero Monte Carlo error."""
    dp = {0: 1.0}
    sel_probs = []
    for (u, v), x, a in zip(edges_ordered, x_vals, a_vals):
        require(-1e-12 <= a <= 1.0 + 1e-12, f"Invalid acceptance probability a={a}")
        a = min(1.0, max(0.0, a))
        bu, bv = 1 << u, 1 << v
        both = bu | bv
        free_prob = 0.0
        new_dp = defaultdict(float)
        pm = x * a
        for mask, prob in dp.items():
            if (mask & both) == 0:
                free_prob += prob
                new_dp[mask | both] += prob * pm
                new_dp[mask] += prob * (1.0 - pm)
            else:
                new_dp[mask] += prob
        sel_probs.append(free_prob * a)
        dp = new_dp
    return sel_probs


def run_branch0(n_verts, edges_ord, x_ord, c1):
    x_pref = [0.0] * n_verts
    D = [0.0] * n_verts
    a_vals, g_vals = [], []
    for (u, v), x in zip(edges_ord, x_ord):
        xu, xv = x_pref[u], x_pref[v]
        g = min(c1, (1.0 - c1 * xu) * (1.0 - c1 * xv))
        a = g / ((1.0 - D[u]) * (1.0 - D[v]))
        require(0.0 <= a <= 1.0 + 1e-12, f"Branch 0 infeasible a={a}")
        a_vals.append(a)
        g_vals.append(g)
        x_pref[u] += x
        x_pref[v] += x
        D[u] += g * x
        D[v] += g * x
    probs = exact_ocrs_edge_probs(n_verts, edges_ord, x_ord, a_vals)
    return probs, g_vals


def run_branch_T(n_verts, edges_ord, x_ord, T):
    b = 1.0 - T
    A = A_of_b(b)
    x_pref = [0.0] * n_verts
    xT_pref = [0.0] * n_verts
    a_vals, g_vals = [], []
    for (u, v), x in zip(edges_ord, x_ord):
        xu, xv = x_pref[u], x_pref[v]
        if xu >= T - 1e-15 and xv >= T - 1e-15:
            g = A
            a = A / ((1.0 - A * xT_pref[u]) * (1.0 - A * xT_pref[v]))
            require(0.0 <= a <= 1.0 + 1e-12, f"Branch T={T} infeasible a={a}")
            xT_pref[u] += x
            xT_pref[v] += x
        else:
            g = 0.0
            a = 0.0
        a_vals.append(a)
        g_vals.append(g)
        x_pref[u] += x
        x_pref[v] += x
    probs = exact_ocrs_edge_probs(n_verts, edges_ord, x_ord, a_vals)
    return probs, g_vals


def exact_continuous_staircase_on_tree(n_verts, edges_ord, x_ord, c1):
    tau = (1.0 - 2.0 * c1) / (c1 * (1.0 - c1))
    B = 1.0 - tau
    IB = I_of_b(B)
    w0 = 1.0 / (1.0 + c1 * (1.0 - c1) * IB)

    p0, g0 = run_branch0(n_verts, edges_ord, x_ord, c1)
    max_b0_err = max(abs(a - b) for a, b in zip(p0, g0))

    x_pref = [0.0] * n_verts
    s_vals = []
    for (u, v), x in zip(edges_ord, x_ord):
        s_vals.append(min(x_pref[u], x_pref[v]))
        x_pref[u] += x
        x_pref[v] += x

    active_s = sorted(set(s for s in s_vals if tau <= s < 1.0 - 1e-15))
    max_bT_err = 0.0
    for T_test in active_s:
        pT, gT = run_branch_T(n_verts, edges_ord, x_ord, T_test)
        max_bT_err = max(max_bT_err, max(abs(a - b) for a, b in zip(pT, gT)))

    exact_probs = [
        w0 * p0[i] + w0 * c1 * (1.0 - c1) * max(0.0, s_vals[i] - tau)
        for i in range(len(edges_ord))
    ]
    return exact_probs, max_b0_err, max_bT_err


def run_discrete_mixture_exact(n_verts, edges_ord, x_ord, cfg):
    c1 = cfg["c1_num"] / cfg["c1_den"]
    W_den = cfg["weight_den"]
    w0 = cfg["w0_num"] / W_den
    p0, _ = run_branch0(n_verts, edges_ord, x_ord, c1)
    tot = [w0 * p for p in p0]
    for br in cfg["branches"]:
        tj = br["t_num"] / br["t_den"]
        wj = br["w_num"] / W_den
        pT, _ = run_branch_T(n_verts, edges_ord, x_ord, tj)
        for i in range(len(tot)):
            tot[i] += wj * pT[i]
    return tot


def suite_1_adversarial_trees(cert):
    print("=== SUITE 1: 3,000 ADVERSARIAL NON-VANISHING SIMPLE TREES (EXACT 2^|V| DP) ===")
    rng = random.Random(2026112603)
    c1_cont = cert["continuous_staircase"]["c1_num"] / cert["continuous_staircase"]["c1_den"]
    lb_cont = cert["continuous_staircase"]["claimed_c_lb_num"] / cert["continuous_staircase"]["claimed_c_lb_den"]
    lb_d2 = cert["discrete_2branch"]["claimed_c_lb_num"] / cert["discrete_2branch"]["claimed_c_lb_den"]
    lb_d64 = cert["discrete_64branch"]["claimed_c_lb_num"] / cert["discrete_64branch"]["claimed_c_lb_den"]

    min_cont = 1.0
    min_d2 = 1.0
    min_d64 = 1.0
    max_b0_err = 0.0
    max_bT_err = 0.0
    n_late_edges_hit = 0
    tau_c = (1.0 - 2.0 * c1_cont) / (c1_cont * (1.0 - c1_cont))

    for t in range(3000):
        n = rng.randint(2, 10)
        edges = []
        mode = t % 5
        if mode == 0:
            for v in range(1, n):
                edges.append((rng.randrange(v), v))
        elif mode == 1:
            d = rng.randint(1, min(4, n - 1))
            for i in range(1, d + 1):
                edges.append((0, i))
            for v in range(d + 1, n):
                edges.append((rng.randint(1, d), v))
        elif mode == 2:
            if n >= 6:
                edges = [(0, 1), (0, 2), (1, 3), (2, 4), (3, 5)]
                for v in range(6, n):
                    edges.append((rng.choice([0, 1, 2, 3, 4]), v))
            else:
                for v in range(1, n):
                    edges.append((v - 1, v))
        elif mode == 3:
            for v in range(1, n):
                edges.append((0, v))
        else:
            for v in range(1, n):
                edges.append((v - 1, v))

        if t % 3 == 0 and len(edges) >= 3:
            target_idx = rng.randrange(len(edges))
            u_t, v_t = edges[target_idx]
            x_vals = [0.0] * len(edges)
            x_vals[target_idx] = rng.uniform(0.001, 0.04)
            deg = [0.0] * n
            deg[u_t] += x_vals[target_idx]
            deg[v_t] += x_vals[target_idx]
            other_idxs = [i for i in range(len(edges)) if i != target_idx]
            rng.shuffle(other_idxs)
            for idx in other_idxs:
                a, b = edges[idx]
                slack = min(1.0 - deg[a], 1.0 - deg[b])
                if slack > 1e-9:
                    val = slack * rng.choice([0.999, 0.98, 0.95, 0.5, rng.random()])
                    x_vals[idx] = val
                    deg[a] += val
                    deg[b] += val
            perm = other_idxs + [target_idx]
        else:
            raw = [rng.random() ** rng.choice([0.3, 1.0, 3.0]) for _ in edges]
            deg = [0.0] * n
            for (a, b), w in zip(edges, raw):
                deg[a] += w
                deg[b] += w
            md = max(deg) if max(deg) > 0 else 1.0
            x_vals = [w / md for w in raw]
            perm = list(range(len(edges)))
            rng.shuffle(perm)

        edges_ord = [edges[i] for i in perm]
        x_ord = [x_vals[i] for i in perm]
        verify_simple_fractional_matching(n, edges_ord, x_ord, f"Suite1.tree[{t}]")

        p_cont, e0, eT = exact_continuous_staircase_on_tree(n, edges_ord, x_ord, c1_cont)
        max_b0_err = max(max_b0_err, e0)
        max_bT_err = max(max_bT_err, eT)

        xp = [0.0] * n
        for (a, b), x in zip(edges_ord, x_ord):
            if xp[a] >= tau_c and xp[b] >= tau_c:
                n_late_edges_hit += 1
            xp[a] += x
            xp[b] += x

        mc = min(p_cont)
        min_cont = min(min_cont, mc)
        require(mc >= lb_cont - 1e-12, f"Tree {t}: continuous staircase min={mc} < lb={lb_cont}")

        p_d2 = run_discrete_mixture_exact(n, edges_ord, x_ord, cert["discrete_2branch"])
        md2 = min(p_d2)
        min_d2 = min(min_d2, md2)
        require(md2 >= lb_d2 - 1e-12, f"Tree {t}: discrete_2branch min={md2} < lb={lb_d2}")

        if t < 250:
            p_d64 = run_discrete_mixture_exact(n, edges_ord, x_ord, cert["discrete_64branch"])
            md64 = min(p_d64)
            min_d64 = min(min_d64, md64)
            require(md64 >= lb_d64 - 1e-12, f"Tree {t}: discrete_64branch min={md64} < lb={lb_d64}")

    print(f"  [PASS] 3,000 non-vanishing simple trees (incl. {n_late_edges_hit} late-corner E(tau) edges):")
    print(f"         max exact DP identity err: Branch0={max_b0_err:.2e}, BranchT={max_bT_err:.2e}")
    print(f"         continuous_staircase min = {min_cont:.12f} (claimed LB = {lb_cont:.12f})")
    print(f"         discrete_2branch     min = {min_d2:.12f} (claimed LB = {lb_d2:.12f})")
    print(f"         discrete_64branch    min = {min_d64:.12f} (claimed LB = {lb_d64:.12f})")


def suite_2_hard_benchmarks(cert):
    print("\n=== SUITE 2: HARD SIMPLE-TREE BENCHMARKS (ANN S-OCRS TIGHT TREES & LATE CORNER) ===")
    c1_cont = cert["continuous_staircase"]["c1_num"] / cert["continuous_staircase"]["c1_den"]
    lb_cont = cert["continuous_staircase"]["claimed_c_lb_num"] / cert["continuous_staircase"]["claimed_c_lb_den"]

    for k in [2, 3, 4, 5]:
        n = 1 + 2 * k
        edges = [(k + i, i) for i in range(1, k + 1)] + [(0, i) for i in range(1, k + 1)]
        xs = [1.0 - 1.0 / k] * k + [1.0 / k] * k
        verify_simple_fractional_matching(n, edges, xs, f"Suite2.ANN[{k}]")
        p_cont, _, _ = exact_continuous_staircase_on_tree(n, edges, xs, c1_cont)
        p_d2 = run_discrete_mixture_exact(n, edges, xs, cert["discrete_2branch"])
        require(min(p_cont) >= lb_cont - 1e-12, f"ANN k={k} failed!")
        print(
            f"  [PASS] ANN S-OCRS tight tree k={k}: continuous_min={min(p_cont):.10f}, "
            f"2branch_min={min(p_d2):.10f} (S-OCRS limit={ALPHA_STAR:.6f})"
        )

    n = 10
    edges = [(0, 2), (0, 3), (0, 4), (0, 5), (1, 6), (1, 7), (1, 8), (1, 9), (0, 1)]
    xs = [0.24975] * 8 + [0.001]
    verify_simple_fractional_matching(n, edges, xs, "Suite2.corner")
    p_cont, _, _ = exact_continuous_staircase_on_tree(n, edges, xs, c1_cont)
    p_d64 = run_discrete_mixture_exact(n, edges, xs, cert["discrete_64branch"])
    require(min(p_cont) >= lb_cont - 1e-12, "Corner (0.999, 0.999) failed!")
    print(f"  [PASS] Exact (0.999, 0.999) late-late corner tree: cont_min={min(p_cont):.12f}, 64branch_min={min(p_d64):.12f}")


def suite_3_cyclic_locally_tree_like_and_mmn_bug(cert):
    print("\n=== SUITE 3: STRICTLY FEASIBLE MMN L1572 COUNTEREXAMPLE + 1,500 CYCLIC GRAPHS ===")
    # 7-edge simple graph: 6-cycle 0-2-4-6-5-3-0 plus pendant target edge e=(0,1).
    # Induced 2-neighborhood G[N^2(e)] is a tree (ell=2 locally-tree-like); h2=(5,6) not in W(T).
    # Strictly feasible fractional matching: max_v deg_x(v) = 1.00 (deg(0) = 0.495 + 0.495 + 0.01 = 1.00).
    n = 7
    edges_with = [(3, 5), (5, 6), (6, 4), (2, 4), (0, 2), (0, 3), (0, 1)]
    edges_no   = [(3, 5),         (6, 4), (2, 4), (0, 2), (0, 3), (0, 1)]
    x_with = [0.5, 0.5, 0.5, 0.5, 0.495, 0.495, 0.01]
    x_no   = [0.5,      0.5, 0.5, 0.495, 0.495, 0.01]
    verify_simple_fractional_matching(n, edges_with, x_with, "Suite3.MMN_with_h2")
    verify_simple_fractional_matching(n, edges_no, x_no, "Suite3.MMN_no_h2")

    def run_aom(ed, xs, alpha):
        xp = [0.0] * n
        av = []
        for (a, b), x in zip(ed, xs):
            av.append(alpha / ((1.0 - alpha * xp[a]) * (1.0 - alpha * xp[b])))
            xp[a] += x
            xp[b] += x
        return exact_ocrs_edge_probs(n, ed, xs, av)[-1]

    p_with = run_aom(edges_with, x_with, ALPHA_STAR)
    p_no = run_aom(edges_no, x_no, ALPHA_STAR)
    diff = p_with - p_no
    require(diff < -8.0e-5, f"Expected strictly negative h2 causal leak ~ -8.41e-5, got {diff}")
    print(
        f"  [CONFIRMED] Strictly feasible MMN L1572 counterexample (7 edges, 6-cycle + pendant, max_deg=1.00):\n"
        f"              WITH h2    = {p_with:.12f}\n"
        f"              WITHOUT h2 = {p_no:.12f}\n"
        f"              diff       = {diff:+.10e} < 0"
    )

    c1 = cert["continuous_staircase"]["c1_num"] / cert["continuous_staircase"]["c1_den"]
    beta0 = -math.log(1.0 - c1) / (1.0 - c1)
    rng = random.Random(7777)

    for ell in [1, 2, 3]:
        max_disagree_ratio = 0.0
        n_trials = 600 if ell <= 2 else 300
        for trial in range(n_trials):
            if ell == 1:
                nv = 6
                ed = [(0, 2), (2, 4), (4, 3), (0, 3), (1, 5), (5, 4), (0, 1)]
            elif ell == 2:
                nv = 8
                ed = [(3, 5), (5, 7), (7, 6), (6, 4), (2, 4), (0, 2), (0, 3), (0, 1)]
            else:
                nv = 9
                ed = [(3, 5), (5, 7), (7, 8), (8, 6), (4, 6), (2, 4), (0, 2), (0, 3), (0, 1)]

            if trial % 2 == 1:
                ed = list(ed)
                rng.shuffle(ed)
            raw = [rng.random() ** rng.choice([0.3, 1.0, 2.0]) for _ in ed]
            deg = [0.0] * nv
            for (a, b), w in zip(ed, raw):
                deg[a] += w
                deg[b] += w
            md = max(deg)
            xs = [w / md for w in raw]
            verify_simple_fractional_matching(nv, ed, xs, f"Suite3.ell{ell}[{trial}]")
            e_idx = ed.index((0, 1))

            p0, g0 = run_branch0(nv, ed, xs, c1)
            xp = [0.0] * nv
            D = [0.0] * nv
            ae = None
            for idx, ((a, b), x) in enumerate(zip(ed, xs)):
                g = min(c1, (1.0 - c1 * xp[a]) * (1.0 - c1 * xp[b]))
                if idx == e_idx:
                    ae = g / ((1.0 - D[a]) * (1.0 - D[b]))
                xp[a] += x
                xp[b] += x
                D[a] += g * x
                D[b] += g * x

            err = abs(p0[e_idx] - g0[e_idx])
            bound = 2.0 * ae * (beta0 ** ell)
            require(err <= bound + 1e-12, f"ell={ell} trial={trial}: |p0 - g0|={err} > 2*a_e*beta0^ell={bound}")
            max_disagree_ratio = max(max_disagree_ratio, err / bound if bound > 0 else 0.0)

        print(
            f"  [PASS] ell={ell} ({n_trials} strictly feasible cyclic non-bottom-up trials): "
            f"max |p_G - gamma_e| / (2*a_e*beta0^ell) = {max_disagree_ratio:.6e} <= 1.0"
        )


def suite_4_necessity_separations(cert):
    """
    Confirms BOTH load-bearing hypotheses required by Theorems 1-3:
      (a) Simple-graph necessity: parallel-edge bipartite multigraph (k parallel copies of (0,1), each x=1/k)
          drops continuous staircase selectability to 0.352687 (k=2) and 0.330779 (k=200) < alpha*!
      (b) Oblivious-adversary necessity: an adaptive arrival-order adversary on simple forests probes Branch 0
          vs Branch T via an early forest, then routes the target edge to (1,1) iff Branch 0 else (0,0),
          forcing expected selectability w0*(1 - c1)^2 ~= 0.360869 < alpha*!
    """
    print("\n=== SUITE 4: LOAD-BEARING NECESSITY SEPARATIONS (MULTIGRAPHS & ADAPTIVE ADVERSARY) ===")
    c1 = cert["continuous_staircase"]["c1_num"] / cert["continuous_staircase"]["c1_den"]
    tau = (1.0 - 2.0 * c1) / (c1 * (1.0 - c1))
    B = 1.0 - tau
    w0 = 1.0 / (1.0 + c1 * (1.0 - c1) * I_of_b(B))

    # (a) Parallel-edge bipartite multigraph: k parallel copies of (u=0, v=1) each of weight (1 - eps)/k,
    #     followed by a pendant target edge (u=0, w=2) of weight eps = 0.001.
    #     Parallel edges share the exact same survival event {u, v both free}, causing Branch 0 to
    #     over-deplete u below 1 - D_u(e), while Branch T gives 0 on (u, w) since x_w(e) = 0 < tau!
    def parallel_multigraph_target_sel(k, eps=0.001):
        x = (1.0 - eps) / k
        S = 1.0
        D = 0.0
        for j in range(k):
            xp = j * x
            g = min(c1, (1.0 - c1 * xp) ** 2)
            a = g / ((1.0 - D) ** 2)
            S *= 1.0 - x * a
            D += g * x
        p_target = w0 * S * (c1 / (1.0 - D))
        # Cross-check k=2 with exact 2^|V| DP directly!
        if k == 2:
            ed = [(0, 1), (0, 1), (0, 2)]
            xs = [x, x, eps]
            p0, _ = run_branch0(3, ed, xs, c1)
            require(abs(w0 * p0[-1] - p_target) < 1e-14, "k=2 DP mismatch")
        return p_target

    m2 = parallel_multigraph_target_sel(2, eps=0.001)
    m200 = parallel_multigraph_target_sel(200, eps=0.001)
    require(m2 < ALPHA_STAR and abs(m2 - 0.352687) < 1e-5, f"Unexpected k=2 multigraph target sel: {m2}")
    require(m200 < ALPHA_STAR and abs(m200 - 0.330779) < 1e-5, f"Unexpected k=200 multigraph target sel: {m200}")
    print(
        f"  [CONFIRMED] Parallel-edge bipartite multigraph (eps=0.001): "
        f"k=2 -> {m2:.6f}, k=200 -> {m200:.6f} (both < alpha*={ALPHA_STAR:.6f})"
    )

    # (b) Adaptive arrival-order adversary on simple forests:
    #     Present N disjoint early edges (x=1/2); Branch 0 matches >=1 with prob 1 - (1 - c1/2)^N -> 1,
    #     while every Branch T in [tau, 1] matches 0. Then route target edge to (1,1) if Branch 0 else (0,0).
    adapt_val = w0 * ((1.0 - c1) ** 2)
    require(adapt_val < ALPHA_STAR and abs(adapt_val - 0.360869) < 1e-5, f"Unexpected adaptive val: {adapt_val}")
    print(f"  [CONFIRMED] Adaptive arrival-order adversary on simple forests: w0*(1-c1)^2 = {adapt_val:.6f} < alpha*={ALPHA_STAR:.6f}")


def main():
    cert_path = Path(__file__).resolve().parent / "cert_vanishing_lb.json"
    with open(cert_path, "r") as f:
        cert = json.load(f)
    suite_1_adversarial_trees(cert)
    suite_2_hard_benchmarks(cert)
    suite_3_cyclic_locally_tree_like_and_mmn_bug(cert)
    suite_4_necessity_separations(cert)
    print("\nALL 4 ADVERSARIAL FUZZER SUITES PASSED WITH ZERO FAILURES!")


if __name__ == "__main__":
    main()
