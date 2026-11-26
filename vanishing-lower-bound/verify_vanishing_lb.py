#!/usr/bin/env python3
"""
Standalone Cleanroom Verifier for Strongly-Online Global-Coin Matching OCRS Lower Bounds.
Zero shared code with generator; 100% python3 -O safe (zero assert-based checks).

Stage-2 audit patches incorporated:
  [V1] Exact branch-count binding (`len(branches) == 1` for `discrete_2branch`,
       `len(branches) == 64` for `discrete_64branch`).
  [V2] All mpmath formatting runs inside `with mp.workdps(100):` (no 53-bit truncation).
  [V3] `corner_iv` clamps Branch 0's contribution via `min(c1_iv, (1 - c1_iv)*(1 - c1_iv*t_hi_iv))`
       so sub-tau cells cannot inflate `corner_iv`.
  [V4] Precision isolation (`iv.dps = 100` and `mp.workdps(100)`) enforced inside
       `verify_certificate_dict()` and `verify_beta_contraction()` regardless of caller state.

Verifies three strongly-online global-coin OCRSs beating alpha* = (3 - sqrt(5))/2:
  1. `continuous_staircase`: BOTH pure exact integer/Fraction quadratic-irrational proof
     AND 100-dps outward-rounded interval arithmetic (`mpmath.iv`).
  2. `discrete_2branch`: 100-dps outward-rounded interval arithmetic + zero-gap 2D box partition
     (8 sub-strips per threshold cell).
  3. `discrete_64branch`: 100-dps outward-rounded interval arithmetic + zero-gap 2D box partition
     (8 sub-strips per threshold cell).
"""
from contextlib import contextmanager
from fractions import Fraction
import json
from pathlib import Path
import sys
import mpmath as mp
from mpmath import iv

TARGET_DPS = 100


class VerificationError(Exception):
    pass


@contextmanager
def high_precision(dps=TARGET_DPS):
    """[V4] Isolate both `mpmath.mp.dps` and `mpmath.iv.dps` against caller mutation."""
    old_iv_dps = iv.dps
    iv.dps = dps
    try:
        with mp.workdps(dps):
            yield
    finally:
        iv.dps = old_iv_dps


def require(cond, msg):
    if not cond:
        raise VerificationError(msg)


def require_exact_keys(obj, expected_keys, ctx):
    require(isinstance(obj, dict), f"{ctx}: expected dict")
    got = set(obj.keys())
    exp = set(expected_keys)
    if got != exp:
        raise VerificationError(
            f"{ctx}: key mismatch (missing={sorted(exp - got)}, extra={sorted(got - exp)})"
        )


def require_int(val, ctx, min_val=None):
    require(type(val) is int, f"{ctx}: expected strict int, got {type(val).__name__}")
    if min_val is not None:
        require(val >= min_val, f"{ctx}: {val} < {min_val}")
    return val


def frac_to_iv(fr: Fraction):
    return iv.mpf(fr.numerator) / iv.mpf(fr.denominator)


def A_of_b_iv(b_iv):
    """Outward-rounded 100-dps enclosure of A(b) = 2 / (2b + 1 + sqrt(4b + 1)),
    independently verified against the defining quadratic A = (1 - b*A)^2 and 0 < A <= 1."""
    require(b_iv.a > 0 and b_iv.b <= 1, f"b_iv={b_iv} out of (0, 1]")
    A_iv = iv.mpf(2) / (iv.mpf(2) * b_iv + iv.mpf(1) + iv.sqrt(iv.mpf(4) * b_iv + iv.mpf(1)))
    require(A_iv.a > 0 and A_iv.b <= 1, f"A_iv={A_iv} out of (0, 1]")
    resid = A_iv - (iv.mpf(1) - b_iv * A_iv) ** 2
    require(resid.a <= 0 <= resid.b, f"A(b) quadratic residual {resid} does not enclose 0!")
    return A_iv


def verify_beta_contraction(D_max_frac: Fraction, beta_ub_frac: Fraction, ctx: str):
    """
    Verifies beta(D_max) = (1/(1 - D_max)) * ln(1/(1 - D_max)) < beta_ub < 1
    via TWO independent methods:
      (M1) 100-dps outward-rounded interval arithmetic (`mpmath.iv`)
      (M2) Pure exact rational Taylor upper bound (80 terms + geometric remainder, zero mpmath).
    """
    with high_precision(TARGET_DPS):
        require(Fraction(0) < D_max_frac < Fraction(1, 2), f"{ctx}: D_max={D_max_frac} not in (0, 1/2)")
        require(Fraction(0) < beta_ub_frac < Fraction(1), f"{ctx}: beta_ub={beta_ub_frac} not in (0, 1)")

        # Method 1: 100-dps outward-rounded interval arithmetic
        D_iv = frac_to_iv(D_max_frac)
        one_m_D = iv.mpf(1) - D_iv
        beta_iv = -iv.log(one_m_D) / one_m_D
        beta_ub_iv = frac_to_iv(beta_ub_frac)
        require(beta_iv.b < beta_ub_iv.a, f"{ctx}: [iv] beta_iv.b={beta_iv.b} >= beta_ub={beta_ub_frac}")
        require(beta_iv.b < 1, f"{ctx}: [iv] beta_iv.b={beta_iv.b} >= 1")

        # Method 2: Pure exact rational Taylor upper bound (N = 80 terms)
        N_terms = 80
        d = D_max_frac
        d_pow = d
        ln_ub = Fraction(0)
        for k in range(1, N_terms + 1):
            ln_ub += d_pow / k
            d_pow *= d
        rem_ub = d_pow / ((N_terms + 1) * (1 - d))
        ln_ub += rem_ub
        beta_rat_ub = ln_ub / (1 - d)
        require(
            beta_rat_ub < beta_ub_frac,
            f"{ctx}: [exact rational Taylor] beta_rat_ub={float(beta_rat_ub):.12f} >= beta_ub={float(beta_ub_frac):.12f}",
        )
        require(beta_rat_ub < Fraction(1), f"{ctx}: [exact rational Taylor] beta_rat_ub >= 1")
        return beta_iv


def verify_gt_alpha_star_exact(c_frac: Fraction, ctx: str):
    """Pure exact integer check that c_frac > (3 - sqrt(5))/2, cross-checked at 100-dps iv."""
    rhs = Fraction(3) - 2 * c_frac
    if rhs > 0:
        require(Fraction(5) > rhs * rhs, f"{ctx}: c={c_frac} is NOT strictly greater than (3-sqrt(5))/2!")
    alpha_iv = (iv.mpf(3) - iv.sqrt(iv.mpf(5))) / iv.mpf(2)
    c_iv = frac_to_iv(c_frac)
    require(c_iv.a > alpha_iv.b, f"{ctx}: [iv] c_iv.a={c_iv.a} <= alpha*.b={alpha_iv.b}")
    return c_iv - alpha_iv


def verify_continuous_staircase(cfg: dict):
    ctx = "continuous_staircase"
    require_exact_keys(
        cfg,
        ["c1_num", "c1_den", "claimed_c_lb_num", "claimed_c_lb_den", "claimed_beta_ub_num", "claimed_beta_ub_den"],
        ctx,
    )
    c1 = Fraction(require_int(cfg["c1_num"], f"{ctx}.c1_num", 1), require_int(cfg["c1_den"], f"{ctx}.c1_den", 1))
    c_lb = Fraction(
        require_int(cfg["claimed_c_lb_num"], f"{ctx}.claimed_c_lb_num", 1),
        require_int(cfg["claimed_c_lb_den"], f"{ctx}.claimed_c_lb_den", 1),
    )
    beta_ub = Fraction(
        require_int(cfg["claimed_beta_ub_num"], f"{ctx}.claimed_beta_ub_num", 1),
        require_int(cfg["claimed_beta_ub_den"], f"{ctx}.claimed_beta_ub_den", 1),
    )

    require(c1 < Fraction(1, 2), f"{ctx}: c1={c1} >= 1/2")
    verify_gt_alpha_star_exact(c1, f"{ctx}.c1")
    tau = (1 - 2 * c1) / (c1 * (1 - c1))
    B = 1 - tau
    require(Fraction(0) < tau < Fraction(1), f"{ctx}: tau={tau} not in (0, 1)")
    require(Fraction(0) < B < Fraction(1), f"{ctx}: B={B} not in (0, 1)")
    require((1 - c1 * tau) * (1 - c1) == c1, f"{ctx}: boundary identity failed!")

    # Method A: Pure exact integer/Fraction quadratic-irrational proof
    R = 1 + 4 * B
    P = 1 + c1 * (1 - c1) * (6 * B * B + 6 * B - 1) / 12
    Q = c1 * (1 - c1) * R / 12
    require(Q > 0 and R > 1, f"{ctx}: invalid Q, R")
    lhs = c1 / c_lb - P
    require(lhs > 0, f"{ctx}: [exact int] c_lb={c_lb} exceeds exact continuous staircase value (lhs <= 0)")
    require(lhs * lhs > Q * Q * R, f"{ctx}: [exact int] c_lb={c_lb} exceeds exact continuous staircase value!")

    # Method B: Independent 100-dps outward-rounded interval arithmetic
    c1_iv = frac_to_iv(c1)
    B_iv = frac_to_iv(B)
    A_B_iv = A_of_b_iv(B_iv)
    I1_iv = (
        iv.mpf(1) / (iv.mpf(2) * A_B_iv ** 2)
        - iv.mpf(1) / (iv.mpf(3) * A_B_iv * iv.sqrt(A_B_iv))
        - iv.mpf(1) / iv.mpf(6)
    )
    I2_iv = (
        iv.mpf(6) * B_iv ** 2
        + iv.mpf(6) * B_iv
        - iv.mpf(1)
        + (iv.mpf(1) + iv.mpf(4) * B_iv) * iv.sqrt(iv.mpf(1) + iv.mpf(4) * B_iv)
    ) / iv.mpf(12)
    diff_I = I1_iv - I2_iv
    require(diff_I.a <= 0 <= diff_I.b, f"{ctx}: two closed forms for I(B) disagree: {diff_I}")

    D_late_iv = B_iv * A_B_iv
    require(D_late_iv.b < c1_iv.a, f"{ctx}: late branch D_max {D_late_iv} >= c1 {c1_iv}")

    w0_iv = iv.mpf(1) / (iv.mpf(1) + c1_iv * (iv.mpf(1) - c1_iv) * I2_iv)
    c_exact_iv = w0_iv * c1_iv
    c_lb_iv = frac_to_iv(c_lb)
    require(c_exact_iv.a > c_lb_iv.b, f"{ctx}: [iv] c_exact.a={c_exact_iv.a} <= claimed_c_lb={c_lb}")

    margin_iv = verify_gt_alpha_star_exact(c_lb, ctx)
    beta_iv = verify_beta_contraction(c1, beta_ub, ctx)

    return {
        "c_exact_iv": c_exact_iv,
        "c_lb": c_lb,
        "margin_iv": margin_iv,
        "beta_iv": beta_iv,
        "w0_iv": w0_iv,
        "tau": tau,
        "B": B,
    }


def verify_discrete_mixture(cfg: dict, ctx: str, expected_k: int):
    require_exact_keys(
        cfg,
        [
            "c1_num",
            "c1_den",
            "weight_den",
            "w0_num",
            "branches",
            "claimed_c_lb_num",
            "claimed_c_lb_den",
            "claimed_beta_ub_num",
            "claimed_beta_ub_den",
        ],
        ctx,
    )
    c1 = Fraction(require_int(cfg["c1_num"], f"{ctx}.c1_num", 1), require_int(cfg["c1_den"], f"{ctx}.c1_den", 1))
    W_den = require_int(cfg["weight_den"], f"{ctx}.weight_den", 1)
    w0_num = require_int(cfg["w0_num"], f"{ctx}.w0_num", 1)
    c_lb = Fraction(
        require_int(cfg["claimed_c_lb_num"], f"{ctx}.claimed_c_lb_num", 1),
        require_int(cfg["claimed_c_lb_den"], f"{ctx}.claimed_c_lb_den", 1),
    )
    beta_ub = Fraction(
        require_int(cfg["claimed_beta_ub_num"], f"{ctx}.claimed_beta_ub_num", 1),
        require_int(cfg["claimed_beta_ub_den"], f"{ctx}.claimed_beta_ub_den", 1),
    )
    branches = cfg["branches"]
    require(isinstance(branches, list), f"{ctx}: branches must be a list")
    # [V1] Exact branch-count binding
    require(
        len(branches) == expected_k,
        f"{ctx}: expected exactly {expected_k} late branch(es), got {len(branches)}",
    )

    require(c1 < Fraction(1, 2), f"{ctx}: c1 >= 1/2")
    verify_gt_alpha_star_exact(c1, f"{ctx}.c1")
    tau = (1 - 2 * c1) / (c1 * (1 - c1))
    require(Fraction(0) < tau < Fraction(1), f"{ctx}: tau={tau} out of (0, 1)")
    require((1 - c1 * tau) * (1 - c1) == c1, f"{ctx}: tau identity failed")

    w_sum = w0_num
    t_list = []
    w_list = []
    for j, br in enumerate(branches):
        bctx = f"{ctx}.branches[{j}]"
        require_exact_keys(br, ["t_num", "t_den", "w_num"], bctx)
        tj = Fraction(require_int(br["t_num"], f"{bctx}.t_num", 0), require_int(br["t_den"], f"{bctx}.t_den", 1))
        wj = require_int(br["w_num"], f"{bctx}.w_num", 1)
        require(Fraction(0) < tj < Fraction(1), f"{bctx}: threshold {tj} not in (0, 1)")
        if j == 0:
            require(tj <= tau, f"{bctx}: first threshold {tj} > tau={tau} leaves gap [tau, t_0) uncovered!")
        else:
            require(tj > t_list[-1], f"{bctx}: thresholds not strictly increasing ({tj} <= {t_list[-1]})")
            if j == 1:
                # [V3] Prevent dead sub-tau cell [t_0, t_1] <= tau from entering threshold partition
                require(tj > tau, f"{bctx}: second threshold {tj} <= tau={tau} creates invalid sub-tau cell!")
        t_list.append(tj)
        w_list.append(Fraction(wj, W_den))
        w_sum += wj
    require(w_sum == W_den, f"{ctx}: exact weight sum {w_sum} != weight_den {W_den}")

    w0_iv = frac_to_iv(Fraction(w0_num, W_den))
    c1_iv = frac_to_iv(c1)
    c_lb_iv = frac_to_iv(c_lb)

    # Region 0: s = min(x, y) in [0, t_0) subset [0, tau) -> gamma_0(x, y) = c1
    reg0_iv = w0_iv * c1_iv
    require(reg0_iv.a > c_lb_iv.b, f"{ctx}: Region 0 value {reg0_iv.a} <= claimed_c_lb {c_lb}")
    worst_iv = reg0_iv

    # Regions j = 0..K-1: s in [t_j, t_{j+1}] (with t_K = 1)
    cum_boost_iv = iv.mpf(0)
    K = len(t_list)
    for j in range(K):
        tj = t_list[j]
        tj_hi = t_list[j + 1] if j + 1 < K else Fraction(1)
        b_iv = frac_to_iv(1 - tj)
        Aj_iv = A_of_b_iv(b_iv)
        require((b_iv * Aj_iv).b < c1_iv.a, f"{ctx}: branch {j} D_max >= c1")
        cum_boost_iv = cum_boost_iv + frac_to_iv(w_list[j]) * Aj_iv

        # Check (1) [V3]: Exact monotone corner (t_{j+1}, 1) clamped by min(c1, ...)
        t_hi_iv = frac_to_iv(tj_hi)
        raw_corner_iv = (iv.mpf(1) - c1_iv) * (iv.mpf(1) - c1_iv * t_hi_iv)
        g0_corner_iv = iv.mpf([min(c1_iv.a, raw_corner_iv.a), min(c1_iv.b, raw_corner_iv.b)])
        corner_iv = w0_iv * g0_corner_iv + cum_boost_iv
        require(corner_iv.a > c_lb_iv.b, f"{ctx}: cell {j} corner {corner_iv.a} <= claimed_c_lb {c_lb}")

        # Check (2): Independent 2D interval box check over [t_j, t_{j+1}] x [t_j, 1]
        # subdivided into 8 sub-strips in s
        n_sub = 8
        for m in range(n_sub):
            sa = tj + (tj_hi - tj) * Fraction(m, n_sub)
            sb = tj + (tj_hi - tj) * Fraction(m + 1, n_sub)
            s_box = iv.mpf([frac_to_iv(sa).a, frac_to_iv(sb).b])
            r_box = iv.mpf([frac_to_iv(sa).a, iv.mpf(1)])
            prod_box = (iv.mpf(1) - c1_iv * s_box) * (iv.mpf(1) - c1_iv * r_box)
            g0_lb = min(c1_iv.a, prod_box.a)
            box_lb = w0_iv.a * g0_lb + cum_boost_iv.a
            require(box_lb > c_lb_iv.b, f"{ctx}: cell {j} sub-box {m} lb {box_lb} <= claimed_c_lb {c_lb}")

        if corner_iv.a < worst_iv.a:
            worst_iv = corner_iv

    margin_iv = verify_gt_alpha_star_exact(c_lb, ctx)
    beta_iv = verify_beta_contraction(c1, beta_ub, ctx)
    return {
        "worst_iv": worst_iv,
        "c_lb": c_lb,
        "margin_iv": margin_iv,
        "beta_iv": beta_iv,
        "num_branches": 1 + K,
    }


def verify_certificate_dict(cert: dict):
    """[V4] Enforces 100-dps precision internally regardless of caller precision state."""
    with high_precision(TARGET_DPS):
        require_exact_keys(
            cert,
            ["schema_version", "continuous_staircase", "discrete_2branch", "discrete_64branch"],
            "root",
        )
        require(
            cert["schema_version"] == "ocrs-vanishing-lb-r3-v1",
            f"unexpected schema_version: {cert['schema_version']!r}",
        )
        r_cont = verify_continuous_staircase(cert["continuous_staircase"])
        r_d2 = verify_discrete_mixture(cert["discrete_2branch"], "discrete_2branch", expected_k=1)
        r_d64 = verify_discrete_mixture(cert["discrete_64branch"], "discrete_64branch", expected_k=64)
        return {"continuous_staircase": r_cont, "discrete_2branch": r_d2, "discrete_64branch": r_d64}


def fmt_100dps(x, digits=18):
    """[V2] Format ivmpf endpoint at full 100-dps precision without 53-bit truncation."""
    with mp.workdps(TARGET_DPS):
        return mp.nstr(mp.mpf(x), digits)


def main():
    default_path = Path(__file__).resolve().parent / "cert_vanishing_lb.json"
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else default_path
    with open(path, "r") as f:
        cert = json.load(f)
    res = verify_certificate_dict(cert)
    print(f"VERIFIED CLEANROOM CERTIFICATE: {path} (100-dps outward-rounded iv + exact integer cross-check)")
    rc = res["continuous_staircase"]
    print(
        f"  [PASS] continuous_staircase (c1=49/125): c_exact in "
        f"[{fmt_100dps(rc['c_exact_iv'].a, 20)}, {fmt_100dps(rc['c_exact_iv'].b, 20)}]\n"
        f"         c_lb = {float(rc['c_lb']):.12f} > alpha* (margin > +{fmt_100dps(rc['margin_iv'].a, 15)}), "
        f"beta0 in [{fmt_100dps(rc['beta_iv'].a, 12)}, {fmt_100dps(rc['beta_iv'].b, 12)}] < 1"
    )
    for k in ["discrete_2branch", "discrete_64branch"]:
        rk = res[k]
        print(
            f"  [PASS] {k:18s} ({rk['num_branches']:2d} branches): worst >= {fmt_100dps(rk['worst_iv'].a, 18)}, "
            f"c_lb = {float(rk['c_lb']):.12f} > alpha* (margin > +{fmt_100dps(rk['margin_iv'].a, 15)}), "
            f"beta0 <= {fmt_100dps(rk['beta_iv'].b, 12)} < 1"
        )


if __name__ == "__main__":
    main()
