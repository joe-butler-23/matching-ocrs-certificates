#!/usr/bin/env python3
"""
Standalone generator for the strongly-online global-coin matching OCRS certificate (`cert_vanishing_lb.json`).
Zero shared code with `verify_vanishing_lb.py`.

Generates exact rational parameters for three strongly-online global-coin OCRSs
that strictly beat alpha* = (3 - sqrt(5))/2 = 0.38196601125010515...:
  1. `continuous_staircase`: 1 rational parameter c1 = 49/125 -> c_lb = 0.382673618862
  2. `discrete_2branch`:     1 biased coin flip (2 branches)  -> c_lb = 0.382309325242
  3. `discrete_64branch`:    64 rational late thresholds      -> c_lb = 0.382661888422
"""
import json
from fractions import Fraction
from pathlib import Path
import mpmath as mp

mp.mp.dps = 120


def A_mp(b):
    if b == 0:
        return mp.mpf(1)
    return mp.mpf(2) / (2 * b + 1 + mp.sqrt(4 * b + 1))


def build_continuous_staircase():
    c1 = Fraction(49, 125)
    tau = (1 - 2 * c1) / (c1 * (1 - c1))
    B = 1 - tau
    c1_m = mp.mpf(c1.numerator) / mp.mpf(c1.denominator)
    B_m = mp.mpf(B.numerator) / mp.mpf(B.denominator)
    I_B = (6 * B_m * B_m + 6 * B_m - 1 + (1 + 4 * B_m) * mp.sqrt(1 + 4 * B_m)) / 12
    w0 = 1 / (1 + c1_m * (1 - c1_m) * I_B)
    c_val = w0 * c1_m
    beta0 = -mp.log(1 - c1_m) / (1 - c1_m)
    scale = 10 ** 12
    c_lb_num = int(mp.floor(c_val * scale)) - 1
    beta_ub_num = int(mp.ceil(beta0 * 10 ** 9)) + 1
    return {
        "c1_num": c1.numerator,
        "c1_den": c1.denominator,
        "claimed_c_lb_num": c_lb_num,
        "claimed_c_lb_den": scale,
        "claimed_beta_ub_num": beta_ub_num,
        "claimed_beta_ub_den": 10 ** 9,
    }


def build_discrete_k_branch(K_late, c1_frac):
    tau = (1 - 2 * c1_frac) / (c1_frac * (1 - c1_frac))
    B = 1 - tau
    ts = [tau + B * Fraction(j, K_late) for j in range(K_late)]
    ts_next = [tau + B * Fraction(j + 1, K_late) for j in range(K_late)]
    c1_m = mp.mpf(c1_frac.numerator) / mp.mpf(c1_frac.denominator)

    r_m = []
    for j in range(K_late):
        dt = mp.mpf((ts_next[j] - ts[j]).numerator) / mp.mpf((ts_next[j] - ts[j]).denominator)
        b_j = mp.mpf((1 - ts[j]).numerator) / mp.mpf((1 - ts[j]).denominator)
        Aj = A_mp(b_j)
        r_m.append(c1_m * (1 - c1_m) * dt / Aj)
    w0_m = 1 / (1 + sum(r_m))
    W_DEN = 10 ** 15
    w_late_nums = [int(mp.ceil(w0_m * rj * W_DEN)) for rj in r_m]
    w0_num = W_DEN - sum(w_late_nums)
    if w0_num <= 0:
        raise RuntimeError("Non-positive w0_num")

    w0_exact = mp.mpf(w0_num) / W_DEN
    min_val = w0_exact * c1_m
    cum_boost = mp.mpf(0)
    for j in range(K_late):
        b_j = mp.mpf((1 - ts[j]).numerator) / mp.mpf((1 - ts[j]).denominator)
        Aj = A_mp(b_j)
        cum_boost += (mp.mpf(w_late_nums[j]) / W_DEN) * Aj
        t_hi = mp.mpf(ts_next[j].numerator) / mp.mpf(ts_next[j].denominator)
        val_j = w0_exact * min(c1_m, (1 - c1_m) * (1 - c1_m * t_hi)) + cum_boost
        if val_j < min_val:
            min_val = val_j

    scale = 10 ** 12
    c_lb_num = int(mp.floor(min_val * scale)) - 1
    beta0 = -mp.log(1 - c1_m) / (1 - c1_m)
    beta_ub_num = int(mp.ceil(beta0 * 10 ** 9)) + 1

    return {
        "c1_num": c1_frac.numerator,
        "c1_den": c1_frac.denominator,
        "weight_den": W_DEN,
        "w0_num": w0_num,
        "branches": [
            {
                "t_num": ts[j].numerator,
                "t_den": ts[j].denominator,
                "w_num": w_late_nums[j],
            }
            for j in range(K_late)
        ],
        "claimed_c_lb_num": c_lb_num,
        "claimed_c_lb_den": scale,
        "claimed_beta_ub_num": beta_ub_num,
        "claimed_beta_ub_den": 10 ** 9,
    }


def main():
    cert = {
        "schema_version": "ocrs-vanishing-lb-r3-v1",
        "continuous_staircase": build_continuous_staircase(),
        "discrete_2branch": build_discrete_k_branch(1, Fraction(386767, 1000000)),
        "discrete_64branch": build_discrete_k_branch(64, Fraction(49, 125)),
    }
    out_path = Path(__file__).resolve().parent / "cert_vanishing_lb.json"
    with open(out_path, "w") as f:
        json.dump(cert, f, indent=2)
    alpha_star = (3 - mp.sqrt(5)) / 2
    for key in ["continuous_staircase", "discrete_2branch", "discrete_64branch"]:
        d = cert[key]
        clb = mp.mpf(d["claimed_c_lb_num"]) / d["claimed_c_lb_den"]
        print(
            f"{key:22s}: claimed_c_lb = {mp.nstr(clb, 15)} "
            f"(margin +{mp.nstr(clb - alpha_star, 12)}), "
            f"beta_ub = {d['claimed_beta_ub_num']/d['claimed_beta_ub_den']:.9f}"
        )
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
