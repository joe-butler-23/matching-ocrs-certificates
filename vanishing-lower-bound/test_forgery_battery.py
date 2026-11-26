#!/usr/bin/env python3
"""
Fired-Once Forgery & Stage-2 Audit Regression Battery for `verify_vanishing_lb.py` (22 checks total).
Must be executed under `python3 -O`:
  - 21 distinct planted forgeries (including Stage-2 [V1] branch-count binding & [V3] sub-tau cell clamp)
    must ALL be rejected with `VerificationError`.
  - 1 Stage-2 [V4] caller-precision-isolation test (caller mutates `iv.dps = 6; mp.mp.dps = 6`)
    must still pass cleanly at full 100-dps precision and restore caller precision.
"""
import copy
from fractions import Fraction
import json
from pathlib import Path
import sys
import mpmath as mp
from mpmath import iv
import verify_vanishing_lb as verifier


def run_forgery_battery(cert_path=None):
    if cert_path is None:
        cert_path = Path(__file__).resolve().parent / "cert_vanishing_lb.json"
    with open(cert_path, "r") as f:
        good = json.load(f)

    # Confirm genuine certificate passes first
    verifier.verify_certificate_dict(good)

    forgeries = []

    def add(name, mut_fn):
        c = copy.deepcopy(good)
        mut_fn(c)
        forgeries.append((name, c))

    # F01-F06: continuous_staircase math & contraction attacks
    add(
        "F01_cont_overclaim_1e-9",
        lambda c: c["continuous_staircase"].update(
            claimed_c_lb_num=c["continuous_staircase"]["claimed_c_lb_num"] + 1000
        ),
    )
    add("F02_cont_c1_ge_half", lambda c: c["continuous_staircase"].update(c1_num=1, c1_den=2))
    add(
        "F03_cont_c1_above_W1_overclaim",
        lambda c: c["continuous_staircase"].update(
            c1_num=44, c1_den=100, claimed_c_lb_num=382000000000, claimed_beta_ub_num=999999999
        ),
    )
    add("F04_cont_c1_below_alpha_star", lambda c: c["continuous_staircase"].update(c1_num=38, c1_den=100))
    add("F05_cont_beta_ub_too_small", lambda c: c["continuous_staircase"].update(claimed_beta_ub_num=810000000))
    add("F06_cont_beta_ub_ge_1", lambda c: c["continuous_staircase"].update(claimed_beta_ub_num=1000000000))

    # F07-F10: strict schema & type attacks
    add("F07_extra_root_key", lambda c: c.update(extra_malicious=123))
    add("F08_wrong_schema_version", lambda c: c.update(schema_version="ocrs-vanishing-lb-r3-v999"))
    add("F09_bool_as_int", lambda c: c["continuous_staircase"].update(c1_num=True))
    add("F10_zero_denominator", lambda c: c["continuous_staircase"].update(c1_den=0))

    # F11-F14: discrete_2branch attacks
    add(
        "F11_d2_weights_sum_gt_1",
        lambda c: c["discrete_2branch"].update(w0_num=c["discrete_2branch"]["w0_num"] + 1),
    )
    add("F12_d2_negative_branch_weight", lambda c: c["discrete_2branch"]["branches"][0].update(w_num=-1))
    add("F13_d2_gap_at_tau", lambda c: c["discrete_2branch"]["branches"][0].update(t_num=99, t_den=100))
    add(
        "F14_d2_overclaim_1e-8",
        lambda c: c["discrete_2branch"].update(claimed_c_lb_num=c["discrete_2branch"]["claimed_c_lb_num"] + 10000),
    )

    # F15-F18: discrete_64branch attacks
    def swap_5_6(c):
        b = c["discrete_64branch"]["branches"]
        b[5], b[6] = b[6], b[5]

    add("F15_d64_unsorted_thresholds", swap_5_6)
    add("F16_d64_threshold_ge_1", lambda c: c["discrete_64branch"]["branches"][-1].update(t_num=1, t_den=1))
    add("F17_d64_extra_branch_key", lambda c: c["discrete_64branch"]["branches"][0].update(sneaky=0))
    add(
        "F18_d64_overclaim_1e-9",
        lambda c: c["discrete_64branch"].update(claimed_c_lb_num=c["discrete_64branch"]["claimed_c_lb_num"] + 1000),
    )

    # F20 [V1 Regression]: Stuff 64-branch payload into `discrete_2branch`!
    add("F20_V1_64branch_stuffed_into_d2", lambda c: c.update(discrete_2branch=copy.deepcopy(c["discrete_64branch"])))

    # F21 [V3 Regression]: Sub-tau cell (t_0 = 1/2 < tau, t_1 = 3/4 < tau) with valid claimed_c_lb
    def sub_tau_cell_attack(c):
        d64 = c["discrete_64branch"]
        d64["branches"][0]["t_num"] = 1
        d64["branches"][0]["t_den"] = 2
        d64["branches"][1]["t_num"] = 3
        d64["branches"][1]["t_den"] = 4

    add("F21_V3_sub_tau_cell_rejected", sub_tau_cell_attack)

    rejected = 0
    for name, bad_cert in forgeries:
        try:
            verifier.verify_certificate_dict(bad_cert)
            print(f"  [FAIL — ACCEPTED FORGERY!!] {name}")
            sys.exit(1)
        except verifier.VerificationError as e:
            rejected += 1
            print(f"  [REJECTED {rejected:2d}/21] {name:34s} -> {e}")

    # F19: Direct unit forgery for W(1) contraction boundary D_max = 433/1000 > 1 - W(1)
    try:
        verifier.verify_beta_contraction(Fraction(433, 1000), Fraction(999999999, 1000000000), "F19_direct_W1")
        print("  [FAIL — ACCEPTED W(1) CONTRACTION FORGERY!!]")
        sys.exit(1)
    except verifier.VerificationError as e:
        rejected += 1
        print(f"  [REJECTED {rejected:2d}/21] {'F19_direct_W1_contraction_gt_1':34s} -> {e}")

    if rejected != 21:
        print(f"  [FAIL] Expected 21 rejected forgeries, got {rejected}")
        sys.exit(1)

    # T22 [V4 Regression]: Caller mutates global precision down to 6 dps -> verifier MUST still pass at 100 dps!
    old_iv = iv.dps
    old_mp = mp.mp.dps
    try:
        iv.dps = 6
        mp.mp.dps = 6
        res = verifier.verify_certificate_dict(good)
        # Verify caller precision was restored cleanly AND output string has full 18-digit accuracy!
        if iv.dps != 6 or mp.mp.dps != 6:
            print(f"  [FAIL — T22_V4] Precision context leaked (iv.dps={iv.dps}, mp.dps={mp.mp.dps})")
            sys.exit(1)
        s_val = verifier.fmt_100dps(res["continuous_staircase"]["c_exact_iv"].a, 18)
        if s_val != "0.382673618863030223":
            print(f"  [FAIL — T22_V4/V2] Expected '0.382673618863030223', got {s_val!r}")
            sys.exit(1)
        print(f"  [PASS    22/22] {'T22_V4_V2_caller_dps6_isolation':34s} -> 100-dps c_exact={s_val}, restored dps=6")
    finally:
        iv.dps = old_iv
        mp.mp.dps = old_mp

    print("\nFORGERY & STAGE-2 REGRESSION BATTERY COMPLETE: 21/21 forgeries REJECTED + [V4]/[V2] isolation PASSED (22/22)!")


if __name__ == "__main__":
    run_forgery_battery()
