#!/usr/bin/env python3
"""
Top-level verification runner for all Matching OCRS certificates and test suites.

Usage:
  uv run --with mpmath --with numpy --with scipy python3 -O verify_all.py --fast
  uv run --with mpmath python3 -O verify_all.py --verifiers-only
  uv run --with mpmath --with numpy --with scipy python3 -O verify_all.py
"""
import argparse, ast, glob, hashlib, os, subprocess, sys, time

ROOT = os.path.dirname(os.path.abspath(__file__))


def check_ast_zero_asserts() -> None:
    files = sorted(glob.glob(os.path.join(ROOT, "**", "*.py"), recursive=True))
    for path in files:
        with open(path) as f:
            tree = ast.parse(f.read(), filename=path)
        lines = [n.lineno for n in ast.walk(tree) if isinstance(n, ast.Assert)]
        if lines:
            raise RuntimeError(f"Found forbidden assert statements in {os.path.relpath(path, ROOT)}: {lines}")
    print(f"[PASS] AST check: 0 assert statements across {len(files)} Python files (python3 -O safe)")


def verify_sha256sums_if_present() -> None:
    manifest = os.path.join(ROOT, "SHA256SUMS")
    if not os.path.exists(manifest):
        return
    count = 0
    with open(manifest) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            expected, relpath = line.split(None, 1)
            full = os.path.join(ROOT, relpath)
            h = hashlib.sha256()
            with open(full, "rb") as bf:
                for chunk in iter(lambda: bf.read(65536), b""):
                    h.update(chunk)
            actual = h.hexdigest()
            if actual != expected:
                raise ValueError(f"SHA-256 mismatch for {relpath}: expected {expected}, got {actual}")
            count += 1
    print(f"[PASS] Verified {count} SHA-256 digests from SHA256SUMS")


def run_step(label: str, cmd: list[str], cwd: str) -> None:
    print(f"\n=== {label} ===", flush=True)
    t0 = time.perf_counter()
    env = dict(os.environ)
    env["PYTHONPATH"] = cwd + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    res = subprocess.run(cmd, cwd=cwd, env=env)
    if res.returncode != 0:
        raise RuntimeError(f"Step failed ({res.returncode}): {' '.join(cmd)}")
    print(f"[PASS] {label} [{time.perf_counter() - t0:.1f}s]", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true", help="Fast ~8 s check: Track 1 full + pure-Q Paley + pure-Fraction K6/K8 + Paley 100-dps")
    ap.add_argument("--verifiers-only", action="store_true", help="Run all core verifiers only (mpmath + stdlib only, zero numpy/scipy)")
    args = ap.parse_args()

    t_all = time.perf_counter()
    py = sys.executable
    d1 = os.path.join(ROOT, "vanishing-lower-bound")
    d2 = os.path.join(ROOT, "general-upper-bound")
    d3 = os.path.join(ROOT, "bipartite-upper-bound")

    verify_sha256sums_if_present()
    check_ast_zero_asserts()

    # Track 1: always run verifier (and tests unless verifiers-only)
    run_step("Track 1: Strongly-online global-coin lower bound verifier (100-dps + exact Z)", [py, "-O", "verify_vanishing_lb.py"], d1)
    if not args.verifiers_only:
        run_step("Track 1: 22-check forgery & regression battery", [py, "-O", "test_forgery_battery.py"], d1)
        run_step("Track 1: 4-suite exact 2^|V|-DP adversarial fuzzer", [py, "-O", "adversarial_fuzzer.py"], d1)

    # Track 2: pure-Q Paley & EMA verifiers
    run_step("Track 2: Pure-Q Paley 71/150-box & 9-prime Pratt tree verifier (stdlib)", [py, "-O", "verify_l3_lean4_paley_rat.py"], d2)
    run_step("Track 2: Exact Fraction top-down pyramid polytope check (stdlib)", [py, "-O", "test_ema_coupled_finite_toy.py"], d2)

    if args.fast:
        run_step("Track 2: Paley T2=1 100-dps dual-pass EMA verifier", [py, "-O", "verify_l3_ema_ocrs.py", "cert_l3_ema_paley_t1_lean.json"], d2)
        run_step("Track 2: Pure-Fraction zero-expander K6 & K8 duals", [py, "-O", "verify_pure_finite_k6_k8.py"], d2)
    elif args.verifiers_only:
        run_step(
            "Track 2: All 3 EMA certificates (100-dps monotone + raw non-monotone)",
            [py, "-O", "verify_l3_ema_ocrs.py", "cert_l3_ema_paley_t1_lean.json", "cert_l3_ema_general_t20.json", "cert_l3_ema_general_t40.json"],
            d2,
        )
        run_step(
            "Track 3: All 4 bipartite & single-pool certificates (100-dps + pure-Q + non-monotone)",
            [
                py, "-O", "verify_ocrs.py", "--non-monotone",
                "certs/cert_bipartite_t40_vanishing.json",
                "certs/cert_general_t40_fused.json",
                "certs/cert_general_paley_lean_293box.json",
                "certs/cert_general_paley_fused_446box.json",
            ],
            d3,
        )
    else:
        run_step(
            "Track 2: All 3 EMA certificates (100-dps monotone + raw non-monotone)",
            [py, "-O", "verify_l3_ema_ocrs.py", "cert_l3_ema_paley_t1_lean.json", "cert_l3_ema_general_t20.json", "cert_l3_ema_general_t40.json"],
            d2,
        )
        run_step("Track 2: Pure-Fraction zero-expander K6 & K8 duals", [py, "-O", "verify_pure_finite_k6_k8.py"], d2)
        run_step("Track 2: 18-case forgery battery & 6-suite lemma fuzzer", [py, "-O", "test_l3_ema_forgery_and_fuzzer.py"], d2)
        run_step(
            "Track 3: All 4 bipartite & single-pool certificates (100-dps + pure-Q + non-monotone)",
            [
                py, "-O", "verify_ocrs.py", "--non-monotone",
                "certs/cert_bipartite_t40_vanishing.json",
                "certs/cert_general_t40_fused.json",
                "certs/cert_general_paley_lean_293box.json",
                "certs/cert_general_paley_fused_446box.json",
            ],
            d3,
        )
        run_step("Track 3: 35-case forgery battery & 5-suite verification", [py, "-O", "test_all.py"], d3)

    print(f"\nALL CHECKS PASSED [{time.perf_counter() - t_all:.1f}s]")


if __name__ == "__main__":
    main()
