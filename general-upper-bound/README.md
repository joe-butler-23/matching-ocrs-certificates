# Track 2: General Non-Bipartite Matching OCRS Upper Bounds ($\le 0.3786732$) & Pure-Finite $K_6$/$K_8$

Self-contained certificates, dual-pass 100-digit outward-rounded interval verifier, pure-$\mathbb{Q}$ Paley verifier, and exact-`Fraction` $K_6$/$K_8$ dual verifiers establishing finite-graph upper bounds strictly below $\frac{3-\sqrt{5}}{2}$ via a zero-Azuma cubic-potential EMA recursive cascade.

## Certified Upper Bounds

| Certificate / Instance | Free Dual Numerators* | Prime $q$ | Monotone 100-dps UB | Raw Non-Monotone 100-dps UB | Pure-$\mathbb{Q}$ UB |
|---|---:|---|---:|---:|---|
| `cert_l3_ema_general_t40.json` ($T_2=40$) | 86 | $10^{13} + 313$ | **`0.378673212573978`** | **`0.378684991895836`** ($\le 0.378685$) | — |
| `cert_l3_ema_general_t20.json` ($T_2=20$) | 46 | $10^{13} + 411$ | **`0.378673518083956`** | **`0.378684999403519`** ($\le 0.378685$) | — |
| `cert_l3_ema_paley_t1_lean.json` ($T_2=1$) | 8 | $10^8 + 7 \equiv 3 \pmod 4$ | **`0.379624600475273`** | **`0.379749499490282`** ($\le 0.37975$) | **`0.380980537`** (**71** $\mathbb{Q}$-boxes); **`0.379999139`** (**150** $\mathbb{Q}$-boxes) |
| Zero-expander $K_6$ (fixed $\pi$, $\varepsilon=10^{-9}$) | — | — | — | — | **`0.389536778`** (exact `Fraction`, $< \beta = 0.3895380533$) |
| Zero-expander $K_8$ (fixed $\pi$, $\varepsilon=10^{-9}$) | — | — | — | — | **`0.388559368`** (exact `Fraction`, $< \beta = 0.3895380533$) |

*\*Each EMA JSON specifies its structural parameters (`T2, q, P, K0, N_L, L, beta_den, eps1_den, eps2_den, D_int`) and `claimed_ub` alongside the free dual-weight numerators. Complete proofs and structural obstructions are in [`../SPEC.md`](../SPEC.md) §3.*

## Verify

```bash
# Pure Python stdlib (zero third-party packages): 9-prime Pratt tree + 71/150 pure-Q boxes + pyramid polytope
python3 -O verify_l3_lean4_paley_rat.py
python3 -O test_ema_coupled_finite_toy.py

# 100-dps dual-pass EMA verifier (mpmath only; ~85 s total)
uv run --with mpmath python3 -O verify_l3_ema_ocrs.py \
  cert_l3_ema_paley_t1_lean.json cert_l3_ema_general_t20.json cert_l3_ema_general_t40.json

# Pure-Fraction K6/K8 duals + 18-case forgery battery & 6-suite lemma fuzzer
uv run --with numpy --with scipy python3 -O verify_pure_finite_k6_k8.py
PYTHONPATH=. uv run --with mpmath --with numpy --with scipy python3 -O test_l3_ema_forgery_and_fuzzer.py
```
