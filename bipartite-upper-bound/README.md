# Track 3: Vanishing Bipartite Upper Bound ($\le 0.38602716$) & Single-Pool Suite

Self-contained unified 100-digit outward-rounded interval and pure-$\mathbb{Q}$ verifier certifying the vanishing-probability bipartite upper bound $c^*_{\mathrm{bip, vanishing}} \le 0.3860271555$ ($\max_e x_e = 10^{-12}$, extending algebraically as $\delta \to 0^+$), together with single-pool general and 293-box pure-$\mathbb{Q}$ Paley certificates and a 35-case forgery battery.

## Certified Bounds (`certs/`)

| Certificate | Graph Family | 100-dps UB (`claimed_ub`) | Pure-$\mathbb{Q}$ UB | Benchmark Comparison |
|---|---|---:|---|---|
| **`certs/cert_bipartite_t40_vanishing.json`** | 2-family Birkhoff 2-colored bipartite ($T_2=40$, $q=10^{13}+313$, $\max_e x_e = 10^{-12} \to 0^+$) | **`0.3860271554305342`** (`0.38602716`) | — | Cuts **63.0%** of $[\alpha^*, \gamma]$ (`0.39293532`) & **46.4%** of $[\alpha^*, \beta]$ (`0.38953805`) |
| `certs/cert_general_t40_fused.json` | 4-family single-pool general ($T_2=40$, `(w_i, D)` + atom) | `0.3796296804569556` (`0.37962969`) | — | Single-pool baseline (superseded by Track 2 EMA cascade `0.3786732`) |
| `certs/cert_general_paley_lean_293box.json` | 3-family Paley ($T_2=1$, $q=2{,}000{,}003$, zero `(w_i, D)`) | `0.3806937049651673` | `0.3809993558` (**293** mono / 502 non-mono $\mathbb{Q}$-boxes) | Single-pool minimal Paley baseline |
| `certs/cert_general_paley_fused_446box.json` | 4-family fused Paley ($T_2=1$, $q=2{,}000{,}003$) | `0.3806319201217787` | `0.3809749243` (214 mono / **446** non-mono $\mathbb{Q}$-boxes) | Single-pool fused Paley baseline |

See [`../SPEC.md`](../SPEC.md) §4 for the $\delta \to 0^+$ vanishing-family proof.

## Verify

```bash
# Core verifier (mpmath only; monotone & non-monotone cross-check)
uv run --with mpmath python3 -O verify_ocrs.py --non-monotone \
  certs/cert_bipartite_t40_vanishing.json \
  certs/cert_general_t40_fused.json \
  certs/cert_general_paley_lean_293box.json \
  certs/cert_general_paley_fused_446box.json

# Full 35-case forgery battery + polytope + Paley + K6/K8 suite
uv run --with mpmath --with numpy --with scipy python3 -O test_all.py
```
