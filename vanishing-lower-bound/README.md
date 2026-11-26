# Track 1: Strongly-Online Global-Coin Matching OCRS Lower Bound ($> \frac{3-\sqrt{5}}{2}$)

Self-contained certificate, 100-digit outward-rounded interval verifier, forgery/regression suite, and exact $2^{|V|}$-DP fuzzer proving that a 1-parameter strongly-online global-coin OCRS strictly beats $\alpha^* = \frac{3-\sqrt{5}}{2} = 0.381966011250\dots$ on all simple finite forests (arbitrary non-vanishing $\mathbf{x}$) and all simple vanishing graphs (general and bipartite) against an oblivious arrival-order adversary.

## Certified Bounds (`cert_vanishing_lb.json`)

| Scheme | Global Randomness | Certified Selectability $c_{\mathrm{LB}}$ | Margin Above $\alpha^*$ | Contraction Factor $\beta_0$ |
|---|---|---:|---:|---:|
| **`continuous_staircase`** ($c_1 = 49/125$) | 1 global draw $\omega \in \{0\} \cup [\tau, 1]$ | **`0.382673618862`** (`0.382673618863030223`) | `+0.000707607612` | `0.818388810882 < 1` |
| **`discrete_2branch`** | **1 single biased coin flip** $\omega \in \{0, 1\}$ | **`0.382309325242`** | `+0.000343313992` | `0.797429879647 < 1` |
| **`discrete_64branch`** | 64 rational thresholds (65 branches) | **`0.382661888422`** | `+0.000695877172` | `0.818388810882 < 1` |

See [`../SPEC.md`](../SPEC.md) §2 for the complete self-contained proofs (Theorems 1–3, load-bearing necessity separations for parallel-edge multigraphs and adaptive adversaries, and the strictly feasible 7-edge counterexample to `arXiv:2406.14506` v3 `lem:tree_reduction_game_part_1`, L1572).

## Verify (~2 s)

```bash
uv run --with mpmath python3 -O verify_vanishing_lb.py
uv run --with mpmath python3 -O test_forgery_battery.py
python3 -O adversarial_fuzzer.py
```
