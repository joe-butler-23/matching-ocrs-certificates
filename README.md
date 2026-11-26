# Matching OCRS Separation Certificates ($c^*_{\mathrm{general}} < \frac{3-\sqrt{5}}{2} < c^*_{\mathrm{trees}},\, c^*_{\mathrm{vanishing}} \le 0.38603$)

Machine-checkable certificates and self-contained proofs establishing strict separations on **both sides** of the golden-ratio independence benchmark $\alpha^* = \frac{3 - \sqrt{5}}{2} = 0.381966011250\dots$ for adversarial-order edge-arrival **Online Contention Resolution Schemes (OCRSs)** on the matching polytope of simple graphs against an oblivious arrival-order adversary:

$$c^*_{\mathrm{general}} \;\le\; \mathbf{0.3786732126} \;<\; \frac{3 - \sqrt{5}}{2} = \mathbf{0.3819660113} \;<\; \mathbf{0.3826736189} \;\le\; c^*_{\mathrm{trees}},\, c^*_{\mathrm{vanishing}} \;\le\; \mathbf{0.3860271555}.$$

This refutes both halves of the matching OCRS conjecture of Ma, MacRury, and Nuti ([`arXiv:2406.14506`](https://arxiv.org/abs/2406.14506) v3, Oct 2026; also quoted in Aminian, Niazadeh, and Nuti, [`arXiv:2603.21532`](https://arxiv.org/abs/2603.21532), Oct 2026): general non-bipartite matching cannot reach $\alpha^*$, while simple finite trees (even non-vanishing) and all vanishing simple graphs (general and bipartite) strictly exceed $\alpha^*$.

The repository is organized into **three independent, self-contained tracks**:

| Track | Result | Entry Point | Verification |
|---|---|---|---|
| **Track 1: Strongly-Online Global-Coin Lower Bound** | $c \ge 0.382673618862 > \alpha^*$ on all simple finite forests (arbitrary non-vanishing $\mathbf{x}$) and $0.382673618862 - o(1)$ on all simple vanishing graphs (general & bipartite) under oblivious arrival orders (1-coin 2-branch corollary: $\ge 0.382309325242$). Includes a 6-line disagreement-coupling proof and a strictly feasible counterexample to `arXiv:2406.14506` v3 `lem:tree_reduction_game_part_1` (L1572). | **[`vanishing-lower-bound/`](vanishing-lower-bound/)** ([README](vanishing-lower-bound/README.md) · [verifier](vanishing-lower-bound/verify_vanishing_lb.py)) | 100-dps `mpmath.iv` + exact integer check + 22 forgery/regression checks + 4 exact $2^{|V|}$-DP fuzzer suites (`~2 s`) |
| **Track 2: General Non-Bipartite Upper Bounds & Pure-Finite Gems** | Cubic-potential EMA cascade certifies $c^*_{\mathrm{general}} \le 0.3786732126$ (monotone) / $0.378685$ (raw non-monotone); 8-parameter Paley 2-coset certificate certifies $\le 0.3796246005$ (`< 0.380` in 150 pure-$\mathbb{Q}$ boxes, `< 381/1000` in **71** pure-$\mathbb{Q}$ boxes with full 9-prime Pratt tree); pure-`Fraction` $2^n$ duals certify fixed-order $K_6 \le 0.389536778$ and $K_8 \le 0.388559368$. | **[`general-upper-bound/`](general-upper-bound/)** ([README](general-upper-bound/README.md) · [verifier](general-upper-bound/verify_l3_ema_ocrs.py)) | Dual-pass 100-dps `mpmath.iv` + pure `Fraction` + 18 forgery guards + 6 lemma fuzzer suites (`~3 min`) |
| **Track 3: Vanishing Bipartite Upper Bound & Single-Pool Suite** | 2-family Birkhoff 2-colored bipartite simple graph certifies $c^*_{\mathrm{bip, vanishing}} \le 0.3860271555$ ($\max_e x_e = 10^{-12}$, extending monotonically as $\delta \to 0^+$), plus hardened single-pool verifier and 293-box pure-$\mathbb{Q}$ Paley certificate. | **[`bipartite-upper-bound/`](bipartite-upper-bound/)** ([README](bipartite-upper-bound/README.md) · [verifier](bipartite-upper-bound/verify_ocrs.py)) | 100-dps `mpmath.iv` (monotone & `--non-monotone`) + pure `Fraction` + 35 forgery guards (`~1.5 min`) |

---

## Summary of Certified Bounds

Decimals are outward-rounded bounds (floor for lower bounds, ceil for upper bounds) verified in 100-digit `mpmath.iv` interval arithmetic or exact Python `fractions.Fraction`. Exact definitions and self-contained proofs are in [`SPEC.md`](SPEC.md).

| Regime (Simple Graphs, Oblivious Adversary) | Certified Lower Bound | Certified Upper Bound | Status vs. $\alpha^* = \frac{3-\sqrt{5}}{2}$ |
|---|---:|---:|---|
| **General Non-Bipartite (Non-Vanishing)** | `0.349` (MacRury–Ma 2025) | **`0.3786732126`** (`0.378685` raw non-mono) | Strictly below $\alpha^*$ (margin `0.003293`) |
| **General Non-Bipartite (8-param Paley, $T_2=1$)** | — | **`0.3796246005`** (`0.380981` in 71 $\mathbb{Q}$-boxes) | Strictly below $\alpha^*$ (zero `ℂ`, zero `exp`) |
| **Concentrated OCRS / Stationary Bipartite S-OCRS** | $\alpha^* = 0.3819660113$ | $\alpha^* = 0.3819660113$ | Tight (MMN 2024/26; ANN 2026) |
| **Finite Simple Forests (Non-Vanishing)** | **`0.382673618862`** (`0.382309325242` 1-coin) | `0.38667` (exact MMN depth-2 tree) | Strictly above $\alpha^*$ (margin `+0.000708`) |
| **All Vanishing Simple Graphs (General & Bipartite)** | **`0.382673618862 - o(1)`** | **`0.3860271555`** (bipartite, $\max_e x_e \le 10^{-12} \to 0^+$) | Strictly above $\alpha^*$ (bracket `[0.38267, 0.38603]`) |
| **Small Zero-Expander Graphs ($K_6$, $K_8$, fixed $\pi$)** | — | **`0.389536778`** ($K_6$), **`0.388559368`** ($K_8$) | Both strictly below `arXiv:2406.14506`'s exact $\beta = 0.3895380533$ |

---

## Scope, Necessary Hypotheses & Trust Boundary

- **Necessary hypotheses for Track 1 (`vanishing-lower-bound/`):** Both **simple graphs (no parallel edges)** and an **oblivious arrival-order adversary** (the exact model of `arXiv:2406.14506`) are strictly load-bearing and tested explicitly in [`adversarial_fuzzer.py`](vanishing-lower-bound/adversarial_fuzzer.py):
  - On bipartite multigraphs with parallel edges, the continuous-staircase selectability drops to `0.352687` ($k=2$) and `0.330779` ($k=200$), below $\alpha^*$.
  - Against an *adaptive* (online) arrival-order adversary on simple forests, probing early edges reveals the global coin and drops selectability to $w_0(1 - c_1)^2 \approx 0.360869 < \alpha^*$ (contrasting with stationary Gibbs S-OCRS in `arXiv:2603.21532`, which achieves $\alpha^*$ against adaptive adversaries).
- **Open boundary:** Whether non-vanishing *cyclic* bipartite graphs admit an oblivious OCRS with selectability $> \alpha^*$ remains open. For general non-bipartite graphs, Track 2 closes 97.0% of the hedging gap to the `0.3786435369` single-layer deterministic triangle-bridge wall; crossing below `0.3786435369` remains open (three structural obstructions are recorded in [`SPEC.md`](SPEC.md) §3.4).
- **Trust boundary:** Every core verifier (`verify_vanishing_lb.py`, `verify_l3_ema_ocrs.py`, `verify_ocrs.py`) uses only Python stdlib + 100-digit outward-rounded `mpmath.iv` interval arithmetic (zero LP solvers, zero floating-point inequalities, zero `assert` statements across all 12 `.py` files under `python3 -O`). Pure-rational checks (`verify_l3_lean4_paley_rat.py`, `test_ema_coupled_finite_toy.py`, and odd-set / Bellman checks) use exact Python `fractions.Fraction`. None of these certificates are formalized in Lean yet (the 71-box and 293-box pure-$\mathbb{Q}$ Paley certificates are structured for zero-`ℂ`/zero-`Real.exp` formalization).

---

## Comparison with Prior Work

| Reference | Prior Result | Relation Here |
|---|---|---|
| **Ezra, Feldman, Gravin & Tang (2020/2022)** | $0.337$-selectable general matching OCRS; identified $\alpha^* = \frac{3-\sqrt{5}}{2}$ independence benchmark. | Strictly separated on both sides: general $\le 0.378673 < \alpha^* < 0.382674 \le$ trees & vanishing. |
| **MacRury, Ma & Grammel (2023/2025)** | $0.349$-selectable general matching OCRS; $0.400$ finite-graph general upper bound. | **Track 2:** Finite simple-graph upper bound improved to $\le 0.3786732$ (and $K_6 \le 0.389537$, $K_8 \le 0.388559$ without expanders). |
| **Ma, MacRury & Nuti (2024/2026)**, [`arXiv:2406.14506`](https://arxiv.org/abs/2406.14506) v3 | $\alpha^*$-selectable vanishing OCRS (`thm:vanishing_to_local`); upper bounds $\gamma = 0.39293532$ (trees/bipartite), $\beta = 0.38953805$ (general vanishing), $\alpha^*$ (concentrated); conjectured optimality of $\alpha^*$. | **Track 1:** Refutes $\alpha^*$ upper bound on simple trees & all vanishing simple graphs ($\ge 0.3826736 > \alpha^*$); finds & bypasses bug in `lem:tree_reduction_game_part_1` (L1572). **Track 3:** Tightens bipartite/vanishing upper bound to $\le 0.3860272$ (cuts 63.0% of $[\alpha^*, \gamma]$ and 46.4% of $[\alpha^*, \beta]$). |
| **Aminian, Niazadeh & Nuti (2026)**, [`arXiv:2603.21532`](https://arxiv.org/abs/2603.21532) | Exact $\alpha^*$-selectable stationary bipartite S-OCRS (surviving adaptive adversaries); quotes MMN's general-matching achievability & tree impossibility conjecture. | **Track 2:** Refutes general achievability of $\alpha^*$ ($\le 0.3786732$). **Track 1:** Separates oblivious global-coin OCRS ($> 0.38267$) from stationary S-OCRS ($= \alpha^*$) and adaptive-order global-coin OCRS ($\le 0.36087$). |

---

## Independent Verification

```bash
# Quick check (~8 s): SHA256SUMS + AST 0-assert check + Track 1 full + pure-Q 71-box Paley + pure-Fraction K6/K8
uv run --with mpmath --with numpy --with scipy python3 -O verify_all.py --fast

# Core verifiers only (zero numpy/scipy; ~1.5 min): all Track 1, Track 2, and Track 3 certificates at 100-dps + pure-Q
uv run --with mpmath python3 -O verify_all.py --verifiers-only

# Full suite (~4.5 min): every verifier + 75 forgery/regression checks + 14 adversarial exact-DP fuzzer suites
uv run --with mpmath --with numpy --with scipy python3 -O verify_all.py
```

---

## Repository Layout & Integrity

- **[`SPEC.md`](SPEC.md)** — Exact problem definitions, complete self-contained proofs for Tracks 1–3, counterexample to `arXiv:2406.14506` v3 L1572, and structural obstruction notes.
- **[`vanishing-lower-bound/`](vanishing-lower-bound/)**, **[`general-upper-bound/`](general-upper-bound/)**, **[`bipartite-upper-bound/`](bipartite-upper-bound/)**, **[`verify_all.py`](verify_all.py)** — Standalone certificates, 100-dps interval and exact-rational verifiers, forgery suites, and top-level runner.
- **[`SHA256SUMS`](SHA256SUMS)**, **[`SHA256SUMS.ots`](SHA256SUMS.ots)**, **[`CITATION.cff`](CITATION.cff)**, **[`LICENSE`](LICENSE)** — Combined SHA-256 integrity manifest (`sha256sum -c SHA256SUMS`), OpenTimestamps Bitcoin-calendar timestamp receipt (`ots verify SHA256SUMS.ots`), citation metadata, and MIT license.
