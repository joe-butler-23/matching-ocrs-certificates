# Formal Specification & Self-Contained Mathematical Proofs

This document states the exact definitions, proof chains, and trust boundaries for all three tracks in this repository.

---

## 1. Model & Definitions

Let $G = (V, E)$ be a finite **simple** graph (no self-loops, no parallel edges) and let $\mathbf{x} = (x_e)_{e \in E} \in [0, 1]^E$ belong to the matching polytope $\mathcal{P}_{\mathrm{match}}(G)$ (convex hull of indicator vectors of matchings of $G$; in particular $\sum_{e \in \partial(v)} x_e \le 1$ for every $v \in V$, and $\sum_{e \in E(U)} x_e \le \lfloor |U|/2 \rfloor$ for every odd $U \subseteq V$). Each edge $e \in E$ is independently active ($X_e = 1$) with probability $x_e$.

- **Adversarial-Order Oblivious Matching OCRS:** Before execution, an oblivious adversary fixes $(G, \mathbf{x})$ and a deterministic arrival order $\pi$ on $E$ (knowing the algorithm's description, not its internal randomness or $(X_e)$). Edges arrive in order $\pi$; upon observing $X_e$, if $X_e = 1$ and neither endpoint of $e$ is matched, the algorithm immediately and irrevocably decides whether to add $e$ to the matching $\mathcal{M}$. An OCRS is **$c$-selectable** if $\Pr[e \in \mathcal{M} \mid X_e = 1] \ge c$ for every $e \in E$ (with $x_e > 0$).
- **Strongly-Online OCRS (`arXiv:2406.14506` v3, L679–683):** The algorithm initially knows only $V$ (and any upfront internal randomness $\omega$) and learns $(x_e, X_e)$ only as edges arrive in oblivious order $\pi$.
- **Vanishing Regime (`arXiv:2406.14506` v3, L340):** An OCRS is $c$-selectable in the vanishing regime if there exists $f(\varepsilon) \to 0$ as $\varepsilon \to 0^+$ such that for all simple $(G, \mathbf{x})$ with $\max_{e \in E} x_e \le \varepsilon$, $\Pr[e \in \mathcal{M} \mid X_e = 1] \ge c - f(\varepsilon)$.
- **Golden-Ratio Independence Benchmark:** $\alpha^* := \frac{3 - \sqrt{5}}{2} = 0.38196601125010515\dots$, the smaller root of $\alpha = (1 - \alpha)^2$.

---

## 2. Track 1: Strongly-Online Global-Coin Lower Bound (`vanishing-lower-bound/`)

### 2.1 Continuous-Staircase Algorithm
Fix any rational $c_1 \in (\alpha^*, 1 - W(1))$ (certified at $c_1 = \frac{49}{125} = 0.392$), and define

$$\tau := \frac{1 - 2c_1}{c_1(1 - c_1)} = \frac{3375}{3724}, \qquad B := 1 - \tau = \frac{349}{3724}, \qquad A(b) := \frac{2}{2b + 1 + \sqrt{4b + 1}} \quad (b \in [0, 1]),$$

so $(1 - c_1\tau)(1 - c_1) = c_1$ and $A(b) = (1 - b A(b))^2$ identically on $[0, 1]$. Since $\frac{1}{A(b)} = b + \frac{1}{2} + \frac{1}{2}\sqrt{1 + 4b}$,

$$I(B) := \int_0^B \frac{db}{A(b)} = \frac{6B^2 + 6B - 1 + (1 + 4B)\sqrt{1 + 4B}}{12}, \qquad w_0 := \frac{1}{1 + c_1(1 - c_1)I(B)},$$

$$c_{\mathrm{exact}} := w_0 c_1 = 0.3826736188630302229\ldots > 0.382673618862 > \alpha^*.$$

Before any edges arrive, draw one global random choice $\omega \in \{0\} \cup [\tau, 1]$: with probability $w_0$ pick Branch $0$; with probability $1 - w_0$ draw $T \in [\tau, 1]$ from density $p(t) = \frac{w_0 c_1(1 - c_1)}{(1 - w_0) A(1 - t)}$. Maintain deterministic prefix degrees $x_v(e) := \sum_{f \in \partial(v), f < e} x_f$:

- **Branch $0$:** Set $\gamma_e^{(0)} := \min\bigl(c_1,\; (1 - c_1 x_u(e))(1 - c_1 x_v(e))\bigr)$ and $D_v^{(0)}(e) := \sum_{f \in \partial(v), f < e} \gamma_f^{(0)} x_f \le c_1 x_v(e) \le c_1$. Accept active feasible $e = (u, v)$ independently with probability $a_e^{(0)} := \frac{\gamma_e^{(0)}}{(1 - D_u^{(0)}(e))(1 - D_v^{(0)}(e))} \in [0, 1]$.
- **Branch $T$ ($T \in [\tau, 1]$):** Let $E(T) := \{e = (u, v) : x_u(e) \ge T \text{ and } x_v(e) \ge T\}$, $b := 1 - T$, $A := A(b)$. Ignore edges outside $E(T)$; for active feasible $e \in E(T)$, accept independently with probability $a_e^{(T)} := \frac{A}{(1 - A x_u^{(T)}(e))(1 - A x_v^{(T)}(e))}$, where $x_v^{(T)}(e) := \sum_{f \in \partial(v) \cap E(T), f < e} x_f$.

### 2.2 Feasibility, Pointwise Guarantee & Theorem 1 (Simple Forests)
1. **Non-Vanishing Residual-Degree Lemma on $E(T)$:** If $\partial(v) \cap E(T) \neq \emptyset$, let $e^*$ be its earliest edge. Then $x_v(e^*) \ge T$, and every $e \in \partial(v) \cap E(T)$ satisfies $e \ge e^*$, so $\sum_{e \in \partial(v) \cap E(T)} x_e \le 1 - x_v(e^*) \le 1 - T = b$. Hence $D_v^{(T)}(e) = A x_v^{(T)}(e) \le b A(b)$ and $(1 - D_u^{(T)}(e))(1 - D_v^{(T)}(e)) \ge (1 - b A(b))^2 = A(b)$, so $a_e^{(T)} \in [0, 1]$ with $D_{\max}^{(T)} \le B A(B) \approx 0.07942 < c_1 = 0.392$.
2. **Pointwise Expectation Bound:** For any edge $e = (u, v)$, let $s = \min(x_u(e), x_v(e))$ and $r = \max(x_u(e), x_v(e))$.
   - If $s < \tau$: $(1 - c_1 s)(1 - c_1 r) > (1 - c_1\tau)(1 - c_1) = c_1$, so $\gamma_e^{(0)} = c_1$ and $\mathbb{E}_\omega[\gamma_e^{(\omega)}] = w_0 c_1 = c_{\mathrm{exact}}$.
   - If $s \in [\tau, 1]$: both $(1 - c_1 s)(1 - c_1 r) \ge (1 - c_1 s)(1 - c_1)$ and $c_1 = (1 - c_1\tau)(1 - c_1) \ge (1 - c_1 s)(1 - c_1)$, so $\gamma_e^{(0)} \ge (1 - c_1)(1 - c_1 s)$. Branch $T$ contributes $A(1 - T)$ whenever $T \in [\tau, s]$, and $A(1 - t)$ cancels inside the integral:
     $$\mathbb{E}_\omega[\gamma_e^{(\omega)}] \;\ge\; w_0(1 - c_1)(1 - c_1 s) + \int_\tau^s \frac{w_0 c_1(1 - c_1)}{A(1 - t)}\,A(1 - t)\,dt \;=\; w_0(1 - c_1)(1 - c_1\tau) \;=\; c_{\mathrm{exact}}.$$
3. **Theorem 1 (Simple Finite Forests, Arbitrary Non-Vanishing $\mathbf{x}$):** For fixed $\omega$, every $a_f^{(\omega)}$ is deterministic given $(G, \mathbf{x}, \pi)$. On a simple forest, deleting $e = (u, v)$ splits the graph into vertex-disjoint subtrees, so $\{u \text{ free at } e^-\}$ and $\{v \text{ free at } e^-\}$ are independent. Strong induction on arrival order gives $\Pr[f \in \mathcal{M}^{(\omega)}] = \gamma_f^{(\omega)} x_f$ and $\Pr[u \text{ free at } e^-] = 1 - D_u^{(\omega)}(e)$ for all $f, e$, hence $\Pr[e \in \mathcal{M} \mid X_e = 1] = \mathbb{E}_\omega[\gamma_e^{(\omega)}] \ge c_{\mathrm{exact}}$. The same proof holds for the 2-branch and 65-branch discrete rational mixtures in `cert_vanishing_lb.json`.

### 2.3 Theorem 2 (6-Line Disagreement Coupling) & Theorem 3 (All Simple Vanishing Graphs)
Let $e = (u, v)$ be $\ell$-locally-tree-like in a simple graph $G$, meaning the **induced subgraph** $T = G[N^\ell_G(e)]$ on vertices at graph distance $\le \ell$ from $\{u, v\}$ is a tree. Fix $\omega$, let $D_{\max} = c_1 = 49/125$, and $\beta_0 := \frac{1}{1 - c_1}\ln\frac{1}{1 - c_1} \approx 0.818389 < 1$. Couple the real run $\mathcal{M}_G$ with a tree run $\mathcal{M}_{\tilde{T}}$ on $\tilde{T} = T \cup \{\tilde{e}_s : s \in D^\ell(e)\}$, sharing $(X_f, A_f)$ and $a_f$ (imported from $G$) on all $f \in E(T)$, while replacing outside edges at each boundary leaf $s \in D^\ell(e)$ with a private dummy edge $\tilde{e}_s$ arriving before $s$'s unique $T$-edge $e_s$ that occupies $s$ independently with probability $D_s(e_s)$:

1. **Exact tree marginal:** Every internal vertex $v \in N^{\ell-1}(e)$ satisfies $\partial_G(v) = \partial_T(v)$, and each leaf $s \in D^\ell(e)$ is free at $e_s^-$ in $\mathcal{M}_{\tilde{T}}$ with probability $1 - D_s(e_s)$. Induction on $\tilde{T}$ yields $\Pr[e \in \mathcal{M}_{\tilde{T}} \mid X_e = 1] = \gamma_e$.
2. **Disagreement tracing:** If $w_0 \in \{u, v\}$ disagrees between $\mathcal{M}_G$ and $\mathcal{M}_{\tilde{T}}$ at $e^-$, let $e_1 = (w_0, w_1) < e$ be the earliest incident edge whose outcome differs. Since $(X_{e_1}, A_{e_1})$ is shared and $w_0$ agreed before $e_1$, necessarily $X_{e_1} A_{e_1} = 1$ and $w_1$ already disagreed before $e_1$.
3. **No backtracking to $D^\ell(e)$:** Inducting yields a strictly decreasing-time edge sequence $e > e_1 > \dots > e_\ell$ with $e_j = (w_{j-1}, w_j)$ and $X_{e_j} A_{e_j} = 1$. Since $T$ is a simple tree and $e_{j+1} < e_j \implies e_{j+1} \neq e_j$, depth increases by $+1$ at every step until $w_\ell \in D^\ell(e)$.
4. **Path independence:** Distinct edges in $E(T)$ have independent $(X_{e_j}, A_{e_j})$, so $\Pr[\bigwedge_{j=1}^\ell (X_{e_j} A_{e_j} = 1)] = \prod_{j=1}^\ell x_{e_j} a_{e_j}$.
5. **Integral branching contraction:** At any $w_{j-1} \in D^{j-1}(e)$ with parent edge $e_{j-1}$ ($e_0 = e$), the child edges $f_1 < \dots < f_k$ arriving before $e_{j-1}$ are all edges incident to $w_{j-1}$ before $f_i$, so $D_{w_{j-1}}(f_i) = \sum_{m < i} d_m$ where $d_m = \gamma_{f_m} x_{f_m}$ and $\sum_{m=1}^k d_m \le D_{\max}$. Since $t \mapsto \frac{1}{1-t}$ is increasing on $[0, 1)$,
   $$\sum_{i=1}^k x_{f_i} a_{f_i} \;\le\; \frac{1}{1 - D_{\max}}\sum_{i=1}^k \frac{d_i}{1 - \sum_{m < i} d_m} \;\le\; \frac{1}{1 - D_{\max}}\int_0^{D_{\max}}\frac{dt}{1 - t} \;=\; \beta_0 < 1.$$
6. **Union bound & Theorem 3:** Backward induction gives $\Pr[w_0 \text{ disagrees at } e^-] \le \beta_0^\ell$. Since $A_e$ is independent of prior edges and $a_e \le \gamma_e / (1 - c_1)^2$, averaging over $\omega$ proves **Theorem 2**:
   $$\Pr[e \in \mathcal{M}_G \mid X_e = 1] \;\ge\; c_{\mathrm{exact}}\!\left(1 - \frac{2\beta_0^\ell}{(1 - c_1)^2}\right).$$
   Because $\mathcal{M}$ is strongly online (drawing $\omega$ upfront independent of the subsampled graph $G'$ and edge states), plugging Theorem 2 into Ma–MacRury–Nuti's black-box `thm:vanishing_to_local` (`arXiv:2406.14506` v3, L688–697, L878–925) proves **Theorem 3**: $(1 - o(1))c_{\mathrm{exact}} \ge 0.382673618862 - o(1)$ on all simple vanishing graphs. $\square$

### 2.4 Counterexample to `arXiv:2406.14506` v3 `lem:tree_reduction_game_part_1` (L1572)
On vertices $\{0,\dots,6\}$ with the 6-cycle `0-2-4-6-5-3-0` plus pendant `e = (0, 1)` arriving in order `[(3,5), (5,6), (6,4), (2,4), (0,2), (0,3), (0,1)]` and strictly feasible $\mathbf{x} = [0.5, 0.5, 0.5, 0.5, 0.495, 0.495, 0.01]$ ($\max_v \deg_{\mathbf{x}}(v) = 1$), $G[N^2_G(e)]$ is cycle-free ($\ell = 2$) and $h_2 = (5, 6) \in \partial T \setminus W(T)$. Exact $2^7$-state evaluation of MMN's `alg:tree_aom` gives `0.381881917636` with $h_2$ vs. `0.381966011250` without $h_2$ (gap `-8.4093614528e-5 < 0`). Theorem 2 avoids `lem:tree_reduction_game_part_1` completely.

---

## 3. Track 2: General Non-Bipartite Upper Bounds (`general-upper-bound/`)

### 3.1 Top-Down Integer Matching Pyramid $G^{\mathrm{EMA}}_{q, P, K_0, N_L, L, \beta, \varepsilon_1, \varepsilon_2}$
Fix $T_2 \ge 1$, prime $q \equiv 1 \pmod{T_2 + 1}$ ($d_0 = (q-1)/(T_2 + 1)$, multiplicative cosets $\mathbb{F}_q^\times = C_0 \sqcup \dots \sqcup C_{T_2}$), integers $P \ge 3, K_0 \ge 2, N_L \ge 2, L \ge 2$, and rationals $\beta = 1/\text{beta\_den} \in (0, 1/2]$, $\varepsilon_1, \varepsilon_2 \in (0, 1/4]$ ($\varepsilon_1 + 2\varepsilon_2 < 1$). Level $k \in \{1,\dots,L\}$ has $N_k := N_L P^{L-k}$ blocks of $q$ triangle-bridge triples $(a^{(k)}_{r,i}, b^{(k)}_{r,i}, w^{(k)}_{r,i})_{i \in \mathbb{F}_q}$ (export pool $D^{(k)} = \{a^{(k)}_{r,i}, b^{(k)}_{r,i}\}$); Level $0$ has $N_0 := K_0 N_L P^L$ blocks paired with private leaves. Set EMA weights $\alpha_{k,0} = (1-\beta)^{k-1}$ and $\alpha_{k,j} = \beta(1-\beta)^{k-1-j}$ ($1 \le j < k$), so $\sum_{j=0}^{k-1} \alpha_{k,j} = 1$.

Sampling matchings top-down ($k = L, \dots, 1$), each still-free block of Level $k$ independently picks:
- $M_1^{\mathrm{bulk}}$ / $M_2^{\mathrm{bulk}}$ (prob. $\frac{1 - \varepsilon_1 - 2\varepsilon_2}{2}$ each): match one side to $w^{(k)}_{r,i}$; pick $j < k$ w.p. $\alpha_{k,j}$, pop a uniformly random still-free block $r'$ of Level $j$, pick side $\sigma \in \{0,1\}$ and bulk shift $s \in C_0 \sqcup \dots \sqcup C_{T_2-1}$, and match the other side to $d^{(j)}_{r',\sigma,i+s}$.
- $M_1^{\mathrm{atom}}$ / $M_2^{\mathrm{atom}}$ (prob. $\varepsilon_2$ each): identical with atom shift $s \in C_{T_2}$.
- $M_3^{\mathrm{dir}}$ (prob. $\varepsilon_1$): match $(a^{(k)}_{r,i}, b^{(k)}_{r,i})$.
- At $k = 0$, match all still-free Level-$0$ blocks to private leaves.

Because $\sum_{k=j+1}^L N_k < N_j / (P - 1) \le N_j / 2 < N_j$ ($j \ge 1$) and $< N_0 / (2 K_0) < N_0$ ($j = 0$), free blocks never exhaust. Every sample is a valid integer matching of a simple graph ($\mathbf{x} \in \mathcal{P}_{\mathrm{match}}(G^{\mathrm{EMA}})$). By $S_{N_j}$-invariance, every Level-$k$ block survives to Level $k$ with uniform probability $\gamma_k \in (1 - \frac{\beta}{P-1}, 1]$.

### 3.2 Frozen-Capacity Relaxation & Pointwise Cubic Potential Identity
Arrive Level $0$, then Levels $k = 1, \dots, L$ in order (within Level $k$: bridge $(a,w), (b,w) \to$ bulk waves $0 \dots T_2 - 1 \to$ atom wave $T_2 \to$ direct $(a,b)$). Relax capacity on surviving vertices $F_j \subseteq D^{(j)}$ to $\infty$ for all future levels $k > j$ (keeping capacity $1$ within Level $j$). Every valid OCRS on $G^{\mathrm{EMA}}$ remains valid in the relaxed instance $\widetilde{G}$ with identical selectabilities, and conditioned on $\mathcal{F}_{k-1}$ all $N_k q$ triples of Level $k$ decouple. Thinning edge probabilities by $\gamma_k \le 1$ can only decrease the optimal Bellman value with nonnegative survival bonus $+\lambda$.

Let $S_j = |F_j|/(2 N_j q)$ and $\tilde{S}_k = (1-\beta)\tilde{S}_{k-1} + \beta S_k = \sum_{j=0}^{k-1} \alpha_{k,j} S_j$. For any $\kappa \ge 0$, setting $\lambda = \kappa \tilde{S}_{k-1}^2 \ge 0$ and $\Delta = \tilde{S}_k - \tilde{S}_{k-1} = \beta(S_k - \tilde{S}_{k-1}) \in [-\beta, \beta]$ satisfies the exact pointwise identity

$$-\kappa \tilde{S}_{k-1}^2 (S_k - \tilde{S}_{k-1}) \;=\; \frac{\kappa}{3\beta}\bigl(\tilde{S}_{k-1}^3 - \tilde{S}_k^3\bigr) + \frac{\kappa}{\beta}\Delta^2\!\left(\tfrac{2}{3}\tilde{S}_{k-1} + \tfrac{1}{3}\tilde{S}_k\right) \;\le\; \frac{\kappa}{3\beta}\bigl(\tilde{S}_{k-1}^3 - \tilde{S}_k^3\bigr) + \kappa\beta.$$

Combining multi-level Gauss-period mixing error $\mathrm{err}_{\mathrm{gauss}} = (W_{\mathrm{p2}} + a_A + a_B)\frac{(T_2+1)\sqrt{q}}{2(q-1)}$, bulk binomial error $\mathrm{err}_{\mathrm{binom}} = W_{\mathrm{p2}}\frac{T_2+1}{8 T_2 N_L (q-1)}$ (zero binomial error on the atom wave via Boole's union bound), and telescoping over $k = 1, \dots, L$ certifies

$$c^*(G^{\mathrm{EMA}}) \;\le\; \frac{\sup_{S \in [0, 1]}\bigl[V_{\mathrm{Bellman}}(S;\, w,\, +\kappa S^2) - \kappa S^3\bigr] + \mathrm{err}_{\mathrm{gauss}} + \mathrm{err}_{\mathrm{binom}} + \kappa\beta + \frac{\kappa}{3\beta L}}{\bigl(1 - \frac{\beta}{P-1}\bigr)\bigl(1 - \varepsilon_1\frac{w_a+w_b}{2} - (\varepsilon_1+2\varepsilon_2)W_{\mathrm{p2}}\bigr)}.$$

For $T_2 = 1$ and $q = 100{,}000{,}007 \equiv 3 \pmod 4$, the Paley matrix identity $A^T A = B^T B = \frac{q+1}{4}I_q + \frac{q-3}{4}J_q$, a complete 9-prime Pratt primality tree down to $2$, and the 5-term Bonferroni polynomial $P_5(x) \ge 1 - e^{-x}$ eliminate complex numbers, square roots, and exponentials (`verify_l3_lean4_paley_rat.py`: 71 pure-$\mathbb{Q}$ boxes for $< 381/1000$, 150 boxes for $< 0.380$).

### 3.3 Pure-Finite $K_6$ & $K_8$ Duals and Deterministic-Wall Obstructions
- **Zero-expander small graphs (`verify_pure_finite_k6_k8.py`):** Exact convex decompositions $K_6 = K_{3,3}(\frac{1-\varepsilon}{3}) \sqcup 2K_3(\frac{\varepsilon}{3})$ and $K_8 = K_{4,4}(\frac{1-\varepsilon}{4}) \sqcup 2K_4(\frac{\varepsilon}{3})$ at $\varepsilon = 10^{-9}$ satisfy all Blossom odd-set inequalities (26 on $K_6$, 120 on $K_8$); exact $2^n$-state pure-`Fraction` backward Bellman duals certify fixed-order $K_6 \le 0.389536778$ and $K_8 \le 0.388559368$.
- **Why crossing below `0.3786435369` remains open:** (1) any vanishing-edge cascade of triangle bridges admits the uncorrelated constant-$S = 1-c$ OCRS solving $\frac{c}{2(1-c)}(1 - ce^{-(1-c)}) = 1 - 2c$ (`c_det = 0.3786435369`); (2) via ghost-matching edge thinning, every $(H, \mathbf{x})$ admitting a $c$-selectable OCRS admits one with $\Pr[v \text{ free}] = 1 - c\deg_{\mathbf{x}}(v)$ for every $v$, ruling out forced over-matching gadgets; (3) replacing vanishing Poisson waves by non-vanishing correlated-pair edges incurs an $O(d)$ discretization loss ($Sd - (1 - e^{-Sd})$) exceeding maximal theoretical anti-correlation across all tested $(r, d)$.

---

## 4. Track 3: Vanishing Bipartite Upper Bound (`bipartite-upper-bound/`)

In `bipartite-upper-bound/certs/cert_bipartite_t40_vanishing.json`, a 2-colored bipartite simple graph $G^{\mathrm{bip}}_{q, K, \varepsilon}$ ($T_2 = 40$, $q = 10^{13}+313$, $K = 10^6$, $\varepsilon = 10^{-12}$) is decomposed via Birkhoff–von Neumann into two integer matching families (Phase-1/Phase-2 bulk coset waves and Phase-2 $t = 1^-$ atom coset waves plus Phase-0 private leaves on pool $D$), with exact $\max_{e \in E} x_e = 10^{-12}$. Fixing $q$ and the dual integers and setting $\varepsilon(\delta) = \delta$, $K(\delta) = \max(10^6, \lceil 1/\delta \rceil)$ for every $\delta \in (0, 10^{-12}]$ satisfies $\max_{e \in E(G_\delta)} x_e(\delta) \le \delta \to 0^+$, non-increasing binomial error, and non-decreasing denominator $Z(q, K(\delta), \delta) \ge Z(q, 10^6, 10^{-12}) > 0$ (verified in exact `Fraction` inside `verify_ocrs.py`), certifying $c^*(G_\delta, x_\delta) \le 0.3860271554305342$ uniformly as $\delta \to 0^+$.
