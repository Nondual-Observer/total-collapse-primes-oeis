# Total-Collapse Primes of Polynomial Maps over GF(p): Eighteen Families and Two Related Sequences

**Author:** Igor M. Zhuk  
**Date:** October 2026  
**MSC 2020:** 11B83, 37P05, 11T06, 11R37.  
**Preprint DOI:** Zenodo, DOI [10.5281/zenodo.23856837](https://doi.org/10.5281/zenodo.23856837).  
**Companion materials:** OEIS sequence cards in [`../OEIS/`](../OEIS/), exact verification code in [`../Code/`](../Code/). Base family preprint (OEIS A401173): Zenodo, DOI [10.5281/zenodo.23420627](https://doi.org/10.5281/zenodo.23420627).

---

## Abstract

Given a polynomial map $T : \mathbb{F}_p \to \mathbb{F}_p$ and a nonempty target set $A \subset \mathbb{F}_p$, a prime $p$ is called a *total-collapse prime* for $(T, A)$ if the forward orbit of every element $x \in \mathbb{F}_p$ eventually enters $A$. When $A$ is $T$-invariant and consists of periodic points, this means that every directed cycle of the functional graph of $T$ on $\mathbb{F}_p$ lies inside $A$. We study 18 explicit polynomial families:

* the odd cubic ladder $\mathcal{K}_n(u) = \bigl(u^3 - (n^2 + n + 1)u\bigr)/\bigl(n(n + 1)\bigr)$ with target $\{-1, 0, 1\}$ for $n \in \{1, 2, 4, 5, 6, 7, 8\}$;
* squares of monic Chebyshev–Dickson polynomials $V_m(x)^2$ for $m \in \{2, 3, 5, 6, 7, 8, 9\}$ (where $m = 2$ gives OEIS sequence A401173);
* three maps governed by Hilbert class fields ($D = -59$ and $D = -251$) and by the octic tower $(x^2 - 2)^4$;
* the diagonal map $x^2(x - 2)^2$ arising from the $\mathrm{SL}_3$ power-trace map.

For each family we establish explicit necessary conditions for total collapse, derived from factorizations of fixed-point and 2-periodic polynomials, quadratic reciprocity, Stickelberger's parity theorem, the dihedral Galois group $D_7$ of the Hilbert class field of $\mathbb{Q}(\sqrt{-251})$, and Dickson's permutation polynomial criterion. Extending the cubic ladder $\mathcal{K}_n$ to arbitrary $n \ge 1$, we prove a Pell–Pythagorean obstruction: whenever $n^2 + (n + 1)^2 = c^2$ (OEIS A001652), total collapse for $p > n + 1$ can occur only if $c$ is prime and $p = c$. Computational tables are provided for all 20 sequences, distinguishing ranges covered by exhaustive full-graph verification from larger ranges probed by bounded-step orbit scans.

---

## 1. Definitions and General Lemmas

Throughout this paper, $p$ denotes a prime, $\left(\frac{a}{p}\right)$ is the Legendre (or Kronecker) symbol, and $V_m \in \mathbb{Z}[x]$ denotes the monic Chebyshev–Dickson polynomial of the first kind defined by $V_m(t + t^{-1}) = t^m + t^{-m}$, with $V_0(x) = 2$, $V_1(x) = x$, and $V_m(x) = x V_{m-1}(x) - V_{m-2}(x)$. Explicitly:

$$
\begin{aligned}
V_2 &= x^2 - 2, \quad V_3 = x(x^2 - 3), \quad V_5 = x(x^4 - 5x^2 + 5), \quad V_6 = V_3^2 - 2, \\
V_7 &= x(x^6 - 7x^4 + 14x^2 - 7), \quad V_8 = V_2(V_2(V_2)), \quad V_9 = V_3(V_3).
\end{aligned}
$$

**Definition 1.1.** Let $T \in \mathbb{Z}[x]$ (or $T \in \mathbb{Z}[1/B][x]$ for an integer $B$ with $p \nmid B$), and let $A \subset \mathbb{F}_p$ be nonempty. A prime $p$ is a *total-collapse prime* for the pair $(T, A)$ if

$$
\forall x_0 \in \mathbb{F}_p \;\; \exists r \ge 1 : \; T^{(r)}(x_0) \in A.
$$

**Lemma 1.2.** Suppose $T(A) \subseteq A$ and every element of $A$ is $T$-periodic. Then $p$ is a total-collapse prime for $(T, A)$ if and only if every directed cycle of the functional graph $(\mathbb{F}_p, T)$ is contained in $A$. Equivalently, $\#\{x \in \mathbb{F}_p : T^{(p)}(x) \in A\} = p$.

*Proof.* Every forward orbit in the finite set $\mathbb{F}_p$ eventually enters a periodic cycle. If all cycles lie in $A$, every orbit enters $A$. Conversely, suppose there is a cycle $C \not\subset A$. If $C \cap A = \varnothing$, points on $C$ never enter $A$. If $C \cap A \ne \varnothing$, invariance $T(A) \subseteq A$ forces $C \subseteq A$, a contradiction. Since every orbit reaches its cycle in at most $p$ steps and $A$ is invariant, the condition is equivalent to $T^{(p)}(x) \in A$ for all $x \in \mathbb{F}_p$. $\square$

**Lemma 1.3 (commutation of composition factors).** Let $h, g : \mathbb{F}_p \to \mathbb{F}_p$, and set $T = g \circ h$ and $F = h \circ g$. Then $h$ maps the set of periodic points of $T$ bijectively onto the set of periodic points of $F$, preserving cycle lengths and satisfying $h \circ T = F \circ h$. In particular, $C$ is the unique cycle of $T$ if and only if $h(C)$ is the unique cycle of $F$.

*Proof.* From $h \circ T = F \circ h$ and $g \circ F = T \circ g$, the maps $h$ and $g$ send $T$-periodic points to $F$-periodic points and vice versa. Restricted to the periodic sets $\mathrm{Per}(T)$ and $\mathrm{Per}(F)$, the compositions $g \circ h = T$ and $h \circ g = F$ are bijections, so $h|_{\mathrm{Per}(T)}$ is a cycle-preserving bijection conjugating $T$ and $F$. $\square$

**Proposition 1.4 (collapse to a single non-fixed value).** For a singleton target $A = \{a\}$, Definition 1.1 holds if and only if the functional graph of $T$ on $\mathbb{F}_p$ is connected and its unique cycle contains $a$. Allowing $r \ge 0$ in Definition 1.1 gives an equivalent condition.

*Proof.* Taking $x_0 = a$ in Definition 1.1 ($r \ge 1$) shows that the forward orbit of $a$ returns to $a$, so $a$ is periodic and lies on some cycle $C_a$. Any other cycle would be closed under $T$ and disjoint from $C_a$, so its elements could never reach $a$; hence $C_a$ is the unique cycle and the functional graph is connected. The converse is immediate. Allowing $r \ge 0$ does not rescue any cycle disjoint from $\{a\}$, so it implies $a$ lies on every cycle and yields the same condition. $\square$

---

## 2. The Eighteen Polynomial Families

**Definition 2.1.** We study the 18 pairs $(T, A)$ listed below, numbered to match the sequence card files in [`../OEIS/`](../OEIS/).

| Card | Polynomial map $T$ | Target $A$ | Admissible primes $p$ |
| :--- | :--- | :--- | :--- |
| 00 (A401173) | $V_2(x)^2 = (x^2 - 2)^2$ | $\{1\}$ | all primes |
| 01 | $(z + 2)(z^3 + 2z - 1) + z$ | $\{-2\}$ | all primes |
| 02 | $V_3(x)^2$ | $\{0\}$ | all primes |
| 03 | $(x^2 - 2)^4$ | $\{1\}$ | all primes |
| 04 | $V_5(x)^2$ | $\{0, 1\}$ | all primes |
| 05 | $V_8(x)^2$ | $\{1\}$ | all primes |
| 06 | $V_9(x)^2$ | $\{0\}$ | all primes |
| 07 | $V_7(x)^2$ | $\{0, 1\}$ | all primes |
| 08 | $(z + 2) P_{251}(z) + z$ | $\{-2\}$ | all primes |
| 09 | $V_6(x)^2$ | value $4$ (Prop. 1.4) | all primes |
| 10 | $D(x) = x^2(x - 2)^2$ | $\mathrm{Fix}(D)$ | all primes |
| 11 / 12 / 13 / 14 / 15 / 16 / 17 | $\mathcal{K}_n$, $n = 2 / 4 / 1 / 5 / 6 / 7 / 8$ | $\{-1, 0, 1\}$ | $p \ge 5$, $p \nmid n(n + 1)$ |

Here

$$
P_{251}(z) = z^7 + 5z^6 + 6z^5 - 2z^4 - 4z^3 + 2z^2 + 9z - 1, \qquad \mathcal{K}_n(u) = \frac{u^3 - A_n u}{B_n}, \quad B_n = n(n + 1), \quad A_n = B_n + 1.
$$

For every family except Card 09, the target set $A$ is $T$-invariant and consists of periodic points, so Lemma 1.2 applies. Indeed: $V_m(1)^2 = 1$ for $m \equiv \pm 2 \pmod{6}$; $V_m(0) = 0$ for odd $m$; $V_m(1) = 1$ for $m \equiv \pm 1 \pmod{6}$; $\mathcal{K}_n(0) = 0$ and $\mathcal{K}_n(\pm 1) = \mp 1$; and $T(-2) = -2$ for Cards 01 and 08.

**Remark 2.2.** The case $m = 1$ ($x \mapsto x^2$ with target $\{0, 1\}$) yields the Fermat primes (OEIS A019434). For comparison, the degree-8 Chebyshev square $V_4(x)^2 = \bigl((x^2 - 2)^2 - 2\bigr)^2$ with target $\{1\}$ differs from both $(x^2 - 2)^2$ and $(x^2 - 2)^4$ and has collapsing primes $3, 5, 11, 13$ below $2 \cdot 10^4$. The combinatorial motivation for grouping the 18 families studied here is recorded in Appendix A.

---

## 3. Necessary Conditions for Total Collapse

### 3.1. The Odd Cubic Ladder (Cards 11–17)

**Theorem 3.1.** For every integer $n \ge 1$, the following polynomial identities hold in $\mathbb{Z}[n, u]$:

$$
u^3 - A_n u - B_n = (u + 1)(u + n)(u - n - 1), \qquad u^3 - A_n u + B_n = (u - 1)(u - n)(u + n + 1),
$$

$$
B_n^4 \bigl(\mathcal{K}_n(\mathcal{K}_n(u)) - u\bigr) = u(u^2 - 1)\bigl(u^2 - D_{\mathrm{fix}}\bigr)\bigl(u^4 - A_n u^2 + B_n^2\bigr),
$$

where $D_{\mathrm{fix}} = 2B_n + 1 = n^2 + (n + 1)^2$. Viewing the quartic factor as a quadratic in $u^2$, its discriminant factors as

$$
A_n^2 - 4B_n^2 = D_2^+ D_2^-, \qquad D_2^+ = 3B_n + 1 = (n + 1)^3 - n^3, \quad D_2^- = 1 - B_n.
$$

Let $p \ge 5$ be a prime with $p \nmid B_n$. Then:

1. if $\left(\frac{D_{\mathrm{fix}}}{p}\right) = 1$, then $\pm\sqrt{D_{\mathrm{fix}}} \in \mathbb{F}_p$ are fixed points outside $\{-1, 0, 1\}$, so $p$ is not a total-collapse prime;
2. if $\left(\frac{D_2^+}{p}\right) = \left(\frac{D_2^-}{p}\right) = 1$, then the four roots $u = \tfrac{1}{2}\bigl(\pm\sqrt{D_2^+} \pm \sqrt{D_2^-}\bigr) \in \mathbb{F}_p$ form two 2-cycles outside $\{-1, 0, 1\}$, so $p$ is not a total-collapse prime.

*Proof.* The polynomial identities hold by direct expansion (verified symbolically in `Code/algebra_checks.jl`).

For (1), $\sqrt{D_{\mathrm{fix}}} \equiv \pm 1 \pmod{p}$ would imply $2B_n \equiv 0 \pmod{p}$, impossible since $p \ge 5$ and $p \nmid B_n$, while $\sqrt{D_{\mathrm{fix}}} \equiv 0 \pmod{p}$ contradicts $\left(\frac{D_{\mathrm{fix}}}{p}\right) = 1$.

For (2), $u = 0$ is not a root of $u^4 - A_n u^2 + B_n^2$ because $p \nmid B_n$, and $u = \pm 1$ is a root only if $p \mid B_n(B_n - 1)$, i.e., $p \mid D_2^-$. A common root with $u^2 - D_{\mathrm{fix}}$ occurs only if $p \mid B_n(3B_n + 1) = B_n D_2^+$. Both divisibilities are excluded when $\left(\frac{D_2^+}{p}\right) = \left(\frac{D_2^-}{p}\right) = 1$. Since $4u^2 = \bigl(\sqrt{D_2^+} \pm \sqrt{D_2^-}\bigr)^2$ in $\mathbb{F}_p$, all four roots lie in $\mathbb{F}_p$ and have exact period $2$ outside $\{-1, 0, 1\}$. $\square$

**Corollary 3.2 (exclusion of level $n = 3$).** Among $1 \le n \le 8$, the discriminant $D_{\mathrm{fix}} = 5, 13, 25, 41, 61, 85, 113, 145$ is an integer square exclusively at $n = 3$ ($D_{\mathrm{fix}} = 25$). For $n = 3$ and every prime $p > 5$, the residues $\pm 5 \bmod p$ are fixed points outside $\{-1, 0, 1\}$. At $p = 5$ we have $\mathcal{K}_3(u) \equiv 3u^3 + u \pmod{5}$, collapsing $\mathbb{F}_5$ into $\{0, \pm 1\}$ in one step. Thus $p = 5$ is the unique total-collapse prime for $\mathcal{K}_3$, which is why $n = 3$ is omitted from the seven non-degenerate cubic families.

### 3.2. Hilbert Class Field Barriers: Cards 00, 01, 08

**Theorem 3.3.**

1. *(Card 00.)* $T(x) - x = (x - 1)f(x)$, where $f(x) = x^3 + x^2 - 3x - 4$ has $\operatorname{disc}(f) = -83$ and $f(1) = -5$. If $p > 5$, $p \ne 83$, and $\left(\frac{-83}{p}\right) = -1$, then by Stickelberger's theorem $f \bmod p$ has a unique linear factor, yielding a fixed point $x_0 \ne 1$ in $\mathbb{F}_p$, so $p$ is not a total-collapse prime. At $p = 83$, $f(x) \equiv (x - 44)^2(x - 77) \pmod{83}$.
2. *(Card 01.)* $T(z) - z = (z + 2) P _ {59}(z)$, where $P _ {59}(z) = z^3 + 2z - 1$ has $\operatorname{disc}(P _ {59}) = -59$ and $P _ {59}(-2) = -13$. If $p \notin \{13, 59\}$ and $\left(\frac{-59}{p}\right) = -1$, there exists a fixed point $z _ 0 \ne -2$ in $\mathbb{F} _ p$. At $p = 59$, $P _ {59}(z) \equiv (z - 45)^2(z - 28) \pmod{59}$.
3. *(Card 08.)* $T(z) - z = (z + 2) P_{251}(z)$, where $\operatorname{disc}(P_{251}) = -2^{10} \cdot 251^3$, $P_{251}(-2) = -11$, and the splitting field of $P_{251}$ is the Hilbert class field $H$ of $\mathbb{Q}(\sqrt{-251})$ with $\mathrm{Gal}(H/\mathbb{Q}) \cong D_7$ [4, 11]. If $p \notin \{2, 251\}$ and $\left(\frac{-251}{p}\right) = -1$, the Frobenius automorphism at $p$ is a reflection in $D_7$ of cycle type $(1)(2)(2)(2)$, so $P_{251} \bmod p$ has exactly one root in $\mathbb{F}_p$. For $p \ne 11$ this root differs from $-2$, preventing total collapse.

*Proof.* All factorizations and discriminants hold by direct calculation (`Code/algebra_checks.jl`).

By Stickelberger's parity theorem [3], for any prime $p \nmid 2 \operatorname{disc}(f)$ one has $\left(\frac{\operatorname{disc}(f)}{p}\right) = (-1)^{d - r}$, where $d = \deg f$ and $r$ is the number of irreducible factors of $f \bmod p$. For $d = 3$ and Legendre symbol $-1$, we obtain $r = 2$, forcing the factorization pattern $1 + 2$ and hence a unique root in $\mathbb{F}_p$.

In (3), any prime inert in $\mathbb{Q}(\sqrt{-251})$ has Frobenius lying in the non-trivial coset $D_7 \setminus C_7$, which consists entirely of reflections. Every reflection of a regular heptagon fixes one vertex and swaps the remaining six in three pairs; by Dedekind's theorem (for $p$ not dividing the index $2^5$), $P_{251} \bmod p$ factors as $(1)(2)(2)(2)$, while $p = 2$ is excluded directly ($T(0) \equiv 0 \not\equiv -2 \pmod{2}$, $T(1) \equiv 1 \pmod{2}$). Across all primes $p < 3000$ ($p \ne 251$), the factorization pattern of $P_{251} \bmod p$ matches the $D_7$ splitting law without exception. $\square$

**Corollary 3.4 (representation by binary quadratic forms).** Both cubics $f$ and $P_{59}$ have fundamental prime discriminants $-83$ and $-59$, whose imaginary quadratic fields have class number $3$. Hence the Galois closure of the cubic field is unramified over $\mathbb{Q}(\sqrt{D})$ and coincides with its Hilbert class field [11]. A prime $p$ splitting in $\mathbb{Q}(\sqrt{D})$ keeps the cubic irreducible modulo $p$ if and only if $p$ is represented by the non-principal binary quadratic form of discriminant $D$. Consequently:

* for Card 00 and $p > 5$, $x = 1$ is the unique fixed point in $\mathbb{F}_p$ if and only if $p = 3x^2 + xy + 7y^2$;
* for Card 01 and $p > 13$, $z = -2$ is the unique fixed point in $\mathbb{F}_p$ if and only if $p = 3x^2 + xy + 5y^2$;
* for Card 08 and $p > 11$, $z = -2$ is the unique fixed point in $\mathbb{F}_p$ if and only if $p$ is represented by one of the non-principal forms $3x^2 + xy + 21y^2$, $7x^2 + xy + 9y^2$, or $5x^2 + 3xy + 13y^2$.

Every listed term satisfies these form conditions (explicit integer representations are recorded in Cards 00, 01, and 08). Furthermore, the substitution $z = x^2 - 3$ maps the roots of $f(x) = x^3 + x^2 - 3x - 4$ to the roots of $z^3 + 2z^2 + 2z - 1$.

### 3.3. Even Chebyshev–Dickson Harmonics: Preimage and Permutation Barriers

**Theorem 3.5.**

1. *(Cards 02, 06.)* Let $p > 3$ with $\left(\frac{3}{p}\right) = -1$. Then $V_3(x) \equiv 0 \pmod{p}$ holds only at $x = 0$, so the preimage basin of $0$ under $V_3^2$ and $V_9^2$ reduces to $\{0\}$, excluding $p$. Thus every odd prime in Cards 02 and 06 satisfies $p \equiv \pm 1 \pmod{12}$.
2. *(Cards 00, 03, 05.)* Let $p \equiv 7 \pmod{12}$, so that $\left(\frac{-1}{p}\right) = \left(\frac{3}{p}\right) = -1$. Then $(u^2 - 2)^2 \equiv 1 \pmod{p}$ implies $u^2 \in \{1, 3\}$, hence $u \equiv \pm 1 \pmod{p}$, while $-1$ has no square or fourth-power root in $\mathbb{F}_p$. By induction along the composition $V_2(V_2(\cdot))$, the backward basin of $1$ equals $\{1, -1\}$ of size $2$, excluding $p$.

**Theorem 3.6 (Dickson permutation barrier, Cards 04, 07).** Let $m \in \{5, 7\}$, $p > m$, and $\gcd(m, p^2 - 1) = 1$. Then $V_m$ is a permutation polynomial on $\mathbb{F}_p$ [2, 5] satisfying $V_m(0) = 0$ and $V_m(\pm 1) = \pm 1$. Consequently $T^{-1}(0) = \{0\}$ and $T^{-1}(1) = \{\pm 1\}$ for $T = V_m^2$. If moreover $p \equiv 3 \pmod{4}$, then $-1$ has no preimage under $T$, so the total backward basin of $\{0, 1\}$ is $\{0, 1, -1\}$ of size $3 < p$, excluding $p$.

**Proposition 3.7 (degree-halving reductions).** By Lemma 1.3:

* Card 00 ($T = s \circ V_2$, $s(y) = y^2$) is equivalent to total collapse of $u \mapsto u^4 - 2$ to $-1$;
* Card 03 ($T = g \circ h$, $h(x) = x^2 - 2$, $g(u) = u^4$) is equivalent to total collapse of $F(u) = u^8 - 2$ to $-1$;
* Card 05 is equivalent to total collapse of $F(u) = V_8(u^2)$ to $-1$.

In particular, for Card 03 and $p > 3$ it is necessary that $(u^8 - 2 - u)/(u + 1) = u^7 - u^6 + u^5 - u^4 + u^3 - u^2 + u - 2$ (discriminant $-26522311$) have no roots in $\mathbb{F}_p$.

**Proposition 3.8 (Cards 09, 10).**

1. For $T_6(x) = V_6(x)^2$, one has $T_6(0) = T_6(\pm 1) = T_6(\pm 2) = 4$. The value $4$ is a fixed point if and only if $p \mid V_6(4)^2 - 4 = 2^6 \cdot 3^3 \cdot 5^2 \cdot 13^2$, and total collapse is characterized by Proposition 1.4.
2. For $D(x) = x^2(x - 2)^2$, one has $D(x) - x = x(x - 1)(x^2 - 3x + 1)$ and $D(2 - x) = D(x)$. For odd $p$, the fixed-point set $\mathrm{Fix}(D)$ equals $\{0, 1\}$ when $\left(\frac{5}{p}\right) = -1$ and $\bigl\{0, 1, \tfrac{3 \pm \sqrt{5}}{2}\bigr\}$ when $\left(\frac{5}{p}\right) \in \{0, 1\}$. Total collapse means absence of periodic cycles of length $\ge 2$. Algebraically, $D$ is the restriction to the diagonal $x = y$ of the coordinatewise square $(x, y) \mapsto \bigl((x^2 - 2y)^2, (y^2 - 2x)^2\bigr)$ of the $\mathrm{SL}_3$ power-trace map $(\operatorname{tr} A, \operatorname{tr} A^{-1}) \mapsto (\operatorname{tr} A^2, \operatorname{tr} A^{-2})$.

---

## 4. Two Related Integer Sequences (Cards 18 and 19)

### 4.1. Least Collapsing Prime Along the Cubic Ladder (Card 18)

**Definition 4.1.** For each integer $n \ge 1$, let $a(n)$ be the least prime $p > n + 1$ ($p \ge 5$) that is a total-collapse prime for $\mathcal{K}_n$ with target $\{-1, 0, 1\}$, or $a(n) = 0$ if no such prime exists.

The bound $p > n + 1$ avoids modular aliasing: $\mathcal{K} _ n \bmod p$ depends only on $n \bmod p$ (for instance, $\mathcal{K} _ {119} \equiv \mathcal{K} _ 2 \pmod{13}$, while $13$ collapses $\mathcal{K} _ 2$).

**Theorem 4.2 (Pell–Pythagorean obstruction).** Suppose $n^2 + (n + 1)^2 = c^2$, so that $n \in \text{A001652} = \{3, 20, 119, 696, 4059, 23660, \dots\}$ and $c \in \text{A001653}$. Then for every prime $p > n + 1$:

* the residues $\pm c \bmod p$ are fixed points of $\mathcal{K}_n$;
* $\pm c \in \{-1, 0, 1\} \pmod{p}$ holds if and only if $c$ is prime and $p = c$, in which case $\mathcal{K}_n(u) \equiv u(1 - 2u^2) \pmod{c}$.

Consequently, $a(n) = 0$ whenever $c$ is composite, and $a(n) \in \{0, c\}$ whenever $c$ is prime.

*Proof.* By Theorem 3.1, $D_{\mathrm{fix}} = c^2$ gives the integer fixed points $\pm c$. Every prime factor of $c^2 - 1 = 2n(n + 1)$ is at most $n + 1$, so $p \nmid (c \mp 1)$. Moreover, $c < \sqrt{2}(n + 1) < 2(n + 1)$. Since $c$ is odd, if $c$ is composite then every prime factor of $c$ is at most $c/3 < n + 1 < p$, so $p \nmid c$. If $c$ is prime, $p \mid c$ holds iff $p = c$. Finally, $2B_n = c^2 - 1 \equiv -1 \pmod{c}$ gives $B_n^{-1} \equiv -2 \pmod{c}$ and $A_n = B_n + 1 \equiv \tfrac{1}{2} \pmod{c}$, whence $\mathcal{K}_n(u) \equiv -2(u^3 - \tfrac{1}{2}u) \equiv u(1 - 2u^2) \pmod{c}$. $\square$

**Examples.**

* $n = 3$, $c = 5$: $a(3) = 5$.
* $n = 20$, $c = 29$: $\mathcal{K}_{20} \bmod 29$ has the 4-cycle $22 \to 12 \to 7 \to 17 \to 22$, hence $a(20) = 0$.
* $n = 119$ ($c = 169 = 13^2$) and $n = 696$ ($c = 985 = 5 \cdot 197$): $c$ is composite, so $a(119) = a(696) = 0$ unconditionally.
* $n = 4059$ ($c = 5741$) and $n = 23660$ ($c = 33461$): direct computation gives $a(4059) = a(23660) = 0$ (e.g., via the 4-cycle $1762 \to 3531 \to 3979 \to 2210 \to 1762 \pmod{5741}$).
* For the 7th, 8th, and 9th Pell legs, $c = 195025 = 5^2 \cdot 29 \cdot 269$, $1136689 = 137 \cdot 8297$, and $6625109 = 37 \cdot 179057$ are all composite, giving $a(n) = 0$.

**Conjecture 4.3.** $a(n) = 0$ for all $n \in \text{A001652}$ with $n > 3$.

For $n \notin \text{A001652}$, existence of at least one collapsing prime $p > n + 1$ remains open. Exact values $a(1), \dots, a(55)$ are given in Card 18 (`Code/collapse_verifier.c`, mode `card18`).

### 4.2. Sums of Binary One-Zero-Block Primes (Card 19)

**Definition 4.4.** For $n \ge 2$, consider the $n$-bit integers

$$
B(n, k, s) = (2^n - 1) - 2^s(2^k - 1), \qquad 1 \le k \le n - 2, \quad 1 \le s \le n - 1 - k,
$$

together with the Mersenne number $2^n - 1$ itself. Equivalently, these are the $n$-bit odd integers whose binary expansion has a single circular run of 1s and at most one internal block of 0s (the odd $n$-bit terms of OEIS A100724). Define $a(n)$ as the sum of all primes of this form.

Exact values $a(2), \dots, a(25)$ and the corresponding prime counts (1, 2, 2, 5, 3, 7, 9, 11, 7, 14, 9, 16, 12, 13, 22, 27, 13, 32, 27, …) are given in Card 19. At $n = 7$ the primes are $67, 71, 79, 97, 103, 113, 127$ (all five single-zero candidates $127 - 2^s$ with $k = 1$ are composite), giving $a(7) = 657$. At $n = 10$ the primes are $769, 911, 967, 991, 1009, 1019, 1021$, giving $a(10) = 6687$. At $n = 16$ there are $22$ such primes summing to $a(16) = 1\,203\,372$, whose two largest summands are the twin primes $65\,519$ and $65\,521 = 2^{16} - 15$ straddling $65\,520 = \operatorname{lcm}\bigl(\{1, \dots, 16\} \setminus \{11\}\bigr)$.

**Proposition 4.5 (divisibility filter modulo 3).** For every candidate $B(n, k, s) = (2^n - 1) - 2^s(2^k - 1)$:

1. if $n \ge 4$ is even, then $3 \mid (2^n - 1)$, and whenever $k$ is even one also has $3 \mid (2^k - 1)$, so $3 \mid B(n, k, s)$ with $B(n, k, s) \ge 2^{n-1} + 1 > 3$ composite; in particular, for even $n \ge 4$ every prime summand of $a(n)$ has odd zero-block length $k$;
2. if $n \ge 5$ is odd, then $2^n - 1 \equiv 1 \pmod{3}$, and whenever $k$ is odd and $s$ is even one has $2^s(2^k - 1) \equiv 1 \pmod{3}$, so again $3 \mid B(n, k, s)$ and $B(n, k, s) > 3$ is composite.

*Proof.* Substituting $2 \equiv -1 \pmod{3}$ gives $2^m - 1 \equiv (-1)^m - 1 \pmod{3}$, which is $0 \pmod{3}$ for even $m$ and $1 \pmod{3}$ for odd $m$, while the smallest candidate for $n \ge 4$ is $B(n, n - 2, 1) = 2^{n-1} + 1 \ge 9 > 3$. $\square$

---

## 5. Computational Verification

### 5.1. Exact Verification and Extended Probe Regimes

For each candidate prime $p$, exact membership is decided in two stages (`Code/collapse_verifier.c`):

1. *Rejection* occurs only upon exhibiting a *certified* extraneous periodic cycle: Brent's cycle-detection algorithm is run from deterministic starting residues with **no step cap** (guaranteed to terminate on the finite set $\mathbb{F}_p$). Any detected cycle violating the target condition certifies that $p$ is not a total-collapse prime.
2. *Acceptance* occurs only after an **exhaustive functional-graph traversal** of all $p$ elements of $\mathbb{F}_p$ using a 3-color state array (`unvisited` / `on current path` / `resolved`), verifying that every cycle satisfies the target condition.

Only the algebraically proven filters of Section 3 (Theorems 3.1, 3.3, 3.5, 3.6) are applied prior to graph traversal.

In Table 5.1 we distinguish two search bounds for each family:

* $P_{\mathrm{exact}}$ — the bound up to which every prime $p \le P_{\mathrm{exact}}$ was tested by the exact unbounded algorithm above ($3 \cdot 10^{10}$ for Card 00; $10^8$ for Cards 01–03; $10^6$ for Cards 04–17), guaranteeing complete term lists up to $P_{\mathrm{exact}}$;
* $P _ {\mathrm{scan}}$ — a larger range ($2 \cdot 10^9$ or $5 \cdot 10^8$ for Cards 01–17) screened by bounded-step orbit probes ($\approx 6\sqrt{p}\text{--}8\sqrt{p}$ steps) to locate large collapsing primes. Every prime discovered in $(P _ {\mathrm{exact}}, P _ {\mathrm{scan}}]$ with $p < 2 \cdot 10^9$ was subsequently certified by an exact full-field traversal of $\mathbb{F} _ p$, while completeness of the lists above $P _ {\mathrm{exact}}$ remains open because an unclosed long transient does not prove the existence of an extraneous cycle.

### 5.2. Summary of Known Terms

**Table 5.1.** Known total-collapse primes across the 18 polynomial families, with exhaustive exact bounds $P_{\mathrm{exact}}$ and bounded-step probe ranges $P_{\mathrm{scan}}$.

| Card | Known collapsing primes | $P_{\mathrm{exact}}$ | $P_{\mathrm{scan}}$ |
| :--- | :--- | :--- | :--- |
| 00 | 3, 5, 37, 59, 61, 109, 359, 5507, 2117573, 503734381, 20900146669 | $3 \cdot 10^{10}$ | $3 \cdot 10^{10}$ |
| 01 | 3, 13, 137, 761, 2719, 150083 | $10^8$ | $2 \cdot 10^9$ |
| 02 | 2, 13, 23, 83, 17929, 18121 | $10^8$ | $2 \cdot 10^9$ |
| 03 | 3, 5, 13, 23, 41, 22381, 45503, 3852347, 1801981681 | $10^8$ | $2 \cdot 10^9$ |
| 04 | 2, 3, 5, 31, 59, 61, 89, 2939, 4261, 11839, 14137859, 599999971 | $10^6$ | $2 \cdot 10^9$ |
| 05 | 3, 5, 13, 17, 23, 47, 71, 89, 193, 313, 863, 1103, 1777, 107183, 84814153, 342037921 | $10^6$ | $2 \cdot 10^9$ |
| 06 | 2, 11, 73, 1297, 6733, 7393, 658871, 1295447 | $10^6$ | $2 \cdot 10^9$ |
| 07 | 2, 3, 5, 43, 421, 1429, 6257, 49279, 1580503, 24107329 | $10^6$ | $5 \cdot 10^8$ |
| 08 | 3, 7, 17, 233, 1373, 1531, 5861, 27427, 52937, 125287 | $10^6$ | $5 \cdot 10^8$ |
| 09 | 2, 3, 5, 7, 11, 13, 17, 37, 47, 79, 131, 167, 241, 257, 129401, 390209, 4926997, 6719077, 23734297 | $10^6$ | $5 \cdot 10^8$ |
| 10 | 2, 3, 5, 7, 11, 23, 29, 31, 59, 181, 521, 1559, 2069, 5261, 6871, 9551, 19207, 1161371, 1574681, 1849721, 85268251 | $10^6$ | $5 \cdot 10^8$ |
| 11 | 5, 7, 13, 19, 31, 47, 59, 73, 151, 331, 401, 1943069, 83584483, 158450797, 1554906533 | $10^6$ | $2 \cdot 10^9$ |
| 12 | 7, 11, 17, 19, 29, 41, 53, 79, 191, 42979, 277829, 35718901 | $10^6$ | $2 \cdot 10^9$ |
| 13 | 5, 7, 13, 23, 83, 97, 383, 14239193, 82411573 | $10^6$ | $2 \cdot 10^9$ |
| 14 | 7, 17, 23, 37, 313, 8753, 72173, 2952797, 3516973, 8879231, 334677433 | $10^6$ | $2 \cdot 10^9$ |
| 15 | 5, 11, 17, 31, 41, 109, 127, 199, 719, 2731, 26317, 307067, 617579, 4348559, 4838117, 30282823, 42166321, 62371087 | $10^6$ | $2 \cdot 10^9$ |
| 16 | 5, 19, 23, 151, 1237, 11941, 761357, 6756377, 1105256599 | $10^6$ | $2 \cdot 10^9$ |
| 17 | 5, 7, 19, 23, 31, 41, 53, 61, 277, 593, 809, 911, 1489, 2053, 4943, 6301, 17393, 229139, 324119, 1262783, 2055551, 228229363 | $10^6$ | $2 \cdot 10^9$ |

Across all 18 families there are 224 prime entries representing 143 distinct primes. Every listed entry with $p < 2 \cdot 10^9$ (all 223 such entries, including $p = 1\,801\,981\,681$ in Card 03 and $p = 1\,554\,906\,533$ in Card 11) is certified by a complete functional-graph traversal of $\mathbb{F}_p$.

For Card 09, the lengths of the unique cycle passing through $4$ along the 19 listed primes are $1, 1, 1, 2, 2, 1, 2, 2, 2, 12, 11, 5, 5, 9, 399, 475, 1125, 1065, 2612$.

### 5.3. Near-Collapse Phenomena

Alongside total collapse ($k = 0$ extraneous cycles), one encounters primes for which all of $\mathbb{F}_p$ outside the tiny basins of one or two short extraneous cycles collapses to the target (reaching target-basin density $> 99.9998\%$). For example, in Card 00 at $p = 150662749$, there are just two extraneous cycles of lengths $3$ and $5$ whose combined basins occupy $234 + 48 = 282$ residues out of $150\,662\,749$ (`./collapse_verifier analyze 0 150662749`).

### 5.4. Finiteness vs. Infinitude

**Question 5.2.** Is any of the 18 prime sequences finite, or are all 18 infinite? Finiteness is proven unconditionally only for the degenerate cubic level $n = 3$ (Corollary 3.2) and for composite-hypotenuse Pell levels (Theorem 4.2).

A naive random-mapping heuristic points toward infinitude: among all functions on a $p$-element set with a distinguished fixed point, the fraction having no other cycles equals $p^{p-2}/p^{p-1} = 1/p$ by Cayley's tree formula [8], and $\sum_p 1/p$ diverges (giving the doubly logarithmic heuristic $O(\log \log X)$ for the prime count up to $X$). However, low-degree polynomials with structured dynatomic factorizations are far from uniform random mappings, so this heuristic is not a proof [6, 9, 10].

---

## 6. Reproducibility

The [`../Code/`](../Code/) directory provides three standalone verification programs:

* `collapse_verifier.c` — multi-threaded exact C verifier for all 18 families and Card 18;
* `reference_implementation.jl` — dependency-free Julia reference implementation of all 20 sequences (`p < 20000` self-test);
* `algebra_checks.jl` — dependency-free exact symbolic and arithmetic verification in Julia of every polynomial identity, discriminant, quadratic-form representation, $D_7$ splitting pattern, and combinatorial count in Sections 3–4 and Appendices A–B.

Build and run instructions are given in [`../Code/README.md`](../Code/README.md).

---

## Appendix A. Projective Combinatorics of PG(3, 2) and Indexing of the 18 Families

This appendix records the projective combinatorics underlying the $7 + 7 + 3 + 1$ indexing of the 18 families.

**Proposition A.1 (book of three Fano planes in PG(3, 2)).** Identify the 15 points of $\mathrm{PG}(3, 2)$ with the nonzero 4-bit integers $1, \dots, 15$. Bit-reversal $\mathrm{Rev}_4$ fixes the central line $L_0 = \{2, 4, 6\}$ pointwise. Exactly three Fano planes in $\mathrm{PG}(3, 2)$ contain $L_0$:

$$
L_0 \cup \{1, 3, 5, 7\}, \qquad L_0 \cup \{8, 10, 12, 14\}, \qquad L_0 \cup \{9, 11, 13, 15\},
$$

defined respectively by the hyperplanes $b_3 = 0$, $b_0 = 0$, and $b_3 = b_0$. The second plane is the image of the first under $\mathrm{Rev}_4$.

On the third sheet $\{9, 11, 13, 15\}$, one has $9 + 4b_2 + 2b_1 \equiv b_2 - b_1 \pmod{3}$. The two primes $11 = 1011_2$ and $13 = 1101_2$ encode the unique pair of reciprocal irreducible cubics over $\mathbb{F}_2$, namely $x^3 + x + 1$ and $x^3 + x^2 + 1$, satisfying

$$
(x^3 + x + 1)(x^3 + x^2 + 1) = \Phi_7(x) + 2x^3, \qquad 11 \otimes 13 = 127, \qquad 11 \cdot 13 = 127 + 16,
$$

where $\otimes$ denotes carryless multiplication over $\mathbb{F}_2$.

**Proposition A.2 (cyclic and projective orbits).** The cyclic group $C_7$ acts freely by bit-rotations on $\mathbb{F}_2^7 \setminus \{0, 127\}$, partitioning the $126$ non-constant vectors into $126 / 7 = 18$ orbits distributed across Hamming weights $1, \dots, 6$ as $1 + 3 + 5 + 5 + 3 + 1$. Under the full collineation group $\mathrm{PSL}_2(7)$ of order $168$, non-empty proper subsets of $\mathrm{PG}(2, 2)$ split into one orbit each in weights $1, 2, 5, 6$ (of sizes $7, 21, 21, 7$) and two orbits each in weights $3, 4$ (of sizes $7 + 28$, lines vs. non-collinear triples).

**Remark A.3 (combinatorial indexing).** We index the 18 polynomial families via the partition $18 = 7 + 7 + 3 + 1$:

* 7 odd cubic maps $\mathcal{K}_n$ — indexed by the 7 points of $\mathrm{PG}(2, 2)$;
* 7 Chebyshev–Dickson squares $V_m^2$ — indexed by the 7 lines of $\mathrm{PG}(2, 2)$;
* 3 maps of Cards 01, 03, 08 — indexed by the 3 Fano planes through $L_0$ in Proposition A.1;
* Card 10 — indexed by the ambient plane.

Note that squaring $h(u) = u^2$ semiconjugates $\mathcal{K}_n$ to $G_n(x) = x(x - A_n)^2 / B_n^2$, which is distinct from the Chebyshev–Dickson squares $V_m^2$; thus the $7 + 7 + 3 + 1$ correspondence serves purely as a combinatorial index.

---

## Appendix B. Auxiliary Arithmetic Identities

* **B.1.** On $\{1, \dots, N\}$, the count of integers with $\Omega(k) \le 1$ (namely $1$ together with the primes $\le N$) equals the count of composites if and only if $N \in \{10, 12, 14\}$. Indeed, for odd $N \ge 15$, the balance $2(1 + \pi(N)) - N$ is non-increasing by [1] and equals $-1$ at $N = 15$.
* **B.2.** Across the first octave $1 \le n \le 8$, the cubic-ladder coefficients satisfy $\sum_{n=1}^8 B_n = 2\binom{10}{3} = 240$ and $\sum_{n=1}^8 A_n = 248$.
* **B.3.** Joining two octahedral graphs $K_{2,2,2}$ via a central apex adjacent to all 12 vertices produces a 13-vertex, 36-edge graph with adjacency spectrum $\{6, 4, 0^{(6)}, (-2)^{(5)}\}$.

---

## References

1. J. B. Rosser, L. Schoenfeld. Approximate formulas for some functions of prime numbers. *Illinois J. Math.* 6 (1962), 64–94.
2. L. E. Dickson. The analytic representation of substitutions on a power of a prime number of letters with a discussion of the linear group. *Ann. of Math.* 11 (1896/97), 65–120.
3. L. Stickelberger. Über eine neue Eigenschaft der Diskriminanten algebraischer Zahlkörper. *Verh. I. Int. Math.-Kongr.*, Zürich (1897), 182–193.
4. B. C. Berndt, H. H. Chan. Ramanujan and the modular $j$-invariant. *Canad. Math. Bull.* 42 (1999), 427–440; E. Konstantinou, A. Kontogeorgis. Computing polynomials of the Ramanujan $t_n$ class invariants. *Canad. Math. Bull.* 52 (2009), 583–597.
5. R. Lidl, H. Niederreiter. *Finite Fields*, 2nd ed. Cambridge Univ. Press, 1997 (Ch. 7: Dickson polynomials).
6. J. H. Silverman. *The Arithmetic of Dynamical Systems*. GTM 241, Springer, 2007.
7. OEIS Foundation Inc. *The On-Line Encyclopedia of Integer Sequences*, https://oeis.org (A000668, A001652, A001653, A001844, A002315, A003215, A014234, A019434, A086383, A100724, A401173).
8. P. Flajolet, A. M. Odlyzko. Random mapping statistics. *EUROCRYPT '89*, LNCS 434 (1990), 329–354.
9. T. Vasiga, J. Shallit. On the iteration of certain quadratic maps over GF(p). *Discrete Math.* 277 (2004), 219–240.
10. S. V. Konyagin, F. Luca, B. Mans, L. Mathieson, M. Sha, I. E. Shparlinski. Functional graphs of polynomials over finite fields. *J. Combin. Theory Ser. B* 116 (2016), 87–122.
11. D. A. Cox. *Primes of the Form* $x^2 + ny^2$, 2nd ed. Wiley, 2013.
