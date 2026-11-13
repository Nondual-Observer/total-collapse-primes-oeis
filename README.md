# Total-Collapse Primes of Polynomial Maps over GF(p): OEIS Sequence Cards, Companion Article, and Exact Verification Code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23856837.svg)](https://doi.org/10.5281/zenodo.23856837)

Publication package accompanying **20 integer sequences** by Igor M. Zhuk: **18 total-collapse prime sequences** for explicit polynomial maps over $\mathbb{F}_p$ (`Cards 00–17`, starting with OEIS [A401173](https://oeis.org/A401173)) and **2 related integer sequences** (`Cards 18–19`).

```text
.
├── README.md
├── CITATION.cff
├── Article/
│   ├── Total_Collapse_Primes_18_Families_Zhuk_2026.md   # Companion mathematical paper (GitHub MathJax)
│   └── Total_Collapse_Primes_18_Families_Zhuk_2026.pdf  # Typeset PDF preprint
├── OEIS/                                                # 20 OEIS sequence cards (%S %N %O %C %H %F %e %o %Y %K %A)
└── Code/                                                # Standalone Julia and C verifiers
```

## Sequence Directory

| Card | Polynomial Map and Target over $\mathbb{F}_p$ | Card File | OEIS Reference |
| :--- | :--- | :--- | :--- |
| **00** | $(x^2 - 2)^2 \to \{1\}$ | [`OEIS/CARD_00_A401173_x2m2_squared_to_1.md`](OEIS/CARD_00_A401173_x2m2_squared_to_1.md) | [A401173](https://oeis.org/A401173) |
| **01** | $(z + 2)(z^3 + 2z - 1) + z \to \{-2\}$ | [`OEIS/CARD_01_tower_D59_to_m2.md`](OEIS/CARD_01_tower_D59_to_m2.md) | [A401959](https://oeis.org/A401959) |
| **02** | $V_3(x)^2 \to \{0\}$ | [`OEIS/CARD_02_V3_squared_to_0.md`](OEIS/CARD_02_V3_squared_to_0.md) | New sequence |
| **03** | $(x^2 - 2)^4 \to \{1\}$ | [`OEIS/CARD_03_x2m2_fourth_power_to_1.md`](OEIS/CARD_03_x2m2_fourth_power_to_1.md) | New sequence |
| **04** | $V_5(x)^2 \to \{0, 1\}$ | [`OEIS/CARD_04_V5_squared_to_0_1.md`](OEIS/CARD_04_V5_squared_to_0_1.md) | New sequence |
| **05** | $V_8(x)^2 \to \{1\}$ | [`OEIS/CARD_05_V8_squared_to_1.md`](OEIS/CARD_05_V8_squared_to_1.md) | New sequence |
| **06** | $V_9(x)^2 \to \{0\}$ | [`OEIS/CARD_06_V9_squared_to_0.md`](OEIS/CARD_06_V9_squared_to_0.md) | New sequence |
| **07** | $V_7(x)^2 \to \{0, 1\}$ | [`OEIS/CARD_07_V7_squared_to_0_1.md`](OEIS/CARD_07_V7_squared_to_0_1.md) | New sequence |
| **08** | $(z + 2) P_{251}(z) + z \to \{-2\}$ | [`OEIS/CARD_08_tower_D251_to_m2.md`](OEIS/CARD_08_tower_D251_to_m2.md) | [A402158](https://oeis.org/A402158) |
| **09** | $V_6(x)^2 \to$ value $4$ (unique cycle through $4$) | [`OEIS/CARD_09_V6_squared_to_4.md`](OEIS/CARD_09_V6_squared_to_4.md) | New sequence |
| **10** | $x^2(x - 2)^2 \to$ fixed points | [`OEIS/CARD_10_diagonal_x2_xm2_squared_fixed_points.md`](OEIS/CARD_10_diagonal_x2_xm2_squared_fixed_points.md) | New sequence |
| **11** | $\mathcal{K}_2(u) = (u^3 - 7u)/6 \to \{-1, 0, 1\}$ | [`OEIS/CARD_11_cubic_K2.md`](OEIS/CARD_11_cubic_K2.md) | [A401667](https://oeis.org/A401667) |
| **12–17** | $\mathcal{K}_n(u) = \bigl(u^3 - (n^2 + n + 1)u\bigr)/\bigl(n(n + 1)\bigr) \to \{-1, 0, 1\}$, $n = 4, 1, 5, 6, 7, 8$ | [`OEIS/CARD_12_cubic_K4.md`](OEIS/CARD_12_cubic_K4.md) … [`CARD_17_cubic_K8.md`](OEIS/CARD_17_cubic_K8.md) | New sequences |
| **18** | Least collapsing prime $p > n + 1$ for $\mathcal{K}_n$ (or $0$ if none) | [`OEIS/CARD_18_least_prime_cubic_Kn.md`](OEIS/CARD_18_least_prime_cubic_Kn.md) | [A398745](https://oeis.org/A398745) |
| **19** | Sum of $n$-bit primes with a single circular run of 1s and at most one internal block of 0s | [`OEIS/CARD_19_binary_one_zero_block_prime_sums.md`](OEIS/CARD_19_binary_one_zero_block_prime_sums.md) | [A401666](https://oeis.org/A401666) |

## Quick Verification

```bash
cd Code
julia reference_implementation.jl              # Self-test of all 20 cards (expects: ALL OK)
julia algebra_checks.jl                        # Symbolic & arithmetic proof checks (all PASS)

clang -O3 -o collapse_verifier collapse_verifier.c -lpthread
./collapse_verifier scan 11 2 1000000 8        # Exhaustive exact scan for Card 11 up to 10^6
```

## Citation & Preprint Links

* **Zenodo Archive & Preprint DOI:** [10.5281/zenodo.23856837](https://doi.org/10.5281/zenodo.23856837)
* **Companion Article (Markdown & PDF):** [`Article/Total_Collapse_Primes_18_Families_Zhuk_2026.md`](Article/Total_Collapse_Primes_18_Families_Zhuk_2026.md) · [`Article/Total_Collapse_Primes_18_Families_Zhuk_2026.pdf`](Article/Total_Collapse_Primes_18_Families_Zhuk_2026.pdf)
* **Base Family Preprint (`Card 00 / A401173`):** Igor M. Zhuk, *Unital Collapse Primes in Arithmetic Dynamics: Total Backward-Orbit Collapse for* $T(x) = (x^2 - 2)^2 \bmod p$, *Galois Stratification, and Fermat-Totient Asymptotics*, Zenodo Preprint (2026), DOI: [10.5281/zenodo.23420627](https://doi.org/10.5281/zenodo.23420627).
