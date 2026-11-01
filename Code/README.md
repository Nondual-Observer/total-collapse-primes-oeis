# Verification Suite (`Code/`)

| File | Purpose | Requirements |
| :--- | :--- | :--- |
| `reference_implementation.jl` | Dependency-free Julia reference implementation of all 20 sequences (full functional-graph traversal, no step cap). Self-test checks Cards 00–17 on all primes $p < 20000$, Card 18 ($n = 1 \dots 55$), and Card 19 ($n = 2 \dots 25$). | Julia ($\ge 1.9$, standard library only) |
| `algebra_checks.jl` | Dependency-free exact polynomial and arithmetic verification of all identities in the companion article (Sections 3–4, Appendices A–B): factorizations, Bareiss Sylvester discriminants and resultants, quadratic-form representations, $D_7$ Frobenius splitting types, Pell–Pythagorean barrier, and $\mathrm{PSL}_2(7)$ orbit counts. | Julia ($\ge 1.9$, `Base` + `LinearAlgebra` only) |
| `collapse_verifier.c` | Multi-threaded exact C verifier for Families 00–17 and Card 18. Rejects a prime *only* upon certifying an extraneous cycle via unbounded Brent cycle detection, and accepts *only* after a complete functional-graph traversal of $\mathrm{GF}(p)$. | C compiler (`clang` or `gcc`), `pthreads` |

## Build & Run

```bash
julia reference_implementation.jl            # Expected output: ALL OK (~2 s)
julia algebra_checks.jl                      # Every output line starts with PASS (~0.5 s)

clang -O3 -o collapse_verifier collapse_verifier.c -lpthread
./collapse_verifier scan 11 2 100000000 8    # Exhaustive exact scan of Card 11 up to 10^8
./collapse_verifier check 11 1554906533      # Certify a single prime (~1.6 GB RAM)
./collapse_verifier analyze 0 150662749      # List all cycles and basin sizes (near-collapse)
./collapse_verifier card18 1 55 20000        # Compute all 55 terms of Card 18
```

The family index `<fam>` in `collapse_verifier` matches the card number `0..17`. Full-field traversal uses $p$ bytes of RAM.

## Verified Search Ranges

* **Exhaustive exact verification (no step cap):** $p \le 3 \cdot 10^{10}$ for Card 00; $p \le 10^8$ for Cards 01–03; $p \le 10^6$ for Cards 04–17.
* **Individual full-field certification:** every prime listed in Cards 00–17 with $p < 2 \cdot 10^9$ (including `1801981681` in Card 03 and `1554906533` in Card 11) is certified by a complete traversal of $\mathrm{GF}(p)$ via `./collapse_verifier check`.
* **Cards 18 and 19:** all 55 terms of Card 18 (including $a(20) = a(4059) = a(23660) = 0$) and all 24 terms of Card 19 are computed directly from the definitions.
* Beyond the exhaustive bounds above, Cards 01–17 were probed up to $2 \cdot 10^9$ or $5 \cdot 10^8$ with bounded-step orbit scans; see Section 5 of the companion article.
