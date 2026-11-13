# Sequence Card 14: K_5(u) = (u^3 - 31u)/30 -> {-1, 0, 1}

* **Companion Article:** [`Total_Collapse_Primes_18_Families_Zhuk_2026.md`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.md) | [`PDF`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.pdf)
* **Verification Code:** [`../Code/`](../Code/)

## DATA (%S)

```text
7, 17, 23, 37, 313, 8753, 72173, 2952797, 3516973, 8879231, 334677433
```

## NAME (%N)

```text
Primes p >= 7 for which the iteration of K_5(u) = (u^3 - 31*u)/30 on GF(p) collapses every element into {-1, 0, 1}.
```

## OFFSET (%O)

```text
1,1
```

## COMMENTS (%C)

```text
Over Z, 30*(K(u) - 1) = (u + 1)*(u + 5)*(u - 6) and 30*(K(u) + 1) = (u - 1)*(u - 5)*(u + 6), so {-1, 0, 1} is invariant: 0 is fixed and {1, -1} is a 2-cycle.
The 2-periodic points satisfy 30^4*(K(K(u)) - u) = u*(u^2 - 1)*(u^2 - 61)*(u^4 - 31*u^2 + 900), with D_fix = 5^2 + 6^2 = 61 and 2-cycle discriminants 6^3 - 5^3 = 91 and 1 - 30 = -29. Hence every term p satisfies Kronecker(61, p) <= 0 and min(Kronecker(91, p), Kronecker(-29, p)) <= 0.
No other terms below 10^6 (exhaustive full functional-graph traversal of GF(p) for every prime p). The terms above 10^6 were found in a bounded-step orbit scan up to 2*10^9 and individually certified by a complete traversal of GF(p); completeness of the list above 10^6 is not established.
```

## LINKS (%H)

```text
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23420627">Unital Collapse Primes in Arithmetic Dynamics: Total Backward-Orbit Collapse for T(x) = (x^2 - 2)^2 mod p, Galois Stratification, and Fermat-Totient Asymptotics</a>, Zenodo preprint (2026).
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23856837">Total-Collapse Primes of Polynomial Maps over GF(p): Eighteen Families and Two Related Sequences</a>, Zenodo preprint (2026); <a href="https://github.com/Nondual-Observer/total-collapse-primes-oeis">GitHub verification suite</a>.
```

## FORMULA (%F)

```text
p is in the sequence iff #{x in {0..p-1} : K_5^(p)(x) in {-1, 0, 1} (mod p)} = p, where K_5(x) = (x^3 - 31*x)/30 mod p and K_5^(k) denotes the k-th iterate of K_5.
```

## EXAMPLE (%e)

```text
For p = 7: 30 == 2 (mod 7) and 31 == 3 (mod 7), so K_5(u) == (u^3 - 3*u)/2 (mod 7), which collapses GF(7) onto {0, 1, 6} within 2 steps; thus 7 is in the sequence.
For p = 17: every element of GF(17) enters {0, 1, 16} within 4 steps under iteration of K_5(u) mod 17; thus 17 is in the sequence.
```

## PROG (%o)

```text
(Julia)
isprime(n::Int) = n > 1 && all(d -> n % d != 0, 2:isqrt(n))
function isok(p::Int)
    (!isprime(p) || p < 7) && return false
    powermod(61, p >> 1, p) == 1 && return false
    (powermod(91, p >> 1, p) == 1 && powermod(-29, p >> 1, p) == 1) && return false
    inv30, st = invmod(30, p), zeros(UInt8, p)
    st[1] = st[2] = st[p] = 2
    f(c) = mod(mod(c * mod(c*c - 31, p), p) * inv30, p)
    for x in 2:(p >> 1)
        st[x + 1] == 2 && continue
        c = x; while st[c + 1] == 0; st[c + 1] = 1; c = f(c); end
        st[c + 1] == 1 && return false
        c = x; while st[c + 1] == 1; st[c + 1] = st[p - c + 1] = 2; c = f(c); end
    end
    return true
end
println([p for p in 7:10000 if isok(p)])
```

## CROSSREFS (%Y)

```text
Cf. A398745, A401173, A401667 (K_2). See the linked article for the companion sequences for K_n, n = 1, 2, 4, 6, 7, 8.
```

## KEYWORDS (%K)

```text
nonn,hard,more
```

## AUTHOR (%A)

```text
_Igor Zhuk_, Nov 01 2026
```
