# Sequence Card 06: V_9(x)^2 -> {0}

* **Companion Article:** [`Total_Collapse_Primes_18_Families_Zhuk_2026.md`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.md) | [`PDF`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.pdf)
* **Verification Code:** [`../Code/`](../Code/)

## DATA (%S)

```text
2, 11, 73, 1297, 6733, 7393, 658871, 1295447
```

## NAME (%N)

```text
Primes p such that the iteration of T_9(x) = V_9(x)^2 (where V_9(x) = V_3(V_3(x)) and V_3(x) = x*(x^2 - 3) is the monic Chebyshev polynomial of degree 3) collapses every element of GF(p) to 0.
```

## OFFSET (%O)

```text
1,1
```

## COMMENTS (%C)

```text
Here V_9(x) = V_3(V_3(x)) = x*(x^2 - 3)*(x^6 - 6*x^4 + 9*x^2 - 3) = 2*T_9(x/2) and T(x) = V_9(x)^2 has degree 18; V_9(0) = 0 and V_9(1) = -2.
For odd p > 3, T(x) == 0 (mod p) requires V_3(x) in {0, +-sqrt(3)}, and 0 has a nonzero preimage iff Kronecker(3, p) = 1; hence every odd term satisfies p == +-1 (mod 12). The prime 2 is a term, while 3 is not.
No other terms below 10^6 (exhaustive full functional-graph traversal of GF(p) for every prime p). The terms above 10^6 were found in a bounded-step orbit scan up to 2*10^9 and individually certified by a complete traversal of GF(p); completeness of the list above 10^6 is not established.
```

## LINKS (%H)

```text
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23420627">Unital Collapse Primes in Arithmetic Dynamics: Total Backward-Orbit Collapse for T(x) = (x^2 - 2)^2 mod p, Galois Stratification, and Fermat-Totient Asymptotics</a>, Zenodo preprint (2026).
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23856837">Total-Collapse Primes of Polynomial Maps over GF(p): Eighteen Families and Two Related Sequences</a>, Zenodo preprint (2026); <a href="https://github.com/Nondual-Observer/total-collapse-primes-oeis">GitHub verification suite</a>.
```

## FORMULA (%F)

```text
p is in the sequence iff #{x in {0..p-1} : T_9^(p)(x) == 0 (mod p)} = p, where T_9(x) = V_3(V_3(x))^2 mod p, V_3(x) = x*(x^2 - 3), and T_9^(k) denotes the k-th iterate of T_9.
```

## EXAMPLE (%e)

```text
For p = 2: T_9(0) = 0, T_9(1) = 4 == 0 (mod 2), so all 2 elements of GF(2) reach 0 in 1 step; thus 2 is in the sequence.
For p = 11: every element of GF(11) reaches 0 within 4 steps under iteration of T_9(x) mod 11; thus 11 is in the sequence.
```

## PROG (%o)

```text
(Julia)
isprime(n::Int) = n > 1 && all(d -> n % d != 0, 2:isqrt(n))
function isok(p::Int)
    (!isprime(p) || (p > 2 && !(p % 12 in (1, 11)))) && return false
    f(c) = (v3 = mod(c * mod(c*c - 3, p), p); v9 = mod(v3 * mod(v3*v3 - 3, p), p); mod(v9*v9, p))
    st = zeros(UInt8, p); st[1] = 2
    for x in 1:(p >> 1)
        st[x + 1] == 2 && continue
        c = x; while st[c + 1] == 0; st[c + 1] = 1; c = f(c); end
        st[c + 1] == 1 && return false
        c = x; while st[c + 1] == 1; st[c + 1] = st[p - c + 1] = 2; c = f(c); end
    end
    return true
end
println([p for p in 2:8000 if isok(p)])
```

## CROSSREFS (%Y)

```text
Cf. A401173. See the linked article for the companion sequence for V_3(x)^2.
```

## KEYWORDS (%K)

```text
nonn,hard,more
```

## AUTHOR (%A)

```text
_Igor Zhuk_, Nov 01 2026
```
