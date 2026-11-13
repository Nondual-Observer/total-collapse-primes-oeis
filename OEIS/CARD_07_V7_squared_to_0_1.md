# Sequence Card 07: V_7(x)^2 -> {0, 1}

* **Companion Article:** [`Total_Collapse_Primes_18_Families_Zhuk_2026.md`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.md) | [`PDF`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.pdf)
* **Verification Code:** [`../Code/`](../Code/)

## DATA (%S)

```text
2, 3, 5, 43, 421, 1429, 6257, 49279, 1580503, 24107329
```

## NAME (%N)

```text
Primes p such that every trajectory of T(x) = x^2*(x^6 - 7*x^4 + 14*x^2 - 7)^2 in GF(p) eventually reaches one of the two fixed points 0 or 1.
```

## OFFSET (%O)

```text
1,1
```

## COMMENTS (%C)

```text
Here V_7(x) = x*(x^6 - 7*x^4 + 14*x^2 - 7) = 2*T_7(x/2); 0 and 1 are fixed points of T(x) = V_7(x)^2 for every p.
By Dickson's theorem V_7 permutes GF(p) when gcd(7, p^2 - 1) = 1; if moreover p > 7 and p == 3 (mod 4), the basin of {0, 1} is exactly {0, 1, -1}, so p is not a term. All known terms p > 5 satisfy p == +-1 (mod 7).
T(4) = 10084^2 = 101687056.
No other terms below 10^6 (exhaustive full functional-graph traversal of GF(p) for every prime p). The terms above 10^6 were found in a bounded-step orbit scan up to 5*10^8 and individually certified by a complete traversal of GF(p); completeness of the list above 10^6 is not established.
```

## LINKS (%H)

```text
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23420627">Unital Collapse Primes in Arithmetic Dynamics: Total Backward-Orbit Collapse for T(x) = (x^2 - 2)^2 mod p, Galois Stratification, and Fermat-Totient Asymptotics</a>, Zenodo preprint (2026).
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23856837">Total-Collapse Primes of Polynomial Maps over GF(p): Eighteen Families and Two Related Sequences</a>, Zenodo preprint (2026); <a href="https://github.com/Nondual-Observer/total-collapse-primes-oeis">GitHub verification suite</a>.
```

## FORMULA (%F)

```text
p is in the sequence iff #{x in {0..p-1} : T^(p)(x) in {0, 1} (mod p)} = p, where T(x) = x^2*(x^6 - 7*x^4 + 14*x^2 - 7)^2 mod p and T^(k) denotes the k-th iterate of T.
```

## EXAMPLE (%e)

```text
For p = 3: T(0) = 0, T(1) = T(2) = 1 (mod 3), so every element of GF(3) reaches 0 or 1 in 1 step; thus 3 is in the sequence.
For p = 43: p - 1 = 42 = 2*3*7, and every element of GF(43) reaches either 0 or 1 within 5 steps under iteration of x^2*(x^6 - 7*x^4 + 14*x^2 - 7)^2 mod 43; thus 43 is in the sequence.
```

## PROG (%o)

```text
(Julia)
isprime(n::Int) = n > 1 && all(d -> n % d != 0, 2:isqrt(n))
function isok(p::Int)
    !isprime(p) && return false
    f(c) = (x2 = mod(c*c, p); t = mod((x2 - 7)*x2 + 14, p); v7 = mod(c * (mod(t*x2, p) - 7), p); mod(v7*v7, p))
    st = zeros(UInt8, p); st[1] = st[2] = st[p] = 2
    for x in 2:(p >> 1)
        st[x + 1] == 2 && continue
        c = x; while st[c + 1] == 0; st[c + 1] = 1; c = f(c); end
        st[c + 1] == 1 && return false
        c = x; while st[c + 1] == 1; st[c + 1] = st[p - c + 1] = 2; c = f(c); end
    end
    return true
end
println([p for p in 2:50000 if isok(p)])
```

## CROSSREFS (%Y)

```text
Cf. A019434, A401173. See the linked article for the companion sequence for V_5(x)^2.
```

## KEYWORDS (%K)

```text
nonn,hard,more
```

## AUTHOR (%A)

```text
_Igor Zhuk_, Nov 01 2026
```
