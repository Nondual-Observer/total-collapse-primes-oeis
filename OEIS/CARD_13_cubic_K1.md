# Sequence Card 13: K_1(u) = (u^3 - 3u)/2 -> {-1, 0, 1}

* **Companion Article:** [`Total_Collapse_Primes_18_Families_Zhuk_2026.md`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.md) | [`PDF`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.pdf)
* **Verification Code:** [`../Code/`](../Code/)

## DATA (%S)

```text
5, 7, 13, 23, 83, 97, 383, 14239193, 82411573
```

## NAME (%N)

```text
Primes p >= 5 for which the iteration of K_1(u) = (u^3 - 3*u)/2 mod p on GF(p) collapses every element into the triad {-1, 0, 1} mod p.
```

## OFFSET (%O)

```text
1,1
```

## COMMENTS (%C)

```text
Over Z, 2*(K(u) - 1) = (u + 1)^2*(u - 2) and 2*(K(u) + 1) = (u - 1)^2*(u + 2), so {-1, 0, 1} is invariant: 0 is fixed and {1, -1} is a 2-cycle.
The 2-periodic points satisfy 2^4*(K(K(u)) - u) = u*(u^2 - 1)*(u^2 - 5)*(u^4 - 3*u^2 + 4), with D_fix = 1^2 + 2^2 = 5 and 2-cycle discriminants 2^3 - 1^3 = 7 and 1 - 2 = -1. Hence every term p satisfies Kronecker(5, p) <= 0 and min(Kronecker(7, p), Kronecker(-1, p)) <= 0.
Here 2*K(u) = V_3(u) = u^3 - 3*u is the monic Chebyshev polynomial of degree 3.
No other terms below 10^6 (exhaustive full functional-graph traversal of GF(p) for every prime p). The terms above 10^6 were found in a bounded-step orbit scan up to 2*10^9 and individually certified by a complete traversal of GF(p); completeness of the list above 10^6 is not established.
```

## LINKS (%H)

```text
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23420627">Unital Collapse Primes in Arithmetic Dynamics: Total Backward-Orbit Collapse for T(x) = (x^2 - 2)^2 mod p, Galois Stratification, and Fermat-Totient Asymptotics</a>, Zenodo preprint (2026).
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23856837">Total-Collapse Primes of Polynomial Maps over GF(p): Eighteen Families and Two Related Sequences</a>, Zenodo preprint (2026); <a href="https://github.com/Nondual-Observer/total-collapse-primes-oeis">GitHub verification suite</a>.
```

## FORMULA (%F)

```text
p is in the sequence iff #{x in {0..p-1} : K_1^(p)(x) in {-1, 0, 1} (mod p)} = p, where K_1(x) = (x^3 - 3*x)/2 mod p and K_1^(k) denotes the k-th iterate of K_1.
```

## EXAMPLE (%e)

```text
For p = 5: 2^(-1) == 3 (mod 5), and every element of GF(5) enters {0, 1, 4} = {0, 1, -1} mod 5 within 2 steps; thus 5 is in the sequence.
For p = 7: 2^(-1) == 4 (mod 7), and every element of GF(7) enters {0, 1, 6} within 2 steps; thus 7 is in the sequence.
```

## PROG (%o)

```text
(PARI)
isok(p) = {
  if (!isprime(p) || p < 5, return(0));
  if (kronecker(5, p) > 0, return(0));
  if (kronecker(7, p) == 1 && kronecker(-1, p) == 1, return(0));
  my(inv2 = (p + 1)\2, st = vector(p));
  st[1] = 2; st[2] = 2; st[p] = 2;
  for (x = 2, p\2,
    if (st[x+1] == 2, next);
    my(c = x, path = List());
    while (st[c+1] == 0, st[c+1] = 1; listput(path, c); c = ((c*(c^2 - 3)) % p * inv2) % p);
    if (st[c+1] == 1, return(0));
    for (i = 1, length(path), my(u = path[i]); st[u+1] = st[p - u + 1] = 2);
  );
  return(1);
}
(Julia)
isprime(n::Int) = n > 1 && all(d -> n % d != 0, 2:isqrt(n))
function isok(p::Int)
    (!isprime(p) || p < 5) && return false
    if p > 5
        powermod(5, p >> 1, p) == 1 && return false
        (powermod(7, p >> 1, p) == 1 && (p & 3) == 1) && return false
    end
    inv2, st = (p + 1) >> 1, zeros(UInt8, p)
    st[1] = st[2] = st[p] = 2
    f(c) = mod(mod(c * mod(c*c - 3, p), p) * inv2, p)
    for x in 2:(p >> 1)
        st[x + 1] == 2 && continue
        c = x; while st[c + 1] == 0; st[c + 1] = 1; c = f(c); end
        st[c + 1] == 1 && return false
        c = x; while st[c + 1] == 1; st[c + 1] = st[p - c + 1] = 2; c = f(c); end
    end
    return true
end
println([p for p in 5:5000 if isok(p)])
(Python)
from sympy import isprime, primerange
def isok(p: int) -> bool:
    if p < 5 or not isprime(p):
        return False
    if p > 5:
        if pow(5, p >> 1, p) == 1:
            return False
        if (p & 3) == 1 and pow(7, p >> 1, p) == 1:
            return False
    inv2 = (p + 1) >> 1
    state = bytearray(p)
    state[0] = state[1] = state[p - 1] = 2
    for x in range(2, (p >> 1) + 1):
        if state[x] == 2:
            continue
        curr, path = x, []
        while state[curr] == 0:
            state[curr] = 1
            path.append(curr)
            curr = ((curr * (curr * curr - 3)) % p * inv2) % p
        if state[curr] == 1:
            return False
        for u in path:
            state[u] = state[p - u] = 2
    return True
print([p for p in primerange(5, 5000) if isok(p)])
```

## CROSSREFS (%Y)

```text
Cf. A398745, A401173, A401667 (K_2). See the linked article for the companion sequences for K_n, n = 2, 4, 5, 6, 7, 8.
```

## KEYWORDS (%K)

```text
nonn,hard,more
```

## AUTHOR (%A)

```text
_Igor Zhuk_, Nov 01 2026
```
