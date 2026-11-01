# Sequence Card 16: K_7(u) = (u^3 - 57u)/56 -> {-1, 0, 1}

* **Companion Article:** [`Total_Collapse_Primes_18_Families_Zhuk_2026.md`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.md) · [`PDF`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.pdf)
* **Verification Code:** [`../Code/`](../Code/)

## DATA (%S)

```text
5, 19, 23, 151, 1237, 11941, 761357, 6756377, 1105256599
```

## NAME (%N)

```text
Primes p >= 5, p != 7, for which the iteration of K_7(u) = (u^3 - 57*u)/56 mod p on GF(p) collapses every element into the triad {-1, 0, 1} mod p.
```

## OFFSET (%O)

```text
1,1
```

## COMMENTS (%C)

```text
Over Z, 56*(K(u) - 1) = (u + 1)*(u + 7)*(u - 8) and 56*(K(u) + 1) = (u - 1)*(u - 7)*(u + 8), so {-1, 0, 1} is invariant: 0 is fixed and {1, -1} is a 2-cycle.
Since 3*56 + 1 = 169 = 13^2, the quartic factor of the 2-periodic points splits over Z: 56^4*(K(K(u)) - u) = u*(u^2 - 1)*(u^2 - 113)*(u^2 - 13*u + 56)*(u^2 + 13*u + 56), the quadratics having discriminant -55. Hence every term p satisfies Kronecker(113, p) <= 0 and Kronecker(-55, p) <= 0 (both equal -1 for all known terms p > 5). For p = 11, +-5 are fixed points outside {-1, 0, 1}; for p = 113 there is the 9-cycle 17 -> 22 -> 83 -> 69 -> 33 -> 27 -> 98 -> 68 -> 49 -> 17.
No other terms below 10^6 (exhaustive full functional-graph traversal of GF(p) for every prime p). The terms above 10^6 were found in a bounded-step orbit scan up to 2*10^9 and individually certified by a complete traversal of GF(p); completeness of the list above 10^6 is not established.
```

## LINKS (%H)

```text
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23420627">Unital Collapse Primes in Arithmetic Dynamics: Total Backward-Orbit Collapse for T(x) = (x^2 - 2)^2 mod p, Galois Stratification, and Fermat-Totient Asymptotics</a>, Zenodo preprint (2026).
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23856837">Total-Collapse Primes of Polynomial Maps over GF(p): Eighteen Families and Two Related Sequences</a>, Zenodo preprint (2026); <a href="https://github.com/Nondual-Observer/total-collapse-primes-oeis">GitHub verification suite</a>.
```

## FORMULA (%F)

```text
p is in the sequence iff #{x in {0..p-1} : K_7^(p)(x) in {-1, 0, 1} (mod p)} = p, where K_7(x) = (x^3 - 57*x)/56 mod p and K_7^(k) denotes the k-th iterate of K_7.
```

## EXAMPLE (%e)

```text
For p = 5: 56 == 1 (mod 5) and 57 == 2 (mod 5), so K_7(u) == u^3 - 2*u (mod 5), collapsing GF(5) onto {0, 1, 4} in 1 step; thus 5 is in the sequence.
For p = 19: 57 == 0 (mod 19), so K_7(u) == 56^(-1)*u^3 == -u^3 (mod 19), collapsing GF(19) onto {0, 1, 18} within 2 steps; thus 19 is in the sequence.
```

## PROG (%o)

```text
(PARI)
isok(p) = {
  if (!isprime(p) || p < 5 || p == 7, return(0));
  if (p > 5 && (kronecker(113, p) != -1 || kronecker(-55, p) != -1), return(0));
  my(inv56 = lift(Mod(1, p)/56), st = vector(p));
  st[1] = 2; st[2] = 2; st[p] = 2;
  for (x = 2, p\2,
    if (st[x+1] == 2, next);
    my(c = x, path = List());
    while (st[c+1] == 0, st[c+1] = 1; listput(path, c); c = ((c*(c^2 - 57)) % p * inv56) % p);
    if (st[c+1] == 1, return(0));
    for (i = 1, length(path), my(u = path[i]); st[u+1] = st[p - u + 1] = 2);
  );
  return(1);
}

(Julia)
isprime(n::Int) = n > 1 && all(d -> n % d != 0, 2:isqrt(n))
function isok(p::Int)
    (!isprime(p) || p < 5 || p == 7) && return false
    (p > 5 && (powermod(113, p >> 1, p) != p - 1 || powermod(-55, p >> 1, p) != p - 1)) && return false
    inv56, st = invmod(56, p), zeros(UInt8, p)
    st[1] = st[2] = st[p] = 2
    f(c) = mod(mod(c * mod(c*c - 57, p), p) * inv56, p)
    for x in 2:(p >> 1)
        st[x + 1] == 2 && continue
        c = x; while st[c + 1] == 0; st[c + 1] = 1; c = f(c); end
        st[c + 1] == 1 && return false
        c = x; while st[c + 1] == 1; st[c + 1] = st[p - c + 1] = 2; c = f(c); end
    end
    return true
end
println([p for p in 5:15000 if isok(p)])

(Python)
from sympy import isprime, primerange

def isok(p: int) -> bool:
    if p < 5 or p == 7 or not isprime(p):
        return False
    if p > 5 and (pow(113, p >> 1, p) != p - 1 or pow(-55, p >> 1, p) != p - 1):
        return False
    inv56 = pow(56, -1, p)
    state = bytearray(p)
    state[0] = state[1] = state[p - 1] = 2
    for x in range(2, (p >> 1) + 1):
        if state[x] == 2:
            continue
        curr, path = x, []
        while state[curr] == 0:
            state[curr] = 1
            path.append(curr)
            curr = ((curr * (curr * curr - 57)) % p * inv56) % p
        if state[curr] == 1:
            return False
        for u in path:
            state[u] = state[p - u] = 2
    return True

print([p for p in primerange(5, 15000) if isok(p)])
```

## CROSSREFS (%Y)

```text
Cf. A401173. See the linked article for the companion sequences for K_n, n = 1, 2, 4, 5, 6, 8.
```

## KEYWORDS (%K)

```text
nonn,hard,more
```

## AUTHOR (%A)

```text
_Igor M. Zhuk_, Oct 18 2026
```
