# Sequence Card 11: K_2(u) = (u^3 - 7u)/6 -> {-1, 0, 1}

* **Companion Article:** [`Total_Collapse_Primes_18_Families_Zhuk_2026.md`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.md) · [`PDF`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.pdf)
* **Verification Code:** [`../Code/`](../Code/)

## DATA (%S)

```text
5, 7, 13, 19, 31, 47, 59, 73, 151, 331, 401, 1943069, 83584483, 158450797, 1554906533
```

## NAME (%N)

```text
Primes p >= 5 for which the iteration of K(u) = (u^3 - 7*u)/6 mod p on GF(p) collapses every element into the triad {-1, 0, 1} mod p.
```

## OFFSET (%O)

```text
1,1
```

## COMMENTS (%C)

```text
Over Z, 6*(K(u) - 1) = (u + 1)*(u + 2)*(u - 3) and 6*(K(u) + 1) = (u - 1)*(u - 2)*(u + 3), so {-1, 0, 1} is invariant: 0 is fixed and {1, -1} is a 2-cycle.
The 2-periodic points satisfy 6^4*(K(K(u)) - u) = u*(u^2 - 1)*(u^2 - 13)*(u^4 - 7*u^2 + 36), with D_fix = 2^2 + 3^2 = 13 and 2-cycle discriminants 3^3 - 2^3 = 19 and 1 - 6 = -5. Hence every term p satisfies Kronecker(13, p) <= 0 and min(Kronecker(19, p), Kronecker(-5, p)) <= 0.
Over Z, K maps {-3, ..., 3} into {-1, 0, 1}.
No other terms below 10^6 (exhaustive full functional-graph traversal of GF(p) for every prime p). The terms above 10^6 were found in a bounded-step orbit scan up to 2*10^9 and individually certified by a complete traversal of GF(p); completeness of the list above 10^6 is not established.
```

## LINKS (%H)

```text
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23420627">Unital Collapse Primes in Arithmetic Dynamics: Total Backward-Orbit Collapse for T(x) = (x^2 - 2)^2 mod p, Galois Stratification, and Fermat-Totient Asymptotics</a>, Zenodo preprint (2026).
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23856799">Total-Collapse Primes of Polynomial Maps over GF(p): Eighteen Families and Two Related Sequences</a>, Zenodo preprint (2026); <a href="https://github.com/Nondual-Observer/total-collapse-primes-oeis">GitHub verification suite</a>.
```

## FORMULA (%F)

```text
p is in the sequence iff #{x in {0..p-1} : K^(p)(x) in {-1, 0, 1} (mod p)} = p, where K(x) = (x^3 - 7*x)/6 mod p and K^(k) denotes the k-th iterate of K.
```

## EXAMPLE (%e)

```text
For p = 5: 6^(-1) == 1 (mod 5), so K(u) == u^3 - 2*u (mod 5). The map is 0 -> 0, 1 -> 4, 4 -> 1, 2 -> 4, 3 -> 1; all elements enter {0, 1, 4} = {0, 1, -1} mod 5 in 1 step, so 5 is in the sequence.
For p = 7: 6^(-1) == -1 (mod 7), so K(u) == -u^3 (mod 7). Since u^3 in {0, 1, 6} mod 7 for all u in GF(7), the entire field collapses in 1 step to {0, 1, 6}, so 7 is in the sequence.
For p = 11: 6^(-1) == 2 (mod 11), so K(u) == 2*u*(u^2 - 7) (mod 11), and 5 -> 4 -> 6 -> 7 -> 5 is an extraneous 4-cycle outside {0, 1, 10}, so 11 is not in the sequence.
```

## PROG (%o)

```text
(PARI)
isok(p) = {
  if (!isprime(p) || p < 5, return(0));
  if (kronecker(13, p) > 0, return(0));
  if (kronecker(19, p) == 1 && kronecker(-5, p) == 1, return(0));
  my(inv6 = lift(Mod(1, p)/6), st = vector(p));
  st[1] = 2; st[2] = 2; st[p] = 2;
  for (x = 2, p\2,
    if (st[x+1] == 2, next);
    my(c = x, path = List());
    while (st[c+1] == 0, st[c+1] = 1; listput(path, c); c = ((c*(c^2 - 7)) % p * inv6) % p);
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
        powermod(13, p >> 1, p) == 1 && return false
        (powermod(19, p >> 1, p) == 1 && powermod(-5, p >> 1, p) == 1) && return false
    end
    inv6, st = invmod(6, p), zeros(UInt8, p)
    st[1] = st[2] = st[p] = 2
    f(c) = mod(mod(c * mod(c*c - 7, p), p) * inv6, p)
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
        if pow(13, p >> 1, p) == 1:
            return False
        if pow(19, p >> 1, p) == 1 and pow(-5, p >> 1, p) == 1:
            return False
    inv6 = pow(6, -1, p)
    state = bytearray(p)
    state[0] = state[1] = state[p - 1] = 2
    for x in range(2, (p >> 1) + 1):
        if state[x] == 2:
            continue
        curr, path = x, []
        while state[curr] == 0:
            state[curr] = 1
            path.append(curr)
            curr = ((curr * (curr * curr - 7)) % p * inv6) % p
        if state[curr] == 1:
            return False
        for u in path:
            state[u] = state[p - u] = 2
    return True

print([p for p in primerange(5, 5000) if isok(p)])
```

## CROSSREFS (%Y)

```text
Cf. A401173. See the linked article for the companion sequences for K_n(u) = (u^3 - (n^2 + n + 1)*u)/(n*(n + 1)), n = 1, 4, 5, 6, 7, 8.
```

## KEYWORDS (%K)

```text
nonn,hard,more
```

## AUTHOR (%A)

```text
_Igor M. Zhuk_, Oct 18 2026
```
