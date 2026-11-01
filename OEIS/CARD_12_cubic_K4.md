# Sequence Card 12: K_4(u) = (u^3 - 21u)/20 -> {-1, 0, 1}

* **Companion Article:** [`Total_Collapse_Primes_18_Families_Zhuk_2026.md`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.md) · [`PDF`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.pdf)
* **Verification Code:** [`../Code/`](../Code/)

## DATA (%S)

```text
7, 11, 17, 19, 29, 41, 53, 79, 191, 42979, 277829, 35718901
```

## NAME (%N)

```text
Primes p >= 7 for which the iteration of K_4(u) = (u^3 - 21*u)/20 mod p on GF(p) collapses every element into the triad {-1, 0, 1} mod p.
```

## OFFSET (%O)

```text
1,1
```

## COMMENTS (%C)

```text
Over Z, 20*(K(u) - 1) = (u + 1)*(u + 4)*(u - 5) and 20*(K(u) + 1) = (u - 1)*(u - 4)*(u + 5), so {-1, 0, 1} is invariant: 0 is fixed and {1, -1} is a 2-cycle.
The 2-periodic points satisfy 20^4*(K(K(u)) - u) = u*(u^2 - 1)*(u^2 - 41)*(u^4 - 21*u^2 + 400), with D_fix = 4^2 + 5^2 = 41 and 2-cycle discriminants 5^3 - 4^3 = 61 and 1 - 20 = -19. Hence every term p satisfies Kronecker(41, p) <= 0 and min(Kronecker(61, p), Kronecker(-19, p)) <= 0.
No other terms below 10^6 (exhaustive full functional-graph traversal of GF(p) for every prime p). The terms above 10^6 were found in a bounded-step orbit scan up to 2*10^9 and individually certified by a complete traversal of GF(p); completeness of the list above 10^6 is not established.
```

## LINKS (%H)

```text
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23420627">Unital Collapse Primes in Arithmetic Dynamics: Total Backward-Orbit Collapse for T(x) = (x^2 - 2)^2 mod p, Galois Stratification, and Fermat-Totient Asymptotics</a>, Zenodo preprint (2026).
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23856837">Total-Collapse Primes of Polynomial Maps over GF(p): Eighteen Families and Two Related Sequences</a>, Zenodo preprint (2026); <a href="https://github.com/Nondual-Observer/total-collapse-primes-oeis">GitHub verification suite</a>.
```

## FORMULA (%F)

```text
p is in the sequence iff #{x in {0..p-1} : K_4^(p)(x) in {-1, 0, 1} (mod p)} = p, where K_4(x) = (x^3 - 21*x)/20 mod p and K_4^(k) denotes the k-th iterate of K_4.
```

## EXAMPLE (%e)

```text
For p = 7: 20 == -1 (mod 7) and 21 == 0 (mod 7), so K_4(u) == -u^3 (mod 7), which maps all of GF(7) onto {0, 1, 6} = {0, 1, -1} mod 7 in one step; thus 7 is in the sequence.
For p = 11: 20^(-1) == 5 (mod 11), so K_4(u) == 5*(u^3 + u) (mod 11). Every element of GF(11) eventually enters {0, 1, 10}, so 11 is in the sequence.
```

## PROG (%o)

```text
(PARI)
isok(p) = {
  if (!isprime(p) || p < 7, return(0));
  if (kronecker(41, p) > 0, return(0));
  if (kronecker(61, p) == 1 && kronecker(-19, p) == 1, return(0));
  my(inv20 = lift(Mod(1, p)/20), st = vector(p));
  st[1] = 2; st[2] = 2; st[p] = 2;
  for (x = 2, p\2,
    if (st[x+1] == 2, next);
    my(c = x, path = List());
    while (st[c+1] == 0, st[c+1] = 1; listput(path, c); c = ((c*(c^2 - 21)) % p * inv20) % p);
    if (st[c+1] == 1, return(0));
    for (i = 1, length(path), my(u = path[i]); st[u+1] = st[p - u + 1] = 2);
  );
  return(1);
}

(Julia)
isprime(n::Int) = n > 1 && all(d -> n % d != 0, 2:isqrt(n))
function isok(p::Int)
    (!isprime(p) || p < 7) && return false
    powermod(41, p >> 1, p) == 1 && return false
    (powermod(61, p >> 1, p) == 1 && powermod(-19, p >> 1, p) == 1) && return false
    inv20, st = invmod(20, p), zeros(UInt8, p)
    st[1] = st[2] = st[p] = 2
    f(c) = mod(mod(c * mod(c*c - 21, p), p) * inv20, p)
    for x in 2:(p >> 1)
        st[x + 1] == 2 && continue
        c = x; while st[c + 1] == 0; st[c + 1] = 1; c = f(c); end
        st[c + 1] == 1 && return false
        c = x; while st[c + 1] == 1; st[c + 1] = st[p - c + 1] = 2; c = f(c); end
    end
    return true
end
println([p for p in 7:3000 if isok(p)])

(Python)
from sympy import isprime, primerange

def isok(p: int) -> bool:
    if p < 7 or not isprime(p) or pow(41, p >> 1, p) == 1:
        return False
    if pow(61, p >> 1, p) == 1 and pow(-19, p >> 1, p) == 1:
        return False
    inv20 = pow(20, -1, p)
    state = bytearray(p)
    state[0] = state[1] = state[p - 1] = 2
    for x in range(2, (p >> 1) + 1):
        if state[x] == 2:
            continue
        curr, path = x, []
        while state[curr] == 0:
            state[curr] = 1
            path.append(curr)
            curr = ((curr * (curr * curr - 21)) % p * inv20) % p
        if state[curr] == 1:
            return False
        for u in path:
            state[u] = state[p - u] = 2
    return True

print([p for p in primerange(7, 3000) if isok(p)])
```

## CROSSREFS (%Y)

```text
Cf. A401173. See the linked article for the companion sequences for K_n, n = 1, 2, 5, 6, 7, 8.
```

## KEYWORDS (%K)

```text
nonn,hard,more
```

## AUTHOR (%A)

```text
_Igor M. Zhuk_, Oct 18 2026
```
