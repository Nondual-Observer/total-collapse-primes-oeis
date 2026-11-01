# Sequence Card 18: Least collapsing prime p > n + 1 for K_n (or 0 if none)

* **Companion Article:** [`Total_Collapse_Primes_18_Families_Zhuk_2026.md`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.md) · [`PDF`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.pdf)
* **Verification Code:** [`../Code/`](../Code/)

## DATA (%S)

```text
5, 5, 5, 7, 7, 11, 19, 19, 41, 13, 13, 17, 127, 19, 23, 19, 23, 31, 31, 0, 23, 31, 53, 29, 37, 89, 89, 31, 41, 67, 37, 41, 53, 41, 53, 41, 47, 83, 211, 53, 61, 59, 991, 47, 8713, 173, 199, 53, 389, 89, 67, 61, 73, 239, 109
```

## NAME (%N)

```text
Least prime p >= 5 with p > n + 1 for which the iteration of K_n(u) = (u^3 - (n^2 + n + 1)*u)/(n*(n + 1)) mod p on GF(p) collapses every element into {-1, 0, 1} mod p, or 0 if no such prime exists.
```

## OFFSET (%O)

```text
1,1
```

## COMMENTS (%C)

```text
For n >= 1 let B_n = n*(n + 1) and A_n = n^2 + n + 1. Then B_n*(K_n(u) - 1) = (u + 1)*(u + n)*(u - n - 1) and B_n*(K_n(u) + 1) = (u - 1)*(u - n)*(u + n + 1), so {-1, 0, 1} is invariant. The condition p > n + 1 guarantees that p does not divide B_n.
B_n^4*(K_n(K_n(u)) - u) = u*(u^2 - 1)*(u^2 - D)*(u^4 - A_n*u^2 + B_n^2) with D = n^2 + (n + 1)^2 = A001844(n); the 2-cycle discriminants are (n + 1)^3 - n^3 = A003215(n) and 1 - n*(n + 1).
If n^2 + (n + 1)^2 = c^2 (n in A001652: 3, 20, 119, 696, 4059, ...), then +-c are fixed points of K_n for every prime p > n + 1. Since c^2 - 1 = 2*n*(n + 1) has only prime factors <= n + 1 and c < 2*(n + 1), these fixed points lie in {-1, 0, 1} (mod p) iff c is prime and p = c, where K_n(u) == u*(1 - 2*u^2) (mod c). Hence a(n) = 0 if c is composite (e.g., n = 119, c = 169; n = 696, c = 985) and a(n) is 0 or c if c is prime: a(3) = 5, while a(20) = 0 (4-cycle 22 -> 12 -> 7 -> 17 -> 22 mod 29), a(4059) = 0 (c = 5741) and a(23660) = 0 (c = 33461) by direct computation. Conjecture: a(n) = 0 for all n in A001652 with n > 3.
K_n mod p depends only on n mod p; e.g., K_119 == K_2 (mod 13) and 13 collapses for n = 2, which is why the definition requires p > n + 1.
For n not in A001652 it is not known whether a suitable prime always exists; the programs below terminate whenever one exists.
```

## LINKS (%H)

```text
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23420627">Unital Collapse Primes in Arithmetic Dynamics: Total Backward-Orbit Collapse for T(x) = (x^2 - 2)^2 mod p, Galois Stratification, and Fermat-Totient Asymptotics</a>, Zenodo preprint (2026).
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23856837">Total-Collapse Primes of Polynomial Maps over GF(p): Eighteen Families and Two Related Sequences</a>, Zenodo preprint (2026); <a href="https://github.com/Nondual-Observer/total-collapse-primes-oeis">GitHub verification suite</a>.
```

## FORMULA (%F)

```text
a(n) = min{p prime : p >= 5, p > n + 1, and #{x in {0..p-1} : K_n^(p)(x) in {-1, 0, 1} (mod p)} = p}, or 0 if no such prime exists (in particular, if n^2 + (n + 1)^2 = c^2 is a square, only p = c with c prime can satisfy the condition).
Every prime p = a(n) > 0 satisfies Kronecker(n^2 + (n + 1)^2, p) <= 0 and min(Kronecker((n + 1)^3 - n^3, p), Kronecker(1 - n*(n + 1), p)) <= 0.
```

## EXAMPLE (%e)

```text
For n = 1: K_1(u) = (u^3 - 3*u)/2; the least prime p >= 5 with p > 2 collapsing GF(p) onto {-1, 0, 1} is p = 5, so a(1) = 5.
For n = 2: K_2(u) = (u^3 - 7*u)/6; the least prime p > 3 collapsing GF(p) onto {-1, 0, 1} is p = 5, so a(2) = 5.
For n = 3: 3^2 + 4^2 = 25 = 5^2; for p = 5 > 4 we have K_3(u) == 3*u^3 + u (mod 5) which sends 0 -> 0, 1 <-> 4, 2 -> 1, 3 -> 4 (collapsing GF(5) onto {-1, 0, 1}), whereas for every prime p > 5 the residues u = +-5 (mod p) are fixed points outside {-1, 0, 1}; thus p = 5 is the unique collapsing prime and a(3) = 5.
For n = 13: K_13(u) = (u^3 - 183*u)/182; the least prime p > 14 collapsing GF(p) onto {-1, 0, 1} is the Mersenne prime a(13) = 127.
For n = 20: 20^2 + 21^2 = 841 = 29^2; for every prime p > 21 with p != 29 (in particular for p = 23, where +-29 == +-6 mod 23), u = +-29 (mod p) are fixed points outside {-1, 0, 1}, while at p = 29 the map K_20(u) == 27*u*(u^2 - 15) (mod 29) has the 4-cycle 22 -> 12 -> 7 -> 17 -> 22 outside {-1, 0, 1}, so a(20) = 0.
```

## PROG (%o)

```text
(PARI)
isok(p, n) = {
  my(bn = (n*(n+1)) % p, dfix = n^2 + (n+1)^2);
  if (bn == 0, return(0));
  if (kronecker(dfix, p) > 0, return(0));
  my(an = (n^2 + n + 1) % p, invb = lift(Mod(1, p)/bn), st = vector(p));
  st[1] = 2; st[2] = 2; st[p] = 2;
  for (x = 2, p\2,
    if (st[x+1] == 2, next);
    my(c = x, path = List());
    while (st[c+1] == 0, st[c+1] = 1; listput(path, c); c = ((c*(c^2 - an)) % p * invb) % p);
    if (st[c+1] == 1, return(0));
    for (i = 1, length(path), my(u = path[i]); st[u+1] = st[p - u + 1] = 2);
  );
  return(1);
}
a(n) = {
  my(dfix = n^2 + (n+1)^2, r = sqrtint(dfix));
  if (r^2 == dfix, return(if(isprime(r) && isok(r, n), r, 0)));
  forprime (p = max(5, n + 2), oo, if (isok(p, n), return(p)));
}

(Julia)
isprime(n::Int) = n > 1 && all(d -> n % d != 0, 2:isqrt(n))
function isok(p::Int, n::Int)
    bn = mod(n * (n + 1), p)
    bn == 0 && return false
    dfix = n * n + (n + 1) * (n + 1)
    (dfix % p != 0 && powermod(dfix, p >> 1, p) == 1) && return false
    an, invb, st = mod(n * n + n + 1, p), invmod(bn, p), zeros(UInt8, p)
    st[1] = st[2] = st[p] = 2
    f(c) = mod(mod(c * mod(c * c - an, p), p) * invb, p)
    for x in 2:(p >> 1)
        st[x + 1] == 2 && continue
        c = x; while st[c + 1] == 0; st[c + 1] = 1; c = f(c); end
        st[c + 1] == 1 && return false
        c = x; while st[c + 1] == 1; st[c + 1] = st[p - c + 1] = 2; c = f(c); end
    end
    return true
end
function a(n::Int)
    dfix = n * n + (n + 1) * (n + 1)
    r = isqrt(dfix)
    r * r == dfix && return (isprime(r) && isok(r, n)) ? r : 0
    p = max(5, n + 2)
    while !(isprime(p) && isok(p, n)); p += 1; end
    return p
end
println([a(n) for n in 1:55])

(Python)
import math
from sympy import isprime

def isok(p: int, n: int) -> bool:
    bn = (n * (n + 1)) % p
    if bn == 0:
        return False
    dfix = n * n + (n + 1) * (n + 1)
    if dfix % p != 0 and pow(dfix, p >> 1, p) == 1:
        return False
    an = (n * n + n + 1) % p
    invb = pow(bn, -1, p)
    st = bytearray(p)
    st[0] = st[1] = st[p - 1] = 2
    for x in range(2, (p >> 1) + 1):
        if st[x] == 2:
            continue
        c, path = x, []
        while st[c] == 0:
            st[c] = 1
            path.append(c)
            c = ((c * (c * c - an)) % p * invb) % p
        if st[c] == 1:
            return False
        for u in path:
            st[u] = st[p - u] = 2
    return True

def a(n: int) -> int:
    dfix = n * n + (n + 1) * (n + 1)
    r = math.isqrt(dfix)
    if r * r == dfix:
        return r if (isprime(r) and isok(r, n)) else 0
    p = max(5, n + 2)
    while not (isprime(p) and isok(p, n)):
        p += 1
    return p

print([a(n) for n in range(1, 56)])
```

## CROSSREFS (%Y)

```text
Cf. A001652, A001653, A001844, A002315, A003215, A086383.
```

## KEYWORDS (%K)

```text
nonn
```

## AUTHOR (%A)

```text
_Igor M. Zhuk_, Oct 19 2026
```
