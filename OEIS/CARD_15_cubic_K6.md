# Sequence Card 15: K_6(u) = (u^3 - 43u)/42 -> {-1, 0, 1}

* **Companion Article:** [`Total_Collapse_Primes_18_Families_Zhuk_2026.md`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.md) | [`PDF`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.pdf)
* **Verification Code:** [`../Code/`](../Code/)

## DATA (%S)

```text
5, 11, 17, 31, 41, 109, 127, 199, 719, 2731, 26317, 307067, 617579, 4348559, 4838117, 30282823, 42166321, 62371087
```

## NAME (%N)

```text
Primes p >= 5, p != 7, for which the iteration of K_6(u) = (u^3 - 43*u)/42 mod p on GF(p) collapses every element into the triad {-1, 0, 1} mod p.
```

## OFFSET (%O)

```text
1,1
```

## COMMENTS (%C)

```text
Over Z, 42*(K(u) - 1) = (u + 1)*(u + 6)*(u - 7) and 42*(K(u) + 1) = (u - 1)*(u - 6)*(u + 7), so {-1, 0, 1} is invariant: 0 is fixed and {1, -1} is a 2-cycle.
The 2-periodic points satisfy 42^4*(K(K(u)) - u) = u*(u^2 - 1)*(u^2 - 85)*(u^4 - 43*u^2 + 1764), with D_fix = 6^2 + 7^2 = 85 and 2-cycle discriminants 7^3 - 6^3 = 127 and 1 - 42 = -41. Hence every term p satisfies Kronecker(85, p) <= 0 and min(Kronecker(127, p), Kronecker(-41, p)) <= 0.
No other terms below 10^6 (exhaustive full functional-graph traversal of GF(p) for every prime p). The terms above 10^6 were found in a bounded-step orbit scan up to 2*10^9 and individually certified by a complete traversal of GF(p); completeness of the list above 10^6 is not established.
```

## LINKS (%H)

```text
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23420627">Unital Collapse Primes in Arithmetic Dynamics: Total Backward-Orbit Collapse for T(x) = (x^2 - 2)^2 mod p, Galois Stratification, and Fermat-Totient Asymptotics</a>, Zenodo preprint (2026).
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23856837">Total-Collapse Primes of Polynomial Maps over GF(p): Eighteen Families and Two Related Sequences</a>, Zenodo preprint (2026); <a href="https://github.com/Nondual-Observer/total-collapse-primes-oeis">GitHub verification suite</a>.
```

## FORMULA (%F)

```text
p is in the sequence iff #{x in {0..p-1} : K_6^(p)(x) in {-1, 0, 1} (mod p)} = p, where K_6(x) = (x^3 - 43*x)/42 mod p and K_6^(k) denotes the k-th iterate of K_6.
```

## EXAMPLE (%e)

```text
For p = 5: 42 == 2 (mod 5) and 43 == 3 (mod 5), so K_6(u) == (u^3 - 3*u)/2 (mod 5); all 5 elements of GF(5) enter {0, 1, 4} in 1 step, so 5 is in the sequence.
For p = 11: 42 == 9 == -2 (mod 11) and 43 == -1 (mod 11), so K_6(u) == 5*(u^3 + u) (mod 11); every element of GF(11) enters {0, 1, 10} within 3 steps, so 11 is in the sequence.
```

## PROG (%o)

```text
(PARI)
isok(p) = {
  if (!isprime(p) || p < 5 || p == 7, return(0));
  if (kronecker(85, p) > 0, return(0));
  if (kronecker(127, p) == 1 && kronecker(-41, p) == 1, return(0));
  my(inv42 = lift(Mod(1, p)/42), st = vector(p));
  st[1] = 2; st[2] = 2; st[p] = 2;
  for (x = 2, p\2,
    if (st[x+1] == 2, next);
    my(c = x, path = List());
    while (st[c+1] == 0, st[c+1] = 1; listput(path, c); c = ((c*(c^2 - 43)) % p * inv42) % p);
    if (st[c+1] == 1, return(0));
    for (i = 1, length(path), my(u = path[i]); st[u+1] = st[p - u + 1] = 2);
  );
  return(1);
}
(Julia)
isprime(n::Int) = n > 1 && all(d -> n % d != 0, 2:isqrt(n))
function isok(p::Int)
    (!isprime(p) || p < 5 || p == 7) && return false
    if p > 5
        powermod(85, p >> 1, p) == 1 && return false
        (powermod(127, p >> 1, p) == 1 && powermod(-41, p >> 1, p) == 1) && return false
    end
    inv42, st = invmod(42, p), zeros(UInt8, p)
    st[1] = st[2] = st[p] = 2
    f(c) = mod(mod(c * mod(c*c - 43, p), p) * inv42, p)
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
    if p < 5 or p == 7 or not isprime(p):
        return False
    if p > 5:
        if pow(85, p >> 1, p) == 1:
            return False
        if pow(127, p >> 1, p) == 1 and pow(-41, p >> 1, p) == 1:
            return False
    inv42 = pow(42, -1, p)
    state = bytearray(p)
    state[0] = state[1] = state[p - 1] = 2
    for x in range(2, (p >> 1) + 1):
        if state[x] == 2:
            continue
        curr, path = x, []
        while state[curr] == 0:
            state[curr] = 1
            path.append(curr)
            curr = ((curr * (curr * curr - 43)) % p * inv42) % p
        if state[curr] == 1:
            return False
        for u in path:
            state[u] = state[p - u] = 2
    return True
print([p for p in primerange(5, 3000) if isok(p)])
```

## CROSSREFS (%Y)

```text
Cf. A398745, A401173, A401667 (K_2). See the linked article for the companion sequences for K_n, n = 1, 2, 4, 5, 7, 8.
```

## KEYWORDS (%K)

```text
nonn,hard,more
```

## AUTHOR (%A)

```text
_Igor Zhuk_, Nov 01 2026
```
