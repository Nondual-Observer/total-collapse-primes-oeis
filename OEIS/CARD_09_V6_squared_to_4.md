# Sequence Card 09: V_6(x)^2 -> value 4 (unique cycle through 4)

* **Companion Article:** [`Total_Collapse_Primes_18_Families_Zhuk_2026.md`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.md) | [`PDF`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.pdf)
* **Verification Code:** [`../Code/`](../Code/)

## DATA (%S)

```text
2, 3, 5, 7, 11, 13, 17, 37, 47, 79, 131, 167, 241, 257, 129401, 390209, 4926997, 6719077, 23734297
```

## NAME (%N)

```text
Primes p such that every trajectory of T(x) = (x^2*(x^2 - 3)^2 - 2)^2 mod p in GF(p) eventually reaches the critical boundary value 4 (mod p).
```

## OFFSET (%O)

```text
1,1
```

## COMMENTS (%C)

```text
Here V_6(x) = V_3(x)^2 - 2 = x^6 - 6*x^4 + 9*x^2 - 2 and T(x) = V_6(x)^2; T(0) = T(1) = T(-1) = T(2) = T(-2) = 4.
p is a term iff the functional graph of T on GF(p) is connected and its unique cycle contains 4. The point 4 is fixed iff p divides V_6(4)^2 - 4 = 2^6*3^3*5^2*13^2, i.e., p in {2, 3, 5, 13}. For the listed terms the length of the cycle through 4 is 1, 1, 1, 2, 2, 1, 2, 2, 2, 12, 11, 5, 5, 9, 399, 475, 1125, 1065, 2612.
No other terms below 10^6 (exhaustive full functional-graph traversal of GF(p) for every prime p). The terms above 10^6 were found in a bounded-step orbit scan up to 5*10^8 and individually certified by a complete traversal of GF(p); completeness of the list above 10^6 is not established.
```

## LINKS (%H)

```text
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23420627">Unital Collapse Primes in Arithmetic Dynamics: Total Backward-Orbit Collapse for T(x) = (x^2 - 2)^2 mod p, Galois Stratification, and Fermat-Totient Asymptotics</a>, Zenodo preprint (2026).
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23856837">Total-Collapse Primes of Polynomial Maps over GF(p): Eighteen Families and Two Related Sequences</a>, Zenodo preprint (2026); <a href="https://github.com/Nondual-Observer/total-collapse-primes-oeis">GitHub verification suite</a>.
```

## FORMULA (%F)

```text
p is in the sequence iff for every x in {0..p-1} there exists k >= 0 such that T^(k)(x) == 4 (mod p), where T(x) = (x^2*(x^2 - 3)^2 - 2)^2 mod p.
```

## EXAMPLE (%e)

```text
For p = 7: 4 mod 7 = 4; T(0) = T(1) = T(2) = T(5) = T(6) = 4 (mod 7), and T(3) = T(4) = 0 -> T(0) = 4 (mod 7), so every element of GF(7) passes through 4 within 2 steps; thus 7 is in the sequence.
For p = 37: every element of GF(37) reaches 4 within 6 steps under iteration of (x^2*(x^2 - 3)^2 - 2)^2 mod 37; thus 37 is in the sequence.
```

## PROG (%o)

```text
(PARI)
isok(p) = {
  if (!isprime(p), return(0));
  my(tgt = 4 % p, st = vector(p)); st[1] = st[tgt + 1] = 2;
  for (x = 1, p\2,
    if (st[x+1] == 2, next);
    my(c = x, path = List());
    while (st[c+1] == 0, st[c+1] = 1; listput(path, c); my(v3 = (c*(c^2 - 3)) % p, v6 = (v3^2 - 2) % p); c = (v6^2) % p);
    if (st[c+1] == 1, return(0));
    for (i = 1, length(path), my(u = path[i]); st[u+1] = st[p - u + 1] = 2);
  );
  return(1);
}
(Julia)
isprime(n::Int) = n > 1 && all(d -> n % d != 0, 2:isqrt(n))
function isok(p::Int)
    !isprime(p) && return false
    f(c) = (v3 = mod(c * mod(c*c - 3, p), p); v6 = mod(v3*v3 - 2, p); mod(v6*v6, p))
    st = zeros(UInt8, p); st[1] = st[mod(4, p) + 1] = 2
    for x in 1:(p >> 1)
        st[x + 1] == 2 && continue
        c = x; while st[c + 1] == 0; st[c + 1] = 1; c = f(c); end
        st[c + 1] == 1 && return false
        c = x; while st[c + 1] == 1; st[c + 1] = st[p - c + 1] = 2; c = f(c); end
    end
    return true
end
println([p for p in 2:3000 if isok(p)])
(Python)
from sympy import isprime
def isok(p: int) -> bool:
    if not isprime(p):
        return False
    state = bytearray(p)
    state[0] = state[4 % p] = 2
    for x in range(1, (p >> 1) + 1):
        if state[x] == 2:
            continue
        curr, path = x, []
        while state[curr] == 0:
            state[curr] = 1
            path.append(curr)
            v3 = (curr * (curr * curr - 3)) % p
            v6 = (v3 * v3 - 2) % p
            curr = (v6 * v6) % p
        if state[curr] == 1:
            return False
        for u in path:
            state[u] = state[p - u] = 2
    return True
print([p for p in range(2, 1000) if isok(p)])
```

## CROSSREFS (%Y)

```text
Cf. A401173. See the linked article for the sequences for V_m(x)^2, m = 3, 5, 7, 8, 9.
```

## KEYWORDS (%K)

```text
nonn,hard,more
```

## AUTHOR (%A)

```text
_Igor Zhuk_, Nov 01 2026
```
