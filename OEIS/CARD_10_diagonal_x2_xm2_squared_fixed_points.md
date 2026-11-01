# Sequence Card 10: x^2(x - 2)^2 -> fixed points

* **Companion Article:** [`Total_Collapse_Primes_18_Families_Zhuk_2026.md`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.md) · [`PDF`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.pdf)
* **Verification Code:** [`../Code/`](../Code/)

## DATA (%S)

```text
2, 3, 5, 7, 11, 23, 29, 31, 59, 181, 521, 1559, 2069, 5261, 6871, 9551, 19207, 1161371, 1574681, 1849721, 85268251
```

## NAME (%N)

```text
Primes p such that every trajectory of D(x) = x^2*(x - 2)^2 mod p in GF(p) eventually reaches a fixed point (i.e., D has no periodic cycles of length >= 2 in GF(p)).
```

## OFFSET (%O)

```text
1,1
```

## COMMENTS (%C)

```text
D(x) - x = x*(x - 1)*(x^2 - 3*x + 1), so 0 and 1 are fixed for every p, and for odd p the roots (3 +- sqrt(5))/2 are additional fixed points iff Kronecker(5, p) >= 0. Also D(2 - x) = D(x).
D is the restriction to the diagonal x = y of (x, y) -> ((x^2 - 2*y)^2, (y^2 - 2*x)^2), the coordinatewise square of the trace map (x, y) -> (x^2 - 2*y, y^2 - 2*x) of SL_3.
The terms with Kronecker(5, p) = -1 are 2, 3, 7, 23, 19207.
For the two-dimensional map on GF(p)^2, the primes p <= 5000 without cycles of length >= 2 are 2, 3, 7, 11, 29.
No other terms below 10^6 (exhaustive full functional-graph traversal of GF(p) for every prime p). The terms above 10^6 were found in a bounded-step orbit scan up to 5*10^8 and individually certified by a complete traversal of GF(p); completeness of the list above 10^6 is not established.
```

## LINKS (%H)

```text
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23420627">Unital Collapse Primes in Arithmetic Dynamics: Total Backward-Orbit Collapse for T(x) = (x^2 - 2)^2 mod p, Galois Stratification, and Fermat-Totient Asymptotics</a>, Zenodo preprint (2026).
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23856799">Total-Collapse Primes of Polynomial Maps over GF(p): Eighteen Families and Two Related Sequences</a>, Zenodo preprint (2026); <a href="https://github.com/Nondual-Observer/total-collapse-primes-oeis">GitHub verification suite</a>.
```

## FORMULA (%F)

```text
p is in the sequence iff #{x in {0..p-1} : D^(p+1)(x) == D^(p)(x) (mod p)} = p, where D(x) = x^2*(x - 2)^2 mod p and D^(k) denotes the k-th iterate of D.
```

## EXAMPLE (%e)

```text
For p = 7: D(0) = 0, D(1) = 1, D(2) = 0, D(3) = 2 -> 0, D(4) = 1, D(5) = 1, D(6) = 2 -> 0 (mod 7), so every element of GF(7) reaches one of the fixed points {0, 1} within 2 steps; thus 7 is in the sequence.
For p = 11: the fixed points of D(x) mod 11 are {0, 1, 5, 9} (where 5^2 - 3*5 + 1 = 11 == 0 mod 11), and every element of GF(11) reaches {0, 1, 5, 9} within 2 steps; thus 11 is in the sequence.
```

## PROG (%o)

```text
(PARI)
isok(p) = {
  if (!isprime(p), return(0));
  my(st = vector(p)); st[1] = st[2] = st[2%p + 1] = 2;
  for (x = 3, (p + 1)\2,
    if (st[x+1] == 2, next);
    my(c = x, path = List());
    while (st[c+1] == 0, st[c+1] = 1; listput(path, c); c = ((c*(c - 2)) % p)^2 % p);
    if (st[c+1] == 1 && (((c*(c - 2)) % p)^2 % p) != c, return(0));
    for (i = 1, length(path), my(u = path[i]); st[u+1] = st[(p + 2 - u)%p + 1] = 2);
  );
  return(1);
}

(Julia)
isprime(n::Int) = n > 1 && all(d -> n % d != 0, 2:isqrt(n))
function isok(p::Int)
    !isprime(p) && return false
    f(c) = mod(mod(c*(c - 2), p)^2, p)
    st = zeros(UInt8, p); st[1] = st[2] = st[mod(2, p) + 1] = 2
    for x in 3:((p + 1) >> 1)
        st[x + 1] == 2 && continue
        c = x; while st[c + 1] == 0; st[c + 1] = 1; c = f(c); end
        (st[c + 1] == 1 && f(c) != c) && return false
        c = x; while st[c + 1] == 1; st[c + 1] = st[mod(2 - c, p) + 1] = 2; c = f(c); end
    end
    return true
end
println([p for p in 2:20000 if isok(p)])

(Python)
from sympy import isprime

def isok(p: int) -> bool:
    if not isprime(p):
        return False
    state = bytearray(p)
    state[0] = state[1] = state[2 % p] = 2
    for x in range(3, ((p + 1) >> 1) + 1):
        if state[x] == 2:
            continue
        curr, path = x, []
        while state[curr] == 0:
            state[curr] = 1
            path.append(curr)
            t = (curr * (curr - 2)) % p
            curr = (t * t) % p
        if state[curr] == 1:
            t = (curr * (curr - 2)) % p
            if (t * t) % p != curr:
                return False
        for u in path:
            state[u] = state[(2 - u) % p] = 2
    return True

print([p for p in range(2, 20000) if isok(p)])
```

## CROSSREFS (%Y)

```text
Cf. A401173.
```

## KEYWORDS (%K)

```text
nonn,hard,more
```

## AUTHOR (%A)

```text
_Igor M. Zhuk_, Oct 17 2026
```
