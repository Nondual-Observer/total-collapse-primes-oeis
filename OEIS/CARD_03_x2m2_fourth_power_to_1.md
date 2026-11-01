# Sequence Card 03: (x^2 - 2)^4 -> {1}

* **Companion Article:** [`Total_Collapse_Primes_18_Families_Zhuk_2026.md`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.md) · [`PDF`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.pdf)
* **Verification Code:** [`../Code/`](../Code/)

## DATA (%S)

```text
3, 5, 13, 23, 41, 22381, 45503, 3852347, 1801981681
```

## NAME (%N)

```text
Primes p such that the iteration of T(x) = (x^2 - 2)^4 mod p collapses every element of GF(p) to 1.
```

## OFFSET (%O)

```text
1,1
```

## COMMENTS (%C)

```text
With h(x) = x^2 - 2 and g(u) = u^4 one has T = g o h, and F = h o g, F(u) = u^8 - 2; the cycles of T and F correspond bijectively via h, so p is a term iff every orbit of F mod p reaches -1 = h(1).
For p > 3 a necessary condition is that (F(u) - u)/(u + 1) = u^7 - u^6 + u^5 - u^4 + u^3 - u^2 + u - 2 (discriminant -26522311) has no root in GF(p). If p == 7 (mod 12), the basin of 1 is {1, -1}, so p is not a term.
The orbit of the critical point 0 over Z is 0 -> 16 -> 254^4 -> ..., and 254^4 - 1 = 3*5*11*17*23*149*433.
No other terms below 10^8 (exhaustive full functional-graph traversal of GF(p) for every prime p). The term 1801981681 was found in a bounded-step orbit scan up to 2*10^9 and certified by a complete traversal of GF(p); completeness of the list above 10^8 is not established.
```

## LINKS (%H)

```text
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23420627">Unital Collapse Primes in Arithmetic Dynamics: Total Backward-Orbit Collapse for T(x) = (x^2 - 2)^2 mod p, Galois Stratification, and Fermat-Totient Asymptotics</a>, Zenodo preprint (2026).
Igor M. Zhuk, <a href="https://github.com/Nondual-Observer/total-collapse-primes-oeis/blob/main/Article/Total_Collapse_Primes_18_Families_Zhuk_2026.pdf">Total-Collapse Primes of Polynomial Maps over GF(p): Eighteen Families and Two Related Sequences</a>, companion preprint and verification suite (2026).
```

## FORMULA (%F)

```text
p is in the sequence iff #{x in {0..p-1} : T^(p)(x) == 1 (mod p)} = p, where T(x) = (x^2 - 2)^4 mod p and T^(k) denotes the k-th iterate of T.
```

## EXAMPLE (%e)

```text
For p = 3: T(0) = 16 == 1 (mod 3), T(1) = T(2) = 1, so all 3 elements of GF(3) reach 1 in 1 step; thus 3 is in the sequence.
For p = 13: T(0) = 16 == 3 (mod 13), T(3) = 7^4 = 2401 == 9 (mod 13), T(9) = 79^4 = 1^4 = 1 (mod 13), and all 13 elements of GF(13) reach 1 within 3 steps; thus 13 is in the sequence.
```

## PROG (%o)

```text
(PARI)
isok(p) = {
  if (!isprime(p) || p == 2 || p%12 == 7, return(0));
  my(st = vector(p)); st[2] = st[p] = 2;
  for (x = 0, p\2,
    if (st[x+1] == 2, next);
    my(c = x, path = List());
    while (st[c+1] == 0, st[c+1] = 1; listput(path, c); c = ((c^2 - 2) % p)^4 % p);
    if (st[c+1] == 1, return(0));
    for (i = 1, length(path), my(u = path[i]); st[u+1] = 2; if (u, st[p - u + 1] = 2));
  );
  return(1);
}

(Julia)
isprime(n::Int) = n > 1 && all(d -> n % d != 0, 2:isqrt(n))
function isok(p::Int)
    (!isprime(p) || p == 2 || p % 12 == 7) && return false
    f(c) = (u = mod(c*c - 2, p); u2 = mod(u*u, p); mod(u2*u2, p))
    st = zeros(UInt8, p); st[2] = st[p] = 2
    for x in 0:(p >> 1)
        st[x + 1] == 2 && continue
        c = x; while st[c + 1] == 0; st[c + 1] = 1; c = f(c); end
        st[c + 1] == 1 && return false
        c = x; while st[c + 1] == 1; st[c + 1] = 2; c > 0 && (st[p - c + 1] = 2); c = f(c); end
    end
    return true
end
println([p for p in 3:50000 if isok(p)])

(Python)
from sympy import isprime

def isok(p: int) -> bool:
    if not isprime(p) or p == 2 or p % 12 == 7:
        return False
    state = bytearray(p)
    state[1] = state[p - 1] = 2
    for x in range((p >> 1) + 1):
        if state[x] == 2:
            continue
        curr, path = x, []
        while state[curr] == 0:
            state[curr] = 1
            path.append(curr)
            curr = pow(curr * curr - 2, 4, p)
        if state[curr] == 1:
            return False
        for u in path:
            state[u] = state[-u % p] = 2
    return True

print([p for p in range(2, 50000) if isok(p)])
```

## CROSSREFS (%Y)

```text
Cf. A019434, A401173.
```

## KEYWORDS (%K)

```text
nonn,hard,more
```

## AUTHOR (%A)

```text
_Igor M. Zhuk_, Oct 17 2026
```
