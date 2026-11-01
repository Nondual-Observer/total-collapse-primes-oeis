# Sequence Card 05: V_8(x)^2 -> {1}

* **Companion Article:** [`Total_Collapse_Primes_18_Families_Zhuk_2026.md`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.md) · [`PDF`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.pdf)
* **Verification Code:** [`../Code/`](../Code/)

## DATA (%S)

```text
3, 5, 13, 17, 23, 47, 71, 89, 193, 313, 863, 1103, 1777, 107183, 84814153, 342037921
```

## NAME (%N)

```text
Primes p such that the iteration of T_8(x) = (((x^2 - 2)^2 - 2)^2 - 2)^2 mod p collapses every element of GF(p) to 1.
```

## OFFSET (%O)

```text
1,1
```

## COMMENTS (%C)

```text
Here V_8(x) = ((x^2 - 2)^2 - 2)^2 - 2 = 2*T_8(x/2) and T(x) = V_8(x)^2 has degree 16; V_8(0) = 2 and V_8(1) = -1, so 1 is a fixed point of T.
With F(u) = V_8(u^2) one has V_8(T(x)) = F(V_8(x)); the cycles of T and F correspond via V_8, so p is a term iff every orbit of F mod p reaches -1.
If p == 7 (mod 12), the basin of 1 is {1, -1}, so p is not a term.
No other terms below 10^6 (exhaustive full functional-graph traversal of GF(p) for every prime p). The terms above 10^6 were found in a bounded-step orbit scan up to 2*10^9 and individually certified by a complete traversal of GF(p); completeness of the list above 10^6 is not established.
```

## LINKS (%H)

```text
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23420627">Unital Collapse Primes in Arithmetic Dynamics: Total Backward-Orbit Collapse for T(x) = (x^2 - 2)^2 mod p, Galois Stratification, and Fermat-Totient Asymptotics</a>, Zenodo preprint (2026).
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23856799">Total-Collapse Primes of Polynomial Maps over GF(p): Eighteen Families and Two Related Sequences</a>, Zenodo preprint (2026); <a href="https://github.com/Nondual-Observer/total-collapse-primes-oeis">GitHub verification suite</a>.
```

## FORMULA (%F)

```text
p is in the sequence iff #{x in {0..p-1} : T_8^(p)(x) == 1 (mod p)} = p, where T_8(x) = (((x^2 - 2)^2 - 2)^2 - 2)^2 mod p and T_8^(k) denotes the k-th iterate of T_8.
```

## EXAMPLE (%e)

```text
For p = 3: T_8(0) = 4 == 1 (mod 3), T_8(1) = T_8(2) = 1, so all 3 elements of GF(3) reach 1 in 1 step; thus 3 is in the sequence.
For p = 17: every element of GF(17) reaches 1 within 4 steps under iteration of T_8(x) mod 17; thus 17 is in the sequence.
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
    while (st[c+1] == 0, st[c+1] = 1; listput(path, c); my(v2 = (c^2 - 2) % p, v4 = (v2^2 - 2) % p, v8 = (v4^2 - 2) % p); c = (v8^2) % p);
    if (st[c+1] == 1, return(0));
    for (i = 1, length(path), my(u = path[i]); st[u+1] = 2; if (u, st[p - u + 1] = 2));
  );
  return(1);
}

(Julia)
isprime(n::Int) = n > 1 && all(d -> n % d != 0, 2:isqrt(n))
function isok(p::Int)
    (!isprime(p) || p == 2 || p % 12 == 7) && return false
    f(c) = (v2 = mod(c*c - 2, p); v4 = mod(v2*v2 - 2, p); v8 = mod(v4*v4 - 2, p); mod(v8*v8, p))
    st = zeros(UInt8, p); st[2] = st[p] = 2
    for x in 0:(p >> 1)
        st[x + 1] == 2 && continue
        c = x; while st[c + 1] == 0; st[c + 1] = 1; c = f(c); end
        st[c + 1] == 1 && return false
        c = x; while st[c + 1] == 1; st[c + 1] = 2; c > 0 && (st[p - c + 1] = 2); c = f(c); end
    end
    return true
end
println([p for p in 3:3000 if isok(p)])

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
            v2 = (curr * curr - 2) % p
            v4 = (v2 * v2 - 2) % p
            v8 = (v4 * v4 - 2) % p
            curr = (v8 * v8) % p
        if state[curr] == 1:
            return False
        for u in path:
            state[u] = state[-u % p] = 2
    return True

print([p for p in range(2, 2000) if isok(p)])
```

## CROSSREFS (%Y)

```text
Cf. A019434, A401173 (V_2(x)^2).
```

## KEYWORDS (%K)

```text
nonn,hard,more
```

## AUTHOR (%A)

```text
_Igor M. Zhuk_, Oct 17 2026
```
