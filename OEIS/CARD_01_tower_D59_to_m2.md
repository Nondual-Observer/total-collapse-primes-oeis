# Sequence Card 01: (z + 2)(z^3 + 2z - 1) + z -> {-2} (D = -59)

* **Companion Article:** [`Total_Collapse_Primes_18_Families_Zhuk_2026.md`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.md) | [`PDF`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.pdf)
* **Verification Code:** [`../Code/`](../Code/)

## DATA (%S)

```text
3, 13, 137, 761, 2719, 150083
```

## NAME (%N)

```text
Primes p such that the iteration of T(z) = (z + 2)*(z^3 + 2*z - 1) + z mod p collapses every element of GF(p) to -2 (mod p).
```

## OFFSET (%O)

```text
1,1
```

## COMMENTS (%C)

```text
T(z) - z = (z + 2)*(z^3 + 2*z - 1) and T(z) + 2 = z*(z + 2)*(z^2 + 2). The cubic z^3 + 2*z - 1 has discriminant -59; its splitting field is the Hilbert class field of Q(sqrt(-59)) (class number 3).
For p = 13 = -P(-2), with P(z) = z^3 + 2*z - 1, the point -2 is a double fixed point. For p > 13, -2 is the only fixed point of T in GF(p) iff P is irreducible mod p, i.e., iff p is represented by the form 3*x^2 + x*y + 5*y^2; in particular every term p > 13 satisfies Kronecker(-59, p) = 1. Representations: 3 = (1,0), 137 = (-3,5), 761 = (12,7), 2719 = (29,4), 150083 = (153,112).
No other terms below 10^8 (exhaustive full functional-graph traversal of GF(p) for every prime p); a bounded-step orbit scan up to 2*10^9 found no further terms.
```

## LINKS (%H)

```text
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23420627">Unital Collapse Primes in Arithmetic Dynamics: Total Backward-Orbit Collapse for T(x) = (x^2 - 2)^2 mod p, Galois Stratification, and Fermat-Totient Asymptotics</a>, Zenodo preprint (2026).
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23856837">Total-Collapse Primes of Polynomial Maps over GF(p): Eighteen Families and Two Related Sequences</a>, Zenodo preprint (2026); <a href="https://github.com/Nondual-Observer/total-collapse-primes-oeis">GitHub verification suite</a>.
```

## FORMULA (%F)

```text
p is in the sequence iff #{z in {0..p-1} : T^(p)(z) == p - 2 (mod p)} = p, where T(z) = (z + 2)*(z^3 + 2*z - 1) + z mod p and T^(k) denotes the k-th iterate of T.
```

## EXAMPLE (%e)

```text
For p = 3: -2 == 1 (mod 3), and T(0) = 1, T(1) = 1, T(2) = 1 (mod 3), so all 3 elements of GF(3) reach -2 in 1 step; thus 3 is in the sequence.
For p = 13: -2 == 11 (mod 13), T(0) = 11, and every element of GF(13) reaches 11 within 6 steps under iteration of T(z) mod 13; thus 13 is in the sequence.
For p = 137: p = 3*(-3)^2 + (-3)*5 + 5*5^2 = 137, P_59(z) = z^3 + 2*z - 1 is irreducible mod 137, and all 137 elements of GF(137) collapse to -2 == 135 (mod 137); thus 137 is in the sequence.
```

## PROG (%o)

```text
(PARI)
isok(p) = {
  if (!isprime(p) || p == 2 || (p > 13 && kronecker(-59, p) != 1), return(0));
  my(tgt = (p - 2) % p, st = vector(p)); st[1] = st[tgt + 1] = 2;
  for (z = 1, p - 1,
    if (st[z+1] == 2, next);
    my(c = z, path = List());
    while (st[c+1] == 0, st[c+1] = 1; listput(path, c); c = ((c*(c + 2) % p)*(c^2 + 2) - 2) % p);
    if (st[c+1] == 1, return(0));
    for (i = 1, length(path), st[path[i]+1] = 2);
  );
  return(1);
}
(Julia)
isprime(n::Int) = n > 1 && all(d -> n % d != 0, 2:isqrt(n))
function isok(p::Int)
    (!isprime(p) || p == 2 || (p > 13 && powermod(-59, p >> 1, p) != 1)) && return false
    f(c) = mod(mod(c*(c + 2), p) * mod(c*c + 2, p) - 2, p)
    st = zeros(UInt8, p); st[1] = st[mod(p - 2, p) + 1] = 2
    for z in 1:p-1
        st[z + 1] == 2 && continue
        c = z; while st[c + 1] == 0; st[c + 1] = 1; c = f(c); end
        st[c + 1] == 1 && return false
        c = z; while st[c + 1] == 1; st[c + 1] = 2; c = f(c); end
    end
    return true
end
println([p for p in 2:3000 if isok(p)])
(Python)
from sympy import isprime
def isok(p: int) -> bool:
    if not isprime(p) or p == 2 or (p > 13 and pow(-59, p >> 1, p) != 1):
        return False
    state = bytearray(p)
    state[0] = state[(p - 2) % p] = 2
    for z in range(1, p):
        if state[z] == 2:
            continue
        curr, path = z, []
        while state[curr] == 0:
            state[curr] = 1
            path.append(curr)
            curr = (curr * (curr + 2) * (curr * curr + 2) - 2) % p
        if state[curr] == 1:
            return False
        for u in path:
            state[u] = 2
    return True
print([p for p in range(2, 3000) if isok(p)])
```

## CROSSREFS (%Y)

```text
Cf. A401173. See the linked article for the companion sequence for (z + 2)*P_251(z) + z.
```

## KEYWORDS (%K)

```text
nonn,hard,more
```

## AUTHOR (%A)

```text
_Igor Zhuk_, Nov 01 2026
```
