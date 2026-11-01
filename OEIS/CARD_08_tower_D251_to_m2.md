# Sequence Card 08: (z + 2) P_251(z) + z -> {-2} (D = -251)

* **Companion Article:** [`Total_Collapse_Primes_18_Families_Zhuk_2026.md`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.md) · [`PDF`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.pdf)
* **Verification Code:** [`../Code/`](../Code/)

## DATA (%S)

```text
3, 7, 17, 233, 1373, 1531, 5861, 27427, 52937, 125287
```

## NAME (%N)

```text
Primes p such that the iteration of T(z) = (z + 2)*(z^7 + 5*z^6 + 6*z^5 - 2*z^4 - 4*z^3 + 2*z^2 + 9*z - 1) + z mod p collapses every element of GF(p) to -2 (mod p).
```

## OFFSET (%O)

```text
1,1
```

## COMMENTS (%C)

```text
P(z) = z^7 + 5*z^6 + 6*z^5 - 2*z^4 - 4*z^3 + 2*z^2 + 9*z - 1 has discriminant -2^10*251^3, P(-2) = -11, and T(z) - z = (z + 2)*P(z). The splitting field of P is the Hilbert class field H of Q(sqrt(-251)) (class number 7), with Gal(H/Q) = D_7.
If Kronecker(-251, p) = -1 (p != 2), Frobenius acts as a reflection, so P has exactly one root mod p; it is a fixed point of T different from -2 when p != 11. Hence every term p > 11 satisfies Kronecker(-251, p) = 1, and then -2 is the only fixed point iff p is represented by one of the non-principal forms 3*x^2 + x*y + 21*y^2 (p = 3, 1373, 52937), 7*x^2 + x*y + 9*y^2 (p = 7, 17, 1531, 5861) or 5*x^2 + 3*x*y + 13*y^2 (p = 233, 27427, 125287).
No other terms below 10^6 (exhaustive full functional-graph traversal of GF(p) for every prime p); a bounded-step orbit scan up to 5*10^8 found no further terms.
```

## LINKS (%H)

```text
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23420627">Unital Collapse Primes in Arithmetic Dynamics: Total Backward-Orbit Collapse for T(x) = (x^2 - 2)^2 mod p, Galois Stratification, and Fermat-Totient Asymptotics</a>, Zenodo preprint (2026).
Igor M. Zhuk, <a href="https://github.com/Nondual-Observer/total-collapse-primes-oeis/blob/main/Article/Total_Collapse_Primes_18_Families_Zhuk_2026.pdf">Total-Collapse Primes of Polynomial Maps over GF(p): Eighteen Families and Two Related Sequences</a>, companion preprint and verification suite (2026).
```

## FORMULA (%F)

```text
p is in the sequence iff #{z in {0..p-1} : T^(p)(z) == p - 2 (mod p)} = p, where T(z) = (z + 2)*(z^7 + 5*z^6 + 6*z^5 - 2*z^4 - 4*z^3 + 2*z^2 + 9*z - 1) + z mod p and T^(k) denotes the k-th iterate of T.
```

## EXAMPLE (%e)

```text
For p = 3: -2 == 1 (mod 3), and every element of GF(3) reaches 1 within 2 steps under iteration of T(z) mod 3; thus 3 is in the sequence.
For p = 233: p = 5*4^2 + 3*4*3 + 13*3^2 = 233 in Cl(-251), P_251(z) is irreducible mod 233, and all 233 elements of GF(233) collapse to -2 == 231 (mod 233) within 27 steps; thus 233 is in the sequence.
```

## PROG (%o)

```text
(PARI)
isok(p) = {
  if (!isprime(p) || p == 2 || (p > 17 && kronecker(-251, p) != 1), return(0));
  my(tgt = (p - 2) % p, st = vector(p)); st[1] = st[tgt + 1] = 2;
  for (z = 1, p - 1,
    if (st[z+1] == 2, next);
    my(c = z, path = List());
    while (st[c+1] == 0, st[c+1] = 1; listput(path, c); c = ((c + 2)*(c^7 + 5*c^6 + 6*c^5 - 2*c^4 - 4*c^3 + 2*c^2 + 9*c - 1) + c) % p);
    if (st[c+1] == 1, return(0));
    for (i = 1, length(path), st[path[i]+1] = 2);
  );
  return(1);
}

(Julia)
isprime(n::Int) = n > 1 && all(d -> n % d != 0, 2:isqrt(n))
function isok(p::Int)
    (!isprime(p) || p == 2 || (p > 17 && powermod(-251, p >> 1, p) != 1)) && return false
    function f(c)
        s = 1
        for a in (5, 6, -2, -4, 2, 9, -1); s = mod(s*c + a, p); end
        return mod((c + 2)*s + c, p)
    end
    st = zeros(UInt8, p); st[1] = st[mod(p - 2, p) + 1] = 2
    for z in 1:p-1
        st[z + 1] == 2 && continue
        c = z; while st[c + 1] == 0; st[c + 1] = 1; c = f(c); end
        st[c + 1] == 1 && return false
        c = z; while st[c + 1] == 1; st[c + 1] = 2; c = f(c); end
    end
    return true
end
println([p for p in 3:6000 if isok(p)])

(Python)
from sympy import isprime

def isok(p: int) -> bool:
    if not isprime(p) or p == 2 or (p > 17 and pow(-251, p >> 1, p) != 1):
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
            p251 = (((((((curr + 5)*curr + 6)*curr - 2)*curr - 4)*curr + 2)*curr + 9)*curr - 1) % p
            curr = ((curr + 2) * p251 + curr) % p
        if state[curr] == 1:
            return False
        for u in path:
            state[u] = 2
    return True

print([p for p in range(2, 6000) if isok(p)])
```

## CROSSREFS (%Y)

```text
Cf. A401173. See the linked article for the companion sequence for (z + 2)*(z^3 + 2*z - 1) + z.
```

## KEYWORDS (%K)

```text
nonn,hard,more
```

## AUTHOR (%A)

```text
_Igor M. Zhuk_, Oct 17 2026
```
