# Sequence Card 02: V_3(x)^2 -> {0}

* **Companion Article:** [`Total_Collapse_Primes_18_Families_Zhuk_2026.md`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.md) · [`PDF`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.pdf)
* **Verification Code:** [`../Code/`](../Code/)

## DATA (%S)

```text
2, 13, 23, 83, 17929, 18121
```

## NAME (%N)

```text
Primes p such that the iteration of T(x) = x^2*(x^2 - 3)^2 mod p collapses every element of GF(p) to 0.
```

## OFFSET (%O)

```text
1,1
```

## COMMENTS (%C)

```text
Here V_3(x) = x*(x^2 - 3) = 2*T_3(x/2) is the monic Chebyshev (Dickson) polynomial of degree 3 and T(x) = V_3(x)^2.
For odd p > 3, T(x) == 0 (mod p) iff x == 0 or x^2 == 3, so 0 has a nonzero preimage iff Kronecker(3, p) = 1; hence every odd term satisfies p == +-1 (mod 12). The prime 2 is a term, while 3 is not.
T(1) = 4 and T(4) = 2704 = 16*13^2.
No other terms below 10^8 (exhaustive full functional-graph traversal of GF(p) for every prime p); a bounded-step orbit scan up to 2*10^9 found no further terms.
```

## LINKS (%H)

```text
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23420627">Unital Collapse Primes in Arithmetic Dynamics: Total Backward-Orbit Collapse for T(x) = (x^2 - 2)^2 mod p, Galois Stratification, and Fermat-Totient Asymptotics</a>, Zenodo preprint (2026).
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23856837">Total-Collapse Primes of Polynomial Maps over GF(p): Eighteen Families and Two Related Sequences</a>, Zenodo preprint (2026); <a href="https://github.com/Nondual-Observer/total-collapse-primes-oeis">GitHub verification suite</a>.
```

## FORMULA (%F)

```text
p is in the sequence iff #{x in {0..p-1} : T^(p)(x) == 0 (mod p)} = p, where T(x) = x^2*(x^2 - 3)^2 mod p and T^(k) denotes the k-th iterate of T.
```

## EXAMPLE (%e)

```text
For p = 2: T(0) = 0, T(1) = 4 == 0 (mod 2), so all 2 elements of GF(2) reach 0 in 1 step; thus 2 is in the sequence.
For p = 13: T(0) = 0, and T(4) = 16*13^2 == 0 (mod 13), with all 13 elements of GF(13) reaching 0 within 4 steps; thus 13 is in the sequence.
For p = 23: all 23 elements of GF(23) reach 0 within 8 steps under iteration of x^2*(x^2 - 3)^2 mod 23; thus 23 is in the sequence.
```

## PROG (%o)

```text
(PARI)
isok(p) = {
  if (!isprime(p) || (p > 2 && p%12 != 1 && p%12 != 11), return(0));
  my(st = vector(p)); st[1] = 2;
  for (x = 1, p\2,
    if (st[x+1] == 2, next);
    my(c = x, path = List());
    while (st[c+1] == 0, st[c+1] = 1; listput(path, c); my(v3 = (c*(c^2 - 3)) % p); c = (v3^2) % p);
    if (st[c+1] == 1, return(0));
    for (i = 1, length(path), my(u = path[i]); st[u+1] = st[p - u + 1] = 2);
  );
  return(1);
}

(Julia)
isprime(n::Int) = n > 1 && all(d -> n % d != 0, 2:isqrt(n))
function isok(p::Int)
    (!isprime(p) || (p > 2 && !(p % 12 in (1, 11)))) && return false
    f(c) = mod(mod(c * mod(c*c - 3, p), p)^2, p)
    st = zeros(UInt8, p); st[1] = 2
    for x in 1:(p >> 1)
        st[x + 1] == 2 && continue
        c = x; while st[c + 1] == 0; st[c + 1] = 1; c = f(c); end
        st[c + 1] == 1 && return false
        c = x; while st[c + 1] == 1; st[c + 1] = st[p - c + 1] = 2; c = f(c); end
    end
    return true
end
println([p for p in 2:20000 if isok(p)])

(Python)
from sympy import isprime

def isok(p: int) -> bool:
    if not isprime(p) or (p > 2 and p % 12 not in (1, 11)):
        return False
    state = bytearray(p)
    state[0] = 2
    for x in range(1, (p >> 1) + 1):
        if state[x] == 2:
            continue
        curr, path = x, []
        while state[curr] == 0:
            state[curr] = 1
            path.append(curr)
            v3 = (curr * (curr * curr - 3)) % p
            curr = (v3 * v3) % p
        if state[curr] == 1:
            return False
        for u in path:
            state[u] = state[p - u] = 2
    return True

print([p for p in range(2, 20000) if isok(p)])
```

## CROSSREFS (%Y)

```text
Cf. A019434, A401173. See the linked article for the companion sequence for V_9(x)^2.
```

## KEYWORDS (%K)

```text
nonn,hard,more
```

## AUTHOR (%A)

```text
_Igor M. Zhuk_, Oct 17 2026
```
