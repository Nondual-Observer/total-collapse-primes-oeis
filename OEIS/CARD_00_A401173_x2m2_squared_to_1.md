# Sequence Card 00: (x^2 - 2)^2 -> {1} (OEIS A401173)

* **OEIS:** [A401173](https://oeis.org/A401173)
* **Companion Article:** [`Total_Collapse_Primes_18_Families_Zhuk_2026.md`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.md) | [`PDF`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.pdf)
* **Verification Code:** [`../Code/`](../Code/)

## DATA (%S)

```text
3, 5, 37, 59, 61, 109, 359, 5507, 2117573, 503734381, 20900146669
```

## NAME (%N)

```text
Primes p such that the iteration of T(x) = (x^2 - 2)^2 mod p collapses every element of GF(p) to 1.
```

## OFFSET (%O)

```text
1,1
```

## COMMENTS (%C)

```text
Equivalently, primes p such that the backward orbit (preimage tree) of 1 under T has size p, i.e., the functional graph of T on GF(p) is connected and its only cycle is the fixed point 1.
T(x) - x = (x - 1)*(x^3 + x^2 - 3*x - 4), and the cubic factor has discriminant -83; its splitting field is the Hilbert class field of Q(sqrt(-83)) (class number 3), and x -> x^2 - 3 maps its roots to the roots of z^3 + 2*z^2 + 2*z - 1. Hence for p > 5, 1 is the only fixed point of T in GF(p) iff p is represented by the non-principal form 3*x^2 + x*y + 7*y^2 (12*p = (6*x + y)^2 + 83*y^2). Every term other than 5 has this form: 3 = (1,0), 37 = (3,1), 59 = (4,1), 61 = (3,2), 109 = (-6,1), 359 = (-11,1), 5507 = (-43,5), 2117573 = (209,518), 503734381 = (-13014,1421), 20900146669 = (24627,-53998).
If p == 7 (mod 12), then T^(-1)(1) = {1, -1} and -1 has no preimage, so p is not a term.
The primitive period-2 factor (T(T(x)) - x)/(T(x) - x) has degree 12 and discriminant -3*641*857^2*3581^2.
No other terms below 3*10^10 (verified by an exact unbounded-Brent and full functional-graph search).
```

## LINKS (%H)

```text
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23420627">Unital Collapse Primes in Arithmetic Dynamics: Total Backward-Orbit Collapse for T(x) = (x^2 - 2)^2 mod p, Galois Stratification, and Fermat-Totient Asymptotics</a>, Zenodo preprint (2026).
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23856837">Total-Collapse Primes of Polynomial Maps over GF(p): Eighteen Families and Two Related Sequences</a>, Zenodo preprint (2026); <a href="https://github.com/Nondual-Observer/total-collapse-primes-oeis">GitHub verification suite</a>.
```

## FORMULA (%F)

```text
p is in the sequence iff #{x in {0..p-1} : T^(p)(x) == 1 (mod p)} = p, where T(x) = (x^2 - 2)^2 mod p and T^(k) denotes the k-th iterate of T.
```

## EXAMPLE (%e)

```text
For p = 3: T(0) = 4 == 1 (mod 3), T(1) = T(2) = 1 (mod 3), so all 3 elements of GF(3) reach 1 in 1 step; thus 3 is in the sequence.
For p = 5: T(0) = 4, T(4) = 1, T(1) = 1, T(2) = T(3) = 4 -> 1 (mod 5), so all 5 elements of GF(5) reach 1 within 2 steps; thus 5 is in the sequence.
For p = 37: p = 3*3^2 + 3*1 + 7*1^2 = 37, Phi_1*(x) = x^3 + x^2 - 3*x - 4 is irreducible mod 37, and all 37 elements of GF(37) collapse to 1 within 7 steps; thus 37 is in the sequence.
```

## PROG (%o)

```text
(PARI)
isok(p) = {
  if (!isprime(p) || p == 2 || (p > 5 && (p%12 == 7 || kronecker(-83, p) != 1)), return(0));
  my(st = vector(p)); st[2] = st[p] = 2;
  for (x = 0, p\2,
    if (st[x+1] == 2, next);
    my(c = x, path = List());
    while (st[c+1] == 0, st[c+1] = 1; listput(path, c); c = ((c^2 - 2) % p)^2 % p);
    if (st[c+1] == 1, return(0));
    for (i = 1, length(path), my(u = path[i]); st[u+1] = 2; if (u, st[p - u + 1] = 2));
  );
  return(1);
}
(Julia)
isprime(n::Int) = n > 1 && all(d -> n % d != 0, 2:isqrt(n))
function isok(p::Int)
    (!isprime(p) || p == 2 || (p > 5 && (p % 12 == 7 || powermod(-83, p >> 1, p) != 1))) && return false
    f(c) = begin
        y = mod(widemul(c, c) - 2, p)
        Int(mod(y*y, p))
    end
    st = zeros(UInt8, p); st[2] = st[p] = 2
    for x in 0:(p >> 1)
        st[x + 1] == 2 && continue
        c = x; while st[c + 1] == 0; st[c + 1] = 1; c = f(c); end
        st[c + 1] == 1 && return false
        c = x; while st[c + 1] == 1; st[c + 1] = 2; c > 0 && (st[p - c + 1] = 2); c = f(c); end
    end
    return true
end
println([p for p in 2:6000 if isok(p)])
(Python)
from sympy import isprime
def isok(p: int) -> bool:
    if not isprime(p) or p == 2:
        return False
    if p > 5 and (p % 12 == 7 or pow(-83, p >> 1, p) != 1):
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
            curr = pow(curr * curr - 2, 2, p)
        if state[curr] == 1:
            return False
        for u in path:
            state[u] = state[-u % p] = 2
    return True
print([p for p in range(2, 6000) if isok(p)])
```

## CROSSREFS (%Y)

```text
Cf. A019434. See the linked article for the companion sequences for (x^2 - 2)^4, V_8(x)^2, (z + 2)*(z^3 + 2*z - 1) + z and (u^3 - 7*u)/6.
```

## KEYWORDS (%K)

```text
nonn,hard,more
```

## AUTHOR (%A)

```text
_Igor Zhuk_, Oct 16 2026
```
