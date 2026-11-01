# Sequence Card 19: Sum of n-bit primes with a single circular run of 1s and at most one internal block of 0s

* **Companion Article:** [`Total_Collapse_Primes_18_Families_Zhuk_2026.md`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.md) | [`PDF`](../Article/Total_Collapse_Primes_18_Families_Zhuk_2026.pdf)
* **Verification Code:** [`../Code/`](../Code/)

## DATA (%S)

```text
3, 12, 24, 119, 167, 657, 1895, 4563, 6687, 22176, 30913, 104134, 175320, 372369, 1203372, 2977563, 3323031, 14296358, 25906687, 43131484, 87119548, 205857709, 498514870, 988349366
```

## NAME (%N)

```text
Sum of n-bit primes whose binary expansion has a single circular run of 1s and at most one internal block of 0s.
```

## OFFSET (%O)

```text
2,1
```

## COMMENTS (%C)

```text
Equivalently, sum of primes of the form 2^n - 1 or (2^n - 1) - 2^s*(2^k - 1) with 1 <= k <= n - 2 and 1 <= s <= n - 1 - k (the odd n-bit prime terms of A100724): binary length n consisting of a block of 1s, an internal block of k zeros, and a block of s ones, together with 2^n - 1 itself.
For even n >= 4, 3 divides both 2^n - 1 and 2^k - 1 for even k, so every prime summand has odd zero-block length k; for odd n >= 5, 3 divides (2^n - 1) - 2^s*(2^k - 1) whenever k is odd and s is even.
The number of such primes of bit length n, for n >= 2, is 1, 2, 2, 5, 3, 7, 9, 11, 7, 14, 9, 16, 12, 13, 22, 27, 13, 32, 27, ...
For n = 7 the primes are 67, 71, 79, 97, 103, 113, 127 and a(7) = 657; for n = 10 they are 769, 911, 967, 991, 1009, 1019, 1021 and a(10) = 6687; for n = 16 there are 22 primes summing to a(16) = 1203372, whose two largest terms are the twin primes 65519 and 65521 = 2^16 - 15 straddling 65520 = lcm(1..10, 12..16).
```

## LINKS (%H)

```text
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23420627">Unital Collapse Primes in Arithmetic Dynamics: Total Backward-Orbit Collapse for T(x) = (x^2 - 2)^2 mod p, Galois Stratification, and Fermat-Totient Asymptotics</a>, Zenodo preprint (2026).
Igor M. Zhuk, <a href="https://doi.org/10.5281/zenodo.23856837">Total-Collapse Primes of Polynomial Maps over GF(p): Eighteen Families and Two Related Sequences</a>, Zenodo preprint (2026); <a href="https://github.com/Nondual-Observer/total-collapse-primes-oeis">GitHub verification suite</a>.
```

## FORMULA (%F)

```text
a(n) = [isprime(2^n - 1)]*(2^n - 1) + Sum_{k=1..n-2} Sum_{s=1..n-1-k} [isprime((2^n - 1) - 2^s*(2^k - 1))] * ((2^n - 1) - 2^s*(2^k - 1)).
```

## EXAMPLE (%e)

```text
For n = 2: the only candidate is 2^2 - 1 = 3 (binary 11_2), which is prime, so a(2) = 3.
For n = 3: the prime candidates are 5 (101_2) and 7 (111_2), so a(3) = 5 + 7 = 12.
For n = 4: the prime candidates are 11 (1011_2) and 13 (1101_2), so a(4) = 11 + 13 = 24.
For n = 5: the 5 prime candidates are 17, 19, 23, 29, 31, so a(5) = 17 + 19 + 23 + 29 + 31 = 119.
For n = 7: the 7 prime candidates are 67, 71, 79, 97, 103, 113, 127, so a(7) = 657.
For n = 10: the 7 prime candidates are 769, 911, 967, 991, 1009, 1019, 1021, so a(10) = 6687.
```

## PROG (%o)

```text
(PARI)
a(n) = {
  my(m = (1<<n) - 1, s = if(isprime(m), m, 0));
  for (k = 1, n - 2,
    my(h = (1<<k) - 1);
    for (sh = 1, n - 1 - k,
      my(p = m - (h << sh));
      if (isprime(p), s += p);
    );
  );
  return(s);
}
(Julia)
isprime(n::Int) = n > 1 && all(d -> n % d != 0, 2:isqrt(n))
function a(n::Int)
    m = (1 << n) - 1
    total = isprime(m) ? m : 0
    for k in 1:(n - 2)
        h = (1 << k) - 1
        for sh in 1:(n - 1 - k)
            p = m - (h << sh)
            isprime(p) && (total += p)
        end
    end
    return total
end
println([a(n) for n in 2:25])
(Python)
from sympy import isprime
def a(n: int) -> int:
    m = (1 << n) - 1
    total = m if isprime(m) else 0
    for k in range(1, n - 1):
        h = (1 << k) - 1
        for sh in range(1, n - k):
            p = m - (h << sh)
            if isprime(p):
                total += p
    return total
print([a(n) for n in range(2, 26)])
```

## CROSSREFS (%Y)

```text
Cf. A000668, A014234, A100724, A401173.
```

## KEYWORDS (%K)

```text
nonn,base,easy
```

## AUTHOR (%A)

```text
_Igor Zhuk_, Nov 01 2026
```
