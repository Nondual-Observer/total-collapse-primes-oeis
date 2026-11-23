#!/usr/bin/env julia
# ==============================================================================
# reference_implementation.jl
# Dependency-free Julia reference implementation for all 20 integer sequences
# in ../OEIS/: Cards 00-17 (total-collapse primes), Card 18 (least collapsing
# prime for K_n), and Card 19 (sums of n-bit primes with at most one internal
# block of 0s).
#
# Definition (companion article, Def. 1.1): p is a term iff every forward orbit
# of T on GF(p) enters the target set A. Traverses the full functional graph
# with no step cap.
#
# Run:  julia reference_implementation.jl      (self-test on all 20 cards)
# ==============================================================================

function isprime_int(n::Int)::Bool
    n < 2 && return false
    n % 2 == 0 && return n == 2
    d = 3
    while d * d <= n
        n % d == 0 && return false
        d += 2
    end
    return true
end

"""Monic Chebyshev-Dickson polynomial V_m(x) mod p."""
function chebyshev_V(m::Int, x::Int, p::Int)::Int
    a, b = mod(2, p), mod(x, p)
    m == 0 && return a
    for _ in 1:(m - 1)
        a, b = b, mod(widemul(x, b) - a, p)
    end
    return b
end

const P251_COEFFS = (1, 5, 6, -2, -4, 2, 9, -1)

function horner_mod(coeffs::Tuple, x::Int, p::Int)::Int
    acc = 0
    for c in coeffs
        acc = Int(mod(widemul(acc, x) + c, p))
    end
    return acc
end

const KN_MAP = Dict(11 => 2, 12 => 4, 13 => 1, 14 => 5, 15 => 6, 16 => 7, 17 => 8)

function eval_family(card::Int, x::Int, p::Int, invB::Int, A::Int)::Int
    if card in (0, 2, 4, 5, 6, 7, 9)
        m = card == 0 ? 2 : card == 2 ? 3 : card == 4 ? 5 : card == 5 ? 8 :
            card == 6 ? 9 : card == 7 ? 7 : 6
        v = chebyshev_V(m, x, p)
        return Int(mod(widemul(v, v), p))
    elseif card == 1
        z3 = Int(mod(widemul(Int(mod(widemul(x, x), p)), x), p))
        p59 = mod(z3 + 2 * x - 1, p)
        return Int(mod(widemul(mod(x + 2, p), p59) + x, p))
    elseif card == 3
        u = mod(Int(mod(widemul(x, x), p)) - 2, p)
        u2 = Int(mod(widemul(u, u), p))
        return Int(mod(widemul(u2, u2), p))
    elseif card == 8
        return Int(mod(widemul(mod(x + 2, p), horner_mod(P251_COEFFS, x, p)) + x, p))
    elseif card == 10
        y = Int(mod(widemul(x, mod(x - 2, p)), p))
        return Int(mod(widemul(y, y), p))
    else
        t = Int(mod(widemul(x, mod(Int(mod(widemul(x, x), p)) - A, p)), p))
        return Int(mod(widemul(t, invB), p))
    end
end

function cycle_ok(card::Int, cyc::Vector{Int}, p::Int)::Bool
    len = length(cyc)
    c0 = cyc[1]
    if card in (0, 3, 5)
        return len == 1 && c0 == mod(1, p)
    elseif card in (1, 8)
        return len == 1 && c0 == mod(-2, p)
    elseif card in (2, 6)
        return len == 1 && c0 == 0
    elseif card in (4, 7)
        return len == 1 && (c0 == 0 || c0 == mod(1, p))
    elseif card == 9
        return mod(4, p) in cyc
    elseif card == 10
        return len == 1
    else
        return (len == 1 && c0 == 0) || (len == 2 && Set(cyc) == Set((1, p - 1)))
    end
end

function is_term(card::Int, p::Int; n_override::Int = 0)::Bool
    !isprime_int(p) && return false
    invB, A = 0, 0
    if card >= 11
        nn = n_override > 0 ? n_override : KN_MAP[card]
        B = nn * (nn + 1)
        (p < 5 || B % p == 0) && return false
        invB = invmod(B, p)
        A = mod(nn * nn + nn + 1, p)
    end
    st = zeros(UInt8, p) # 0 unvisited, 1 on current path, 2 resolved
    ncyc = 0
    cyc = Int[]
    for x in 0:(p - 1)
        st[x + 1] != 0 && continue
        c = x
        while st[c + 1] == 0
            st[c + 1] = 1
            c = eval_family(card, c, p, invB, A)
        end
        if st[c + 1] == 1
            empty!(cyc)
            push!(cyc, c)
            y = eval_family(card, c, p, invB, A)
            while y != c
                push!(cyc, y)
                y = eval_family(card, y, p, invB, A)
            end
            ncyc += 1
            !cycle_ok(card, cyc, p) && return false
            (card == 9 && ncyc > 1) && return false
        end
        c = x
        while st[c + 1] == 1
            st[c + 1] = 2
            c = eval_family(card, c, p, invB, A)
        end
    end
    return true
end

function card18(n::Int; pmax::Int = 100_000)::Int
    d = n * n + (n + 1)^2
    r = isqrt(d)
    if r * r == d # Pell-Pythagorean obstruction (Article, Thm. 4.2)
        return (isprime_int(r) && r > n + 1 && is_term(11, r; n_override = n)) ? r : 0
    end
    for p in max(5, n + 2):pmax
        is_term(11, p; n_override = n) && return p
    end
    error("No collapsing prime found below pmax = $(pmax)")
end

function card19(n::Int)::Int
    cands = Set{Int}([ (1 << n) - 1 ])
    for k in 1:(n - 2), s in 1:(n - 1 - k)
        push!(cands, ((1 << n) - 1) - ((1 << k) - 1) << s)
    end
    return sum(q for q in cands if isprime_int(q); init = 0)
end

const DATA = Dict(
    0  => [3, 5, 37, 59, 61, 109, 359, 5507, 2117573, 503734381, 20900146669],
    1  => [3, 13, 137, 761, 2719, 150083],
    2  => [2, 13, 23, 83, 17929, 18121],
    3  => [3, 5, 13, 23, 41, 22381, 45503, 3852347, 1801981681],
    4  => [2, 3, 5, 31, 59, 61, 89, 2939, 4261, 11839, 14137859, 599999971],
    5  => [3, 5, 13, 17, 23, 47, 71, 89, 193, 313, 863, 1103, 1777, 107183, 84814153, 342037921],
    6  => [2, 11, 73, 1297, 6733, 7393, 658871, 1295447],
    7  => [2, 3, 5, 43, 421, 1429, 6257, 49279, 1580503, 24107329],
    8  => [3, 7, 17, 233, 1373, 1531, 5861, 27427, 52937, 125287],
    9  => [2, 3, 5, 7, 11, 13, 17, 37, 47, 79, 131, 167, 241, 257, 129401, 390209, 4926997, 6719077, 23734297],
    10 => [2, 3, 5, 7, 11, 23, 29, 31, 59, 181, 521, 1559, 2069, 5261, 6871, 9551, 19207, 1161371, 1574681, 1849721, 85268251],
    11 => [5, 7, 13, 19, 31, 47, 59, 73, 151, 331, 401, 1943069, 83584483, 158450797, 1554906533],
    12 => [7, 11, 17, 19, 29, 41, 53, 79, 191, 42979, 277829, 35718901],
    13 => [5, 7, 13, 23, 83, 97, 383, 14239193, 82411573],
    14 => [7, 17, 23, 37, 313, 8753, 72173, 2952797, 3516973, 8879231, 334677433],
    15 => [5, 11, 17, 31, 41, 109, 127, 199, 719, 2731, 26317, 307067, 617579, 4348559, 4838117, 30282823, 42166321, 62371087],
    16 => [5, 19, 23, 151, 1237, 11941, 761357, 6756377, 1105256599],
    17 => [5, 7, 19, 23, 31, 41, 53, 61, 277, 593, 809, 911, 1489, 2053, 4943, 6301, 17393, 229139, 324119, 1262783, 2055551, 228229363],
)

const DATA18 = [
    5, 5, 5, 7, 7, 11, 19, 19, 41, 13, 13, 17, 127, 19, 23, 19, 23, 31, 31, 0,
    23, 31, 53, 29, 37, 89, 89, 31, 41, 67, 37, 41, 53, 41, 53, 41, 47, 83, 211, 53,
    61, 59, 991, 47, 8713, 173, 199, 53, 389, 89, 67, 61, 73, 239, 109,
]

const DATA19 = [
    3, 12, 24, 119, 167, 657, 1895, 4563, 6687, 22176, 30913, 104134,
    175320, 372369, 1203372, 2977563, 3323031, 14296358, 25906687, 43131484,
    87119548, 205857709, 498514870, 988349366,
]

function self_test()
    lim = 20_000
    allok = true
    for card in 0:17
        got = [p for p in 2:(lim - 1) if is_term(card, p)]
        ref = [p for p in DATA[card] if p < lim]
        ok = (got == ref)
        println("Card $(lpad(card, 2, '0')): $(ok ? "OK " : "MISMATCH") $(got)")
        allok &= ok
    end
    g18 = [card18(n) for n in 1:55]
    println("Card 18: $(g18 == DATA18 ? "OK" : "MISMATCH")")
    allok &= (g18 == DATA18)
    g19 = [card19(n) for n in 2:25]
    println("Card 19: $(g19 == DATA19 ? "OK" : "MISMATCH")")
    allok &= (g19 == DATA19)
    println(allok ? "ALL OK" : "SOME MISMATCH")
    exit(allok ? 0 : 1)
end

if abspath(PROGRAM_FILE) == @__FILE__
    self_test()
end
