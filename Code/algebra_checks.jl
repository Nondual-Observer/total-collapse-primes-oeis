#!/usr/bin/env julia
# ==============================================================================
# algebra_checks.jl
# Exact symbolic and arithmetic verification of all polynomial identities,
# discriminants, resultants, quadratic-form representations, D_7 Frobenius
# splitting types, Pell-Pythagorean obstructions, and projective/spectral
# identities in the companion article (Sections 3-4, Appendices A-B).
#
# Zero external dependencies (uses only Julia Base + LinearAlgebra).
# Run:  julia algebra_checks.jl     (every line starts with PASS)
# ==============================================================================

using LinearAlgebra

# --- Exact univariate polynomial arithmetic over BigInt: coeffs [c_0, c_1, ..., c_d] ---
function ptrim(a::Vector{BigInt})::Vector{BigInt}
    i = length(a)
    while i > 1 && iszero(a[i])
        i -= 1
    end
    return isempty(a) ? BigInt[0] : a[1:i]
end

pdeg(a::Vector{BigInt})::Int = (length(a) == 1 && iszero(a[1])) ? -1 : length(a) - 1

function padd(a::Vector{BigInt}, b::Vector{BigInt})::Vector{BigInt}
    m = max(length(a), length(b))
    c = zeros(BigInt, m)
    for i in eachindex(a); c[i] += a[i]; end
    for i in eachindex(b); c[i] += b[i]; end
    return ptrim(c)
end

psub(a::Vector{BigInt}, b::Vector{BigInt})::Vector{BigInt} = padd(a, [-x for x in b])

function pmul(a::Vector{BigInt}, b::Vector{BigInt})::Vector{BigInt}
    c = zeros(BigInt, length(a) + length(b) - 1)
    for i in eachindex(a), j in eachindex(b)
        c[i + j - 1] += a[i] * b[j]
    end
    return ptrim(c)
end

function ppow(a::Vector{BigInt}, e::Int)::Vector{BigInt}
    r = BigInt[1]
    for _ in 1:e
        r = pmul(r, a)
    end
    return r
end

function pcomp(f::Vector{BigInt}, g::Vector{BigInt})::Vector{BigInt}
    acc = BigInt[0]
    for i in length(f):-1:1
        acc = padd(pmul(acc, g), BigInt[f[i]])
    end
    return ptrim(acc)
end

function pderiv(a::Vector{BigInt})::Vector{BigInt}
    length(a) <= 1 && return BigInt[0]
    return ptrim([BigInt(i - 1) * a[i] for i in 2:length(a)])
end

function peval(a::Vector{BigInt}, x::Integer)::BigInt
    acc = BigInt(0)
    bx = BigInt(x)
    for i in length(a):-1:1
        acc = acc * bx + a[i]
    end
    return acc
end

function pdivrem_exact(a::Vector{BigInt}, b::Vector{BigInt})
    r = copy(a)
    db = pdeg(b)
    lead = b[end]
    q = zeros(BigInt, max(1, length(a) - length(b) + 1))
    while pdeg(r) >= db && !(length(r) == 1 && iszero(r[1]))
        d = pdeg(r) - db
        q_coeff, rem_c = divrem(r[end], lead)
        @assert iszero(rem_c) "Non-exact integer polynomial division"
        q[d + 1] = q_coeff
        for j in 0:db
            r[d + j + 1] -= q_coeff * b[j + 1]
        end
        r = ptrim(r)
    end
    return ptrim(q), ptrim(r)
end

"""Exact Bareiss fraction-free determinant over BigInt."""
function det_bareiss(M_in::Matrix{BigInt})::BigInt
    n = size(M_in, 1)
    n == 0 && return BigInt(1)
    M = deepcopy(M_in)
    sign = 1
    prev_pivot = BigInt(1)
    for k in 1:(n - 1)
        if iszero(M[k, k])
            piv = findfirst(r -> !iszero(M[r, k]), (k + 1):n)
            piv === nothing && return BigInt(0)
            r = k + piv
            M[k, :], M[r, :] = M[r, :], M[k, :]
            sign = -sign
        end
        for i in (k + 1):n
            for j in (k + 1):n
                M[i, j] = div(M[i, j] * M[k, k] - M[i, k] * M[k, j], prev_pivot)
            end
            M[i, k] = 0
        end
        prev_pivot = M[k, k]
    end
    return sign * M[n, n]
end

"""Sylvester resultant Res_x(f, g) over Z."""
function resultant_z(f::Vector{BigInt}, g::Vector{BigInt})::BigInt
    m, n = pdeg(f), pdeg(g)
    m == 0 && return f[1]^n
    n == 0 && return g[1]^m
    N = m + n
    S = fill(BigInt(0), N, N)
    # f_m x^m + ... + f_0 in n rows
    for row in 1:n
        for j in 0:m
            S[row, row + j] = f[m + 1 - j]
        end
    end
    # g_n x^n + ... + g_0 in m rows
    for row in 1:m
        for j in 0:n
            S[n + row, row + j] = g[n + 1 - j]
        end
    end
    return det_bareiss(S)
end

function discriminant_z(f::Vector{BigInt})::BigInt
    d = pdeg(f)
    sgn = ((d * (d - 1)) >> 1) % 2 == 1 ? -1 : 1
    return div(sgn * resultant_z(f, pderiv(f)), f[end])
end

# --- Polynomial arithmetic over GF(p) for distinct-degree factorization ---
function fp_trim(a::Vector{Int}, p::Int)::Vector{Int}
    b = [mod(c, p) for c in a]
    i = length(b)
    while i > 1 && b[i] == 0; i -= 1; end
    return isempty(b) ? [0] : b[1:i]
end

fp_deg(a::Vector{Int})::Int = (length(a) == 1 && a[1] == 0) ? -1 : length(a) - 1

function fp_mul(a::Vector{Int}, b::Vector{Int}, p::Int)::Vector{Int}
    c = zeros(Int, length(a) + length(b) - 1)
    for i in eachindex(a), j in eachindex(b)
        c[i + j - 1] = mod(c[i + j - 1] + a[i] * b[j], p)
    end
    return fp_trim(c, p)
end

function fp_rem(a::Vector{Int}, b::Vector{Int}, p::Int)::Vector{Int}
    r = fp_trim(a, p)
    db = fp_deg(b)
    inv_lead = invmod(b[end], p)
    while fp_deg(r) >= db && fp_deg(r) >= 0
        d = fp_deg(r) - db
        coef = mod(r[end] * inv_lead, p)
        for j in 0:db
            r[d + j + 1] = mod(r[d + j + 1] - coef * b[j + 1], p)
        end
        r = fp_trim(r, p)
    end
    return r
end

function fp_gcd(a::Vector{Int}, b::Vector{Int}, p::Int)::Vector{Int}
    x, y = fp_trim(a, p), fp_trim(b, p)
    while fp_deg(y) >= 0
        x, y = y, fp_rem(x, y, p)
    end
    if fp_deg(x) >= 0
        inv_l = invmod(x[end], p)
        x = [mod(c * inv_l, p) for c in x]
    end
    return x
end

function fp_powmod(base::Vector{Int}, exp::Int, modp::Vector{Int}, p::Int)::Vector{Int}
    res = [1]
    b = fp_rem(base, modp, p)
    e = exp
    while e > 0
        if (e & 1) == 1
            res = fp_rem(fp_mul(res, b, p), modp, p)
        end
        b = fp_rem(fp_mul(b, b, p), modp, p)
        e >>= 1
    end
    return res
end

"""Degree partition of a squarefree polynomial over GF(p) via distinct-degree factorization."""
function ddf_degrees(poly_z::Vector{BigInt}, p::Int)::Vector{Int}
    f = fp_trim([Int(mod(c, p)) for c in poly_z], p)
    degs = Int[]
    h = [0, 1] # x
    d = 1
    while fp_deg(f) > 0
        h = fp_powmod(h, p, f, p) # x^(p^d) mod f
        hx = copy(h)
        if length(hx) < 2; push!(hx, 0); end
        hx[2] = mod(hx[2] - 1, p)
        g = fp_gcd(fp_trim(hx, p), f, p)
        dg = fp_deg(g)
        if dg > 0
            cnt = div(dg, d)
            for _ in 1:cnt; push!(degs, d); end
            # divide f by g over GF(p)
            while fp_deg(fp_rem(f, g, p)) == -1 && fp_deg(f) >= dg && dg > 0
                # single exact division
                q = zeros(Int, length(f) - length(g) + 1)
                r = copy(f)
                inv_l = invmod(g[end], p)
                while fp_deg(r) >= dg
                    sh = fp_deg(r) - dg
                    c = mod(r[end] * inv_l, p)
                    q[sh + 1] = c
                    for j in 0:dg
                        r[sh + j + 1] = mod(r[sh + j + 1] - c * g[j + 1], p)
                    end
                    r = fp_trim(r, p)
                end
                f = fp_trim(q, p)
                break
            end
        end
        d += 1
    end
    return sort(degs)
end

function isprime_int(n::Integer)::Bool
    n < 2 && return false
    n % 2 == 0 && return n == 2
    d = oftype(n, 3)
    while d * d <= n
        n % d == 0 && return false
        d += 2
    end
    return true
end

function kronecker_odd_p(a::Int, p::Int)::Int
    m = mod(a, p)
    m == 0 && return 0
    r = powermod(m, (p - 1) >> 1, p)
    return r == 1 ? 1 : -1
end

function rep_form(p::Int, a::Int, b::Int, c::Int)::Bool
    ymax = isqrt(div(4 * a * p, 4 * a * c - b * b)) + 2
    xmax = isqrt(p) + 3
    for Y in 0:ymax, X in (-xmax):xmax
        a * X * X + b * X * Y + c * Y * Y == p && return true
    end
    return false
end

function all_perms(v::Vector{Int})
    length(v) == 1 && return [copy(v)]
    out = Vector{Vector{Int}}()
    for i in eachindex(v)
        rest = [v[j] for j in eachindex(v) if j != i]
        for p in all_perms(rest)
            push!(out, vcat(v[i], p))
        end
    end
    return out
end

function run_checks()
    failures = 0
    function ok(cond::Bool, msg::String)
        if cond
            println("PASS " * msg)
        else
            failures += 1
            println("FAIL " * msg)
        end
    end

    # --- Section 3.1 (Theorem 3.1 & Corollary 3.2): odd cubic ladder K_n ---
    # Polynomial identity of total degree <= 9 in u and <= 10 in n holds in Z[n, u]
    # iff it holds for 12 distinct integer values n = 1..12.
    id1_ok, id2_ok, id3_ok, id4_ok, id5_ok = true, true, true, true, true
    X = BigInt[0, 1]
    for nv in 1:12
        B = BigInt(nv) * (nv + 1)
        A = BigInt(nv)^2 + nv + 1
        Dfix = 2 * B + 1
        numK = BigInt[0, -A, 0, 1] # u^3 - A u
        id1_ok &= (psub(numK, BigInt[B]) == pmul(pmul(BigInt[1, 1], BigInt[nv, 1]), BigInt[-(nv + 1), 1]))
        id2_ok &= (padd(numK, BigInt[B]) == pmul(pmul(BigInt[-1, 1], BigInt[-nv, 1]), BigInt[nv + 1, 1]))
        # B^4*(K(K(u)) - u) = (numK(u))^3 - A*B^2*numK(u) - B^4*u
        lhs = psub(psub(ppow(numK, 3), pmul(BigInt[A * B^2], numK)), pmul(BigInt[B^4], X))
        rhs = pmul(pmul(pmul(X, BigInt[-1, 0, 1]), BigInt[-Dfix, 0, 1]), BigInt[B^2, 0, -A, 0, 1])
        id3_ok &= (lhs == rhs)
        id4_ok &= (A^2 - 4 * B^2 == (3 * B + 1) * (1 - B))
        id5_ok &= (3 * B + 1 == BigInt(nv + 1)^3 - BigInt(nv)^3)
    end
    ok(id1_ok, "Thm 3.1: u^3 - A_n*u - B_n = (u+1)(u+n)(u-n-1)")
    ok(id2_ok, "Thm 3.1: u^3 - A_n*u + B_n = (u-1)(u-n)(u+n+1)")
    ok(id3_ok, "Thm 3.1: B_n^4*(K_n(K_n(u)) - u) = u(u^2-1)(u^2-D_fix)(u^4-A_n*u^2+B_n^2)")
    ok(id4_ok, "Thm 3.1: A_n^2 - 4*B_n^2 = (3*B_n + 1)(1 - B_n)")
    ok(id5_ok, "Thm 3.1: 3*B_n + 1 = (n+1)^3 - n^3")

    dfix_oct = [k * k + (k + 1)^2 for k in 1:8]
    ok(dfix_oct == [5, 13, 25, 41, 61, 85, 113, 145],
       "Cor 3.2: D_fix(1..8) = [5, 13, 25, 41, 61, 85, 113, 145] (square only at n = 3)")

    # --- Section 3.2 (Theorem 3.3 & Corollary 3.4): Cards 00, 01, 08 ---
    T2 = ppow(BigInt[-2, 0, 1], 2)
    f83 = BigInt[-4, -3, 1, 1]
    P83 = BigInt[-1, 2, 2, 1]
    ok(psub(T2, X) == pmul(BigInt[-1, 1], f83), "Thm 3.3(1): (x^2-2)^2 - x = (x-1)(x^3+x^2-3x-4)")
    ok(discriminant_z(f83) == -83 && peval(f83, 1) == -5, "Thm 3.3(1): disc(f) = -83, f(1) = -5")
    f83_mod83 = fp_trim([Int(mod(c, 83)) for c in psub(pmul(ppow(BigInt[-44, 1], 2), BigInt[-77, 1]), f83)], 83)
    ok(f83_mod83 == [0], "Thm 3.3(1): f(x) == (x-44)^2*(x-77) (mod 83)")

    # Resultant_x(f83(x), z - (x^2 - 3)) == P83(z) tested on 5 interpolation nodes z = -2..2
    res_pts_ok = all(resultant_z(f83, BigInt[zv + 3, 0, -1]) == peval(P83, zv) for zv in -2:2)
    ok(res_pts_ok && discriminant_z(P83) == -83,
       "Cor 3.4: z = x^2 - 3 maps roots of f(x) to z^3 + 2z^2 + 2z - 1 (disc = -83)")

    TT2 = pcomp(T2, T2)
    phi2, rem2 = pdivrem_exact(psub(TT2, X), psub(T2, X))
    ok(rem2 == BigInt[0] && pdeg(phi2) == 12 &&
       discriminant_z(phi2) == -BigInt(3) * 641 * BigInt(857)^2 * BigInt(3581)^2,
       "Card 00: deg((T(T(x))-x)/(T(x)-x)) = 12, disc = -3*641*857^2*3581^2")

    reps83 = [
        (BigInt(3), 1, 0), (BigInt(37), 3, 1), (BigInt(59), 4, 1), (BigInt(61), 3, 2),
        (BigInt(109), -6, 1), (BigInt(359), -11, 1), (BigInt(5507), -43, 5),
        (BigInt(2117573), 209, 518), (BigInt(503734381), -13014, 1421),
        (BigInt(20900146669), 24627, -53998),
    ]
    ok(all(isprime_int(p) && 3 * BigInt(a)^2 + BigInt(a) * b + 7 * BigInt(b)^2 == p for (p, a, b) in reps83),
       "Card 00: primality and 3x^2 + xy + 7y^2 representations of all listed terms p != 5")

    P59 = BigInt[-1, 2, 0, 1]
    T59 = padd(pmul(BigInt[2, 1], P59), X)
    ok(discriminant_z(P59) == -59 && peval(P59, -2) == -13,
       "Thm 3.3(2): disc(z^3+2z-1) = -59, P_59(-2) = -13")
    ok(padd(T59, BigInt[2]) == pmul(pmul(X, BigInt[2, 1]), BigInt[2, 0, 1]),
       "Card 01: T(z) + 2 = z(z+2)(z^2+2)")
    p59_mod59 = fp_trim([Int(mod(c, 59)) for c in psub(pmul(ppow(BigInt[-45, 1], 2), BigInt[-28, 1]), P59)], 59)
    ok(p59_mod59 == [0], "Thm 3.3(2): P_59(z) == (z-45)^2*(z-28) (mod 59)")

    reps59 = [(3, 1, 0), (137, -3, 5), (761, 12, 7), (2719, 29, 4), (150083, 153, 112)]
    ok(all(isprime_int(p) && 3 * a * a + a * b + 5 * b * b == p for (p, a, b) in reps59),
       "Card 01: primality and 3x^2 + xy + 5y^2 representations of all terms p != 13")

    P251 = BigInt[-1, 9, 2, -4, -2, 6, 5, 1]
    ok(discriminant_z(P251) == -BigInt(2)^10 * BigInt(251)^3 && peval(P251, -2) == -11,
       "Thm 3.3(3): disc(P_251) = -2^10*251^3, P_251(-2) = -11")

    bad_d7 = Int[]
    for p in 3:2999
        (!isprime_int(p) || p == 251) && continue
        degs = ddf_degrees(P251, p)
        ks = kronecker_odd_p(-251, p)
        if (ks == -1 && degs != [1, 2, 2, 2]) || (ks == 1 && degs != [1, 1, 1, 1, 1, 1, 1] && degs != [7])
            push!(bad_d7, p)
        end
    end
    ok(isempty(bad_d7), "Thm 3.3(3): P_251 splitting types match Gal = D_7 across all primes p < 3000 (also proving Q-irreducibility)")

    bad_f83 = [p for p in 7:2999 if isprime_int(p) &&
               (all(mod(r^3 + r^2 - 3 * r - 4, p) != 0 for r in 0:(p - 1)) != rep_form(p, 3, 1, 7))]
    ok(isempty(bad_f83), "Cor 3.4: f(x) irreducible mod p <=> p = 3x^2 + xy + 7y^2 for all 7 <= p < 3000")

    # --- Section 3.3 (Propositions 3.7 & 3.8): Cards 03, 09, 10 ---
    u7 = BigInt[-2, 1, -1, 1, -1, 1, -1, 1]
    ok(pmul(BigInt[1, 1], u7) == BigInt[-2, -1, 0, 0, 0, 0, 0, 0, 1] && discriminant_z(u7) == -26522311,
       "Prop 3.7 (Card 03): (u^8 - u - 2)/(u + 1) has discriminant -26522311")

    V3 = BigInt[0, -3, 0, 1]
    V6 = psub(ppow(V3, 2), BigInt[2])
    T6 = ppow(V6, 2)
    ok(all(peval(T6, v) == 4 for v in (0, 1, -1, 2, -2)) &&
       peval(T6, 4) - 4 == BigInt(2)^6 * BigInt(3)^3 * BigInt(5)^2 * BigInt(13)^2,
       "Prop 3.8(1) (Card 09): V_6(x)^2 = 4 on {0, +-1, +-2} and V_6(4)^2 - 4 = 2^6*3^3*5^2*13^2")

    Dmap = ppow(pmul(X, BigInt[-2, 1]), 2)
    ok(psub(Dmap, X) == pmul(pmul(X, BigInt[-1, 1]), BigInt[1, -3, 1]) &&
       pcomp(Dmap, BigInt[2, -1]) == Dmap,
       "Prop 3.8(2) (Card 10): D(x) - x = x(x-1)(x^2-3x+1) and D(2-x) = D(x)")

    # --- Section 4.1 (Theorem 4.2): Pell-Pythagorean obstruction (Card 18) ---
    nk = BigInt[0, 3]
    while length(nk) < 10
        push!(nk, 6 * nk[end] - nk[end - 1] + 2)
    end
    mk = [isqrt(k * k + (k + 1)^2) for k in nk[2:end]]
    ok(all(k * k + (k + 1)^2 == c * c for (k, c) in zip(nk[2:end], mk)),
       "Thm 4.2: first 9 Pell legs n_k in A001652 satisfy n_k^2 + (n_k+1)^2 = c_k^2")
    g_pell(v, m) = mod(v * (1 - 2 * v * v), m)
    cyc29 = [22]; for _ in 1:4; push!(cyc29, g_pell(cyc29[end], 29)); end
    cyc5741 = [1762]; for _ in 1:4; push!(cyc5741, g_pell(cyc5741[end], 5741)); end
    ok(cyc29 == [22, 12, 7, 17, 22] && length(Set(cyc29[1:4])) == 4,
       "Thm 4.2: extraneous 4-cycle 22 -> 12 -> 7 -> 17 -> 22 for n = 20, p = 29")
    ok(cyc5741 == [1762, 3531, 3979, 2210, 1762] && length(Set(cyc5741[1:4])) == 4,
       "Thm 4.2: extraneous 4-cycle 1762 -> 3531 -> 3979 -> 2210 -> 1762 for n = 4059, p = 5741")
    ok(195025 == 5^2 * 29 * 269 && 1136689 == 137 * 8297 && 6625109 == 37 * 179057 &&
       isprime_int(269) && isprime_int(137) && isprime_int(8297) && isprime_int(37) && isprime_int(179057),
       "Thm 4.2: composite factorizations of c_7 = 195025, c_8 = 1136689, c_9 = 6625109")
    ok(all(mod(v * (v * v - (119^2 + 119 + 1)) * invmod(119 * 120, 13) -
               v * (v * v - 7) * invmod(6, 13), 13) == 0 for v in 0:12),
       "Sec 4.1: modular aliasing K_119 == K_2 (mod 13)")

    # --- Appendices A & B: projective combinatorics and auxiliary identities ---
    lines_fano = Set([
        Set([0, 1, 3]), Set([1, 2, 4]), Set([2, 3, 5]), Set([3, 4, 6]),
        Set([4, 5, 0]), Set([5, 6, 1]), Set([6, 0, 2]),
    ])
    G = Vector{Vector{Int}}()
    for perm in all_perms(collect(0:6))
        if all(Set([perm[i + 1] for i in L]) in lines_fano for L in lines_fano)
            push!(G, perm)
        end
    end
    orb_by_w = Dict(w => Int[] for w in 1:6)
    seen = Set{UInt8}()
    for mask in UInt8(1):UInt8(126)
        mask in seen && continue
        w = count_ones(mask)
        orb = Set{UInt8}()
        for g in G
            img = UInt8(0)
            for b in 0:6
                if ((mask >> b) & 0x01) == 0x01
                    img |= UInt8(1) << g[b + 1]
                end
            end
            push!(orb, img)
        end
        union!(seen, orb)
        push!(orb_by_w[w], length(orb))
    end
    ok(length(G) == 168 &&
       Dict(w => sort(orb_by_w[w]) for w in 1:6) ==
       Dict(1 => [7], 2 => [21], 3 => [7, 28], 4 => [7, 28], 5 => [21], 6 => [7]),
       "App A.2: |PSL_2(7)| = 168 and subset orbit sizes [7, 21, 7+28, 7+28, 21, 7]")

    Bs = [k * (k + 1) for k in 1:8]
    As = [b + 1 for b in Bs]
    ok(sum(Bs) == 240 && sum(As) == 248,
       "App B.2: sum_{n=1..8} B_n = 240 and sum_{n=1..8} A_n = 248")

    M13 = zeros(Float64, 13, 13)
    for off in (0, 6), i in 0:5, j in (i + 1):5
        if j != xor(i, 1)
            M13[off + i + 1, off + j + 1] = 1.0
            M13[off + j + 1, off + i + 1] = 1.0
        end
    end
    for i in 1:12
        M13[13, i] = 1.0
        M13[i, 13] = 1.0
    end
    spec13 = [round(Int, v) for v in sort(eigvals(Symmetric(M13)))]
    ok(Int(sum(M13)) == 72 && spec13 == vcat(fill(-2, 5), fill(0, 6), [4, 6]),
       "App B.3: 13-vertex joined-octahedra graph has 36 edges and spectrum {6, 4, 0^6, (-2)^5}")

    exit(failures == 0 ? 0 : 1)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_checks()
end
