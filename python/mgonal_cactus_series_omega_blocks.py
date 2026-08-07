"""
mgonal_cactus_series_omega_blocks.py

Computes strict cactus enumeration for a mixed set Omega of cycle lengths,
indexed by NUMBER OF BLOCKS k -- the convention actually used by the
published OEIS sequences for this family (A398033-A398035, A397210,
A397546, etc.), which are all indexed by block count, not vertex count.

This distinction is invisible for a SINGLE fixed m (vertex count n and
block count k are related by the bijection n = 1+(m-1)k, so "the list of
non-zero x^n coefficients in order" already equals the by-k sequence).
It stops being invisible for a mixed Omega: a fixed k admits several
possible n (depending on the split between block sizes), and a fixed n
can arise from several different k. Reindexing by k therefore requires a
genuinely bivariate computation -- x marking vertices, u marking blocks
-- not just a relabeling of the vertex-indexed series already computed by
mgonal_cactus_series_omega.py.

Validated against the real, already-published A398035 (m=6, unrooted, by
blocks): this script reproduces 1, 1, 4, 13, 67, 372, ... exactly when
run with Omega=(6,) -- see the __main__ block below.

Reference: Fr. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
Non-Plane m-Gonal Cactus Graphs via Split-Decomposition" [1], Sections
5.1-5.3, and the Ω-mixed generalization in mgonal_cactus_series_omega.py
from this repository.

Author: Frederic G. Speyser
Run: python3 mgonal_cactus_series_omega_blocks.py --omega 5,6 --terms 10
"""
import argparse
from fractions import Fraction as F

NX = 60  # truncation in x (vertices); raised automatically if needed
NU = 12  # truncation in u (blocks) -- how many block-count terms to compute


def zero():
    return [[F(0)] * (NU + 1) for _ in range(NX + 1)]


def add(a, b):
    return [[a[n][k] + b[n][k] for k in range(NU + 1)] for n in range(NX + 1)]


def sub(a, b):
    return [[a[n][k] - b[n][k] for k in range(NU + 1)] for n in range(NX + 1)]


def scale(a, c):
    return [[a[n][k] * c for k in range(NU + 1)] for n in range(NX + 1)]


def mul(a, b):
    c = zero()
    for n1 in range(NX + 1):
        for k1 in range(NU + 1):
            v1 = a[n1][k1]
            if v1 == 0:
                continue
            for n2 in range(NX + 1 - n1):
                for k2 in range(NU + 1 - k1):
                    v2 = b[n2][k2]
                    if v2 == 0:
                        continue
                    c[n1 + n2][k1 + k2] += v1 * v2
    return c


def stretch(a, i):
    """Substitute x -> x^i AND u -> u^i simultaneously (both the vertex
    count and the block count of a repeated/doubled sub-structure scale
    together by the same factor i)."""
    c = zero()
    for n in range(NX + 1):
        if n * i > NX:
            break
        for k in range(NU + 1):
            if k * i > NU:
                break
            c[n * i][k * i] = a[n][k]
    return c


def power_int(a, p):
    r = zero()
    r[0][0] = F(1)
    base = a
    while p > 0:
        if p & 1:
            r = mul(r, base)
        base = mul(base, base)
        p >>= 1
    return r


def shift_by_x(a):
    c = zero()
    for n in range(NX):
        for k in range(NU + 1):
            c[n + 1][k] = a[n][k]
    return c


def exp_series_in_i(terms_by_i):
    """exp(Sum_i terms_by_i[i]), differentiating w.r.t. x only (u is a
    spectator variable carried along coefficient-wise)."""
    u_total = zero()
    for t in terms_by_i.values():
        u_total = add(u_total, t)
    v = zero()
    v[0][0] = F(1)
    for n in range(1, NX + 1):
        for k in range(NU + 1):
            acc = F(0)
            for a in range(0, n + 1):
                for b in range(0, k + 1):
                    if a == 0 and b == 0:
                        continue
                    if u_total[a][b] == 0:
                        continue
                    acc += a * u_total[a][b] * v[n - a][k - b]
            v[n][k] = acc / n
    return v


def K_C_single(s, m):
    """Bivariate kernel for one block of size m. The leading u^1 shift
    marks exactly this one block; the internal branches (s(x,u) or
    s(x^2,u^2)) propagate block-counts from the recursive sub-cacti
    hanging off this block's other markers."""
    s2 = stretch(s, 2)
    if m % 2 == 1:
        term1 = power_int(s, m - 1)
        term2 = power_int(s2, (m - 1) // 2)
        kernel = scale(add(term1, term2), F(1, 2))
    else:
        term1 = power_int(s, m - 1)
        term2 = mul(s, power_int(s2, (m - 2) // 2))
        kernel = scale(add(term1, term2), F(1, 2))
    out = zero()
    for n in range(NX + 1):
        for k in range(NU):
            out[n][k + 1] = kernel[n][k]
    return out


def K_C(s, omega):
    total = zero()
    for m in omega:
        total = add(total, K_C_single(s, m))
    return total


def outer_symmetrization_terms(kc_base, omega):
    """Stretch the WHOLE K_C result (including its own u^1 block-marker)
    by i, rather than recomputing K_C from a pre-stretched s: the block's
    own marker must scale together with everything else when a block is
    treated as one repeated unit under the outer multiset construction."""
    terms_by_i = {}
    i = 1
    min_deg = min(m - 1 for m in omega)
    while i * min_deg <= NX:
        terms_by_i[i] = scale(stretch(kc_base, i), F(1, i))
        i += 1
    return terms_by_i


def solve_s(omega, iters=None):
    if iters is None:
        iters = NX + 2
    s = zero()
    s[1][0] = F(1)
    for _ in range(iters):
        kc_base = K_C(s, omega)
        terms_by_i = outer_symmetrization_terms(kc_base, omega)
        E = exp_series_in_i(terms_by_i)
        Eminus1 = [row[:] for row in E]
        Eminus1[0][0] -= 1
        s_new = zero()
        s_new[1][0] = F(1)
        s_new = add(s_new, shift_by_x(Eminus1))
        s = s_new
    return s


def phi_totient(n):
    result = n
    p = 2
    nn = n
    while p * p <= nn:
        if nn % p == 0:
            while nn % p == 0:
                nn //= p
            result -= result // p
        p += 1
    if nn > 1:
        result -= result // nn
    return result


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def Z_Dm(s, m):
    p = {i: stretch(s, i) for i in range(1, m + 1)}
    total = zero()
    for d in divisors(m):
        term = power_int(p[d], m // d)
        total = add(total, scale(term, F(phi_totient(d), 2 * m)))
    if m % 2 == 1:
        extra = mul(p[1], power_int(p[2], (m - 1) // 2))
        total = add(total, scale(extra, F(1, 2)))
    else:
        extra1 = mul(power_int(p[1], 2), power_int(p[2], (m - 2) // 2))
        extra2 = power_int(p[2], m // 2)
        total = add(total, scale(add(extra1, extra2), F(1, 4)))
    out = zero()
    for n in range(NX + 1):
        for k in range(NU):
            out[n][k + 1] = total[n][k]
    return out


def solve_G(omega, iters=None):
    s = solve_s(omega, iters=iters)
    kc_base = K_C(s, omega)
    terms_by_i = outer_symmetrization_terms(kc_base, omega)
    E = exp_series_in_i(terms_by_i)
    KC = kc_base
    Eminus1 = [row[:] for row in E]
    Eminus1[0][0] -= 1
    S_X = shift_by_x(Eminus1)
    S_C = sub(Eminus1, KC)
    T_S = shift_by_x(S_C)
    T_SCm = mul(KC, S_X)
    T_Cm = zero()
    for m in omega:
        T_Cm = add(T_Cm, Z_Dm(s, m))
    G = add(sub(T_Cm, T_SCm), T_S)
    return s, G


def total_by_blocks(bivariate_series, max_k):
    """Sum over all n of the coefficient of u^k, for k=1..max_k -- finite
    for each k since only finitely many n are reachable with exactly k
    blocks."""
    return [sum(bivariate_series[n][k] for n in range(NX + 1))
            for k in range(1, max_k + 1)]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=str, default="5,6")
    parser.add_argument("--terms", type=int, default=10)
    args = parser.parse_args()
    omega = tuple(sorted(int(v) for v in args.omega.split(",")))

    NX = 1 + max(omega) * args.terms
    NU = args.terms + 1

    print(f"Omega = {{{', '.join(map(str, omega))}}}")
    print(f"Computing rooted and unrooted series by block count, "
          f"k=1..{args.terms} (NX={NX}, NU={NU})...")
    s, G = solve_G(omega, iters=NX + 2)
    rooted = total_by_blocks(s, args.terms)
    unrooted = total_by_blocks(G, args.terms)

    print(f"\nRooted, by k=1..{args.terms}:")
    print("  " + ", ".join(str(int(v)) for v in rooted))
    print(f"\nUnrooted, by k=1..{args.terms}:")
    print("  " + ", ".join(str(int(v)) for v in unrooted))

    if omega == (6,):
        expected = [1, 1, 4, 13, 67, 372, 2419, 16551]
        got = [int(v) for v in unrooted[:len(expected)]]
        print(f"\nValidation (Omega={{6}} only): matches published A398035 "
              f"{'YES' if got == expected else 'NO -- DO NOT TRUST OUTPUT'}")
