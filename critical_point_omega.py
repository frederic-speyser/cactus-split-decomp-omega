"""
critical_point_omega.py

Numerically solves the critical-point system Phi(rho,tau)=tau,
Phi_y(rho,tau)=1 for Omega = {5, 6}, and checks whether tau_Omega admits
any closed form recognizable by comparison with the pure tau_5, tau_6
values of [1], or by simple pattern search (PSLQ-free: just a handful of
natural algebraic guesses, since no PSLQ implementation is used here).

Structural expectation, worked out by hand before writing this script
(see the accompanying journal entry): for Omega={5,6},

  K_C(x,y) = (1/2)y^4 + (1/2)y^5 + (1/2)s(x^2)^2 * (1+y)

so the criticality condition Phi_y(rho,tau)=1 reduces, after substituting
x*exp(K_C+h)=tau from Phi(rho,tau)=tau, to

  tau * [ 2*tau^3 + (5/2)*tau^4 + (1/2)*s(rho^2)^2 ] = 1

This involves s(rho^2), not expressible from (rho,tau) alone -- the same
structural obstruction as Proposition 1 of [1] for m even, inherited here
because the m=6 branch contributes a term linear in y. No closed form is
therefore expected by the method of [1]'s Theorem 2. This script checks
that expectation computationally and searches for a closed form anyway,
as a sanity check, rather than assuming the prediction is correct.

This script does not prove anything. It only reports what the numbers
show.

Reference: Fr. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
Non-Plane m-Gonal Cactus Graphs via Split-Decomposition" [1], Theorem 2
and Proposition 1 (the closed-form / obstruction results for pure m).

Author: Frederic G. Speyser
Run: python3 critical_point_omega.py --omega 5,6
"""
import argparse
import numpy as np


def solve_s_fast(omega, N):
    """Rooted series s(x) up to x^N, double precision, order-by-order,
    for a mixed Omega. Mirrors the structure of solve_s_fast in
    asymptotic_convergence.py from [1], generalized to a set of sizes."""
    s = np.zeros(N + 1)
    s[1] = 1.0
    for n in range(2, N + 1):
        s_t = s
        sp = np.zeros(N + 1)
        sp[0:N + 1:2] = s_t[:len(sp[0:N + 1:2])]

        KC = np.zeros(N + 1)
        for m in omega:
            if m % 2 == 1:
                spk = np.zeros(N + 1); spk[0] = 1.0
                for _ in range((m - 1) // 2):
                    spk = np.convolve(spk, sp)[:N + 1]
                sk = np.zeros(N + 1); sk[0] = 1.0
                for _ in range(m - 1):
                    sk = np.convolve(sk, s_t)[:N + 1]
                KC += 0.5 * (sk + spk)
            else:
                spk = np.zeros(N + 1); spk[0] = 1.0
                for _ in range((m - 2) // 2):
                    spk = np.convolve(spk, sp)[:N + 1]
                sk = np.zeros(N + 1); sk[0] = 1.0
                for _ in range(m - 1):
                    sk = np.convolve(sk, s_t)[:N + 1]
                skp = np.convolve(s_t, spk)[:N + 1]
                KC += 0.5 * (sk + skp)

        G = np.zeros(N + 1)
        min_deg = min(m - 1 for m in omega)
        i = 1
        while i * min_deg <= n:
            kc_i = np.zeros(N + 1)
            kc_i[0:N + 1:i] = KC[:len(kc_i[0:N + 1:i])]
            G += kc_i / i
            i += 1
        E = np.zeros(N + 1); E[0] = 1.0
        for k in range(1, n + 1):
            E[k] = sum(j * G[j] * E[k - j] for j in range(1, k + 1)) / k
        s[n] = E[n - 1]
    return s


def eval_series(coeffs, x):
    total, xp = 0.0, 1.0
    for c in coeffs:
        total += c * xp
        xp *= x
        if xp < 1e-300:
            break
    return total


def K_C_and_dKCdy(omega, x, y, s_at_x2):
    """K_C(x,y) and d K_C/dy, given s(x^2) already evaluated."""
    KC = 0.0
    dKC = 0.0
    for m in omega:
        if m % 2 == 1:
            KC += 0.5 * (y ** (m - 1) + s_at_x2 ** ((m - 1) // 2))
            dKC += 0.5 * (m - 1) * y ** (m - 2)
        else:
            KC += 0.5 * (y ** (m - 1) + y * s_at_x2 ** ((m - 2) // 2))
            dKC += 0.5 * ((m - 1) * y ** (m - 2) + s_at_x2 ** ((m - 2) // 2))
    return KC, dKC


def find_critical_point(omega, s_coeffs, rho_lo=0.3, rho_hi=0.9, iters=200):
    """Bisection on rho: for a candidate rho, tau is defined implicitly by
    Phi(rho,tau)=tau (solved by fixed-point iteration at fixed rho using
    the already-known series), then check whether Phi_y(rho,tau) is above
    or below 1, and bisect on rho until Phi_y(rho,tau)=1."""

    def tau_at(rho):
        # s(rho) itself, evaluated from the coefficients, IS tau once rho
        # is truly the radius of convergence; for a trial rho inside the
        # true radius, s(rho) is just a finite, well-defined number.
        return eval_series(s_coeffs, rho)

    def phi_y_at(rho):
        tau = tau_at(rho)
        s_rho2 = eval_series(s_coeffs, rho ** 2)
        _, dKC = K_C_and_dKCdy(omega, rho, tau, s_rho2)
        # Phi_y(rho,tau) = rho * exp(KC+h) * dKC/dy = tau * dKC/dy,
        # using rho*exp(KC+h) = tau at any point where s(rho)=tau exactly
        # (true for any rho within the radius of convergence, by definition
        # of s as the sum of its own series).
        return tau * dKC

    lo, hi = rho_lo, rho_hi
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        val = phi_y_at(mid)
        if val < 1:
            lo = mid
        else:
            hi = mid
    rho = 0.5 * (lo + hi)
    tau = tau_at(rho)
    return rho, tau


CLOSED_FORM_CANDIDATES = {
    "2^(-1/9)": 2 ** (-1 / 9),
    "2^(-1/10)": 2 ** (-1 / 10),
    "(2/9)^(1/9)": (2 / 9) ** (1 / 9),
    "(2/10)^(1/10)": (2 / 10) ** (1 / 10),
    "sqrt(rho_5 * rho_6) [geometric mean]": None,   # filled in at runtime
    "(rho_5 + rho_6)/2 [arithmetic mean]": None,     # filled in at runtime
}

# Reference values from [1], Table 3
RHO_5, TAU_5 = 0.604765, 0.840896
RHO_6, TAU_6 = 0.633235, 0.821008


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--omega", type=str, default="5,6")
    parser.add_argument("--N", type=int, default=1200,
                         help="Truncation order for the series used to "
                              "evaluate s(x) numerically")
    args = parser.parse_args()
    omega = tuple(sorted(int(x) for x in args.omega.split(",")))

    print(f"Omega = {{{', '.join(map(str, omega))}}}")
    print(f"Computing s(x) numerically up to N={args.N}...")
    s_coeffs = solve_s_fast(omega, args.N)

    rho, tau = find_critical_point(omega, s_coeffs)
    print(f"\nrho_Omega  = {rho:.6f}")
    print(f"tau_Omega  = {tau:.6f}")

    print(f"\nFor comparison (pure cases, from [1], Table 3):")
    print(f"  rho_5 = {RHO_5:.6f}   tau_5 = {TAU_5:.6f}")
    print(f"  rho_6 = {RHO_6:.6f}   tau_6 = {TAU_6:.6f}")

    print(f"\nDoes rho_Omega interpolate between rho_5 and rho_6?")
    between = RHO_5 < rho < RHO_6
    print(f"  rho_5 < rho_Omega < rho_6 : {between}  "
          f"({RHO_5:.6f} < {rho:.6f} < {RHO_6:.6f})")

    print(f"\nClosed-form pattern search (not a proof -- just checking a")
    print(f"handful of natural guesses against tau_Omega = {tau:.6f}):")
    CLOSED_FORM_CANDIDATES["sqrt(rho_5 * rho_6) [geometric mean]"] = (RHO_5 * RHO_6) ** 0.5
    CLOSED_FORM_CANDIDATES["(rho_5 + rho_6)/2 [arithmetic mean]"] = (RHO_5 + RHO_6) / 2
    found_match = False
    for label, val in CLOSED_FORM_CANDIDATES.items():
        diff = abs(val - tau)
        match = diff < 1e-4
        found_match = found_match or match
        print(f"  {label:38s} = {val:.6f}   |diff| = {diff:.6f}   "
              f"{'MATCH' if match else ''}")

    print(f"\nStructural prediction (worked out by hand before running this")
    print(f"script -- see journal entry): the m=6 branch of K_C reintroduces")
    print(f"a term linear in y that does not cancel, so the substitution")
    print(f"used in Theorem 2 of [1] for m odd should NOT close here, and no")
    print(f"closed form is expected by that method.")
    print(f"Closed-form match found among candidates tested: {found_match}")
    print(f"(A negative result here is not a proof of non-existence -- only")
    print(f"a proof would settle this. See Progress in the README.)")


if __name__ == "__main__":
    main()
