# Numerical Extension of the Enumeration of Non-Plane Cactus Graphs to the Case Ω = {5, 6} - exploratory computations

## About this exploration

A *cactus graph* is a connected graph in which every edge lies on at most
one cycle. A companion paper (Speyser, 2026) enumerates *strict m-gonal
cacti* — cacti in which every block is a cycle of a single fixed length
*m*, in the free (non-plane) setting — for *m* ≥ 5, deriving closed-form
and asymptotic results via split-decomposition. Its concluding remarks
note that the method extends "without difficulty" to a finite set Ω of
admissible cycle lengths — a cactus mixing, say, pentagons and hexagons at
different cut vertices — but this extension was never carried out, for any
Ω, either numerically or analytically.

This repository begins that execution, starting from the case Ω = {5, 6}.
Unlike the paper it extends, nothing here has been proved, submitted, or
peer-reviewed: this is a numerical and computational exploration, not a
completed piece of research. It is independent of a second paper [2],
which proves the original paper's Conjecture 1 (strict monotonicity of
the growth rate) but explicitly excludes the mixed-Ω case from its own
scope.

Exploratory computations extending:

> Fr. G. Speyser, *Enumeration and Asymptotic Analysis of Strict Non-Plane m-Gonal Cactus Graphs via Split-Decomposition*, 2026.
> Submitted to *The Electronic Journal of Combinatorics*.
> Preprint: https://doi.org/10.5281/zenodo.21513753

A written account of this work, once complete, is intended to be deposited
on Zenodo as a working paper (not a journal submission) titled *"Numerical
Extension of the Enumeration of Non-Plane Cactus Graphs to the Case
Ω = {5, 6}: Exploratory Computations"*.

## Core scripts

- **`mgonal_cactus_series_omega.py`** — computes the rooted and unrooted
  enumeration series for strict cactus graphs admitting a finite set Ω of
  cycle lengths, generalizing `mgonal_cactus_series.py` from the original
  repository: the kernel *K_C* becomes a sum of one USEQ term per size in
  Ω (§5.1 of the original paper, generalized), using the same exact
  rational formal power series arithmetic (Python `Fraction`).
  Cross-checked against the pure *m* = 5, *m* = 6 series of [1] at
  degree 21 (see CHANGELOG).
- **`growth_rate_omega.py`** *(planned)* — will estimate the exponential
  growth rate 1/ρ_Ω from the coefficients computed above, via the same
  *n*−3/2-corrected ratio test used in `mgonal_cactus_growth_rate.py`,
  and compare it against the already-published values 1/ρ_5, 1/ρ_6.
- **`critical_point_omega.py`** *(planned)* — will solve the critical-point
  system Φ(ρ,τ)=τ, Φ_y(ρ,τ)=1 numerically for Ω = {5, 6}, and check whether
  τ_Ω admits any recognizable closed form (by comparison with the pure
  τ_5, τ_6 and pattern search), or whether the obstruction of Proposition 1
  in the original paper has an analogue here. Whatever the outcome, it will
  be documented as such — this script is not expected to succeed.

## Supplementary verification scripts

- **`exhaustive_iso_omega.py`** — builds strict cacti for Ω = {5, 6} with
  1 and 2 blocks directly as graphs (no functional equation involved),
  including configurations explicitly mixing a pentagon and a hexagon,
  and deduplicates by graph isomorphism (via `networkx`) — in the same
  spirit as `exhaustive_iso.py` in the original repository. Matches the
  solver's *k* = 1, *k* = 2 unrooted coefficients exactly (see CHANGELOG).
- **`split_tree_omega.py`** *(planned)* — a brute-force split-decomposition
  search (Definition 1 of the original paper), extending `split_tree_v2.py`
  to test whether Theorem 1's characterization still holds unchanged when
  two different cycle sizes can meet at the same cut vertex — a structural
  question the original paper never addresses, since it only ever treats
  a single fixed *m*.
- **`verify_dissymmetry_omega.py`** *(planned)* — will verify the unrooted
  series G(x) via the dissymmetry theorem for mixed Ω, where the re-rooted
  term T_Cm becomes a sum of dihedral cycle indices Z_D5 + Z_D6 rather than
  a single one — extending `verify_dissymmetry_all_m.py` to the mixed case.
- **`verify_pari_omega.gp`** *(planned)* — a second, independent solver in
  PARI/GP, using native truncated power series arithmetic, as a cross-check
  of `mgonal_cactus_series_omega.py` by a different code path entirely —
  extending `verify_pari.gp` to the mixed-Ω case.

## Usage

```bash
python3 mgonal_cactus_series_omega.py --omega 5,6
python3 growth_rate_omega.py --omega 5,6
python3 critical_point_omega.py --omega 5,6
python3 exhaustive_iso_omega.py       # requires: pip install networkx
python3 split_tree_omega.py
python3 verify_dissymmetry_omega.py   # requires: pip install sympy
gp -q verify_pari_omega.gp
```

No dependencies beyond the Python standard library and `numpy`; `networkx`
is required only for `exhaustive_iso_omega.py`, and `sympy` only for
`verify_dissymmetry_omega.py`. `verify_pari_omega.gp` requires PARI/GP.

Actual output of `mgonal_cactus_series_omega.py --omega 5,6 --terms 25`:

```
Omega = {5, 6}
Truncation order used: N = 100

rooted (offset 0)  : x + x^5 + x^6 + 3x^9 + 6x^10 + 4x^11 + 13x^13 +
41x^14 + 49x^15 + 22x^16 + 62x^17 + 278x^18 + 498x^19 + 415x^20 +
473x^21 + 1920x^22 + 4600x^23 + 5693x^24 + 5547x^25 + 14359x^26 +
40326x^27 + 66324x^28 + 74199x^29 + 126743x^30 + 349403x^31

unrooted (offset 1): x^5 + x^6 + x^9 + x^10 + x^11 + 3x^13 + 6x^14 +
6x^15 + 4x^16 + 8x^17 + 25x^18 + 42x^19 + 32x^20 + 44x^21 + 140x^22 +
302x^23 + 357x^24 + 353x^25 + 848x^26 + 2192x^27 + 3391x^28 + 3759x^29 +
6300x^30 + 16348x^31 + 31201x^32
```

Not yet submitted to the OEIS — see Data availability below.

## Relation to prior work

Bahrani and Lumbroso's general split-decomposition template [3] already
covers the mixed-Ω case in principle (their Ω-parameterized grammar), but
neither their paper nor the original paper [1] this repository extends
instantiate it numerically for any non-singleton Ω. As far as independent
literature searches have shown, no numerical enumeration for a mixed
Ω = {m₁, m₂, ...} appears elsewhere for this class.

## Data availability

The 25 rooted and 25 unrooted terms shown above (§ Usage) have been
computed but not yet submitted to the OEIS — no search for a possible
prior match has been done yet. The *k* = 1 and *k* = 2 terms have been
independently cross-checked by direct construction
(`exhaustive_iso_omega.py`); everything beyond *k* = 2 still rests on the
solver alone.

## Progress

- [x] Extend the solver to non-singleton Ω (`mgonal_cactus_series_omega.py`)
- [x] Independent construction of the small cases, mixed sizes
      (`exhaustive_iso_omega.py`)
- [ ] Numerical estimate of ρ_{5,6} and comparison with ρ_5, ρ_6
      (`growth_rate_omega.py`)
- [ ] Search for a closed form for τ_Ω, or evidence of a structural
      obstruction comparable to Proposition 1 of [1]
      (`critical_point_omega.py`)
- [ ] Check whether Theorem 1's split-decomposition characterization
      holds unchanged for mixed cycle sizes (`split_tree_omega.py`)
- [ ] Independent verification of the unrooted series via the dissymmetry
      theorem (`verify_dissymmetry_omega.py`)
- [ ] Second, independent solver in PARI/GP (`verify_pari_omega.gp`)

This list will be updated as the work progresses; nothing above should be
taken as established until its box is checked.

## Acknowledgments

None yet. This section will be updated as this exploration develops and is
discussed with others.

## Citation

No citable version of this work exists yet. Once the working paper
mentioned above is deposited on Zenodo, its DOI will be added here,
alongside guidance for citing this repository's code directly.

## References

[1] Speyser, F. G. *Enumeration and Asymptotic Analysis of Strict Non-Plane
    m-Gonal Cactus Graphs via Split-Decomposition.* Submitted to *The
    Electronic Journal of Combinatorics*, 2026. Preprint: DOI
    [10.5281/zenodo.21513753](https://doi.org/10.5281/zenodo.21513753).

[2] Speyser, F. G. *Strict Monotonicity of the Growth Rate for Non-Plane
    Strict m-Gonal Cactus Graphs.* In preparation for submission, 2026.

[3] Bahrani, M., Lumbroso, J. *Split-Decomposition Trees with Prime Nodes:
    Enumeration and Random Generation of Cactus Graphs.* Proceedings of
    ANALCO 2018, pp. 143–157. DOI
    [10.1137/1.9781611975062.13](https://doi.org/10.1137/1.9781611975062.13).

## Author

Frédéric G. Speyser — Independent Researcher, France - https://orcid.org/0000-0002-1767-5325

## License

MIT (see `LICENSE`), for consistency with the original repository.
