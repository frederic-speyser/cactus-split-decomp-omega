# Numerical Extension of the Enumeration of Non-Plane Cactus Graphs to the Case Ω = {5, 6} - Exploratory computations

**Status: Feature-complete** (all seven planned scripts written and
cross-checked, v1.1) - this refers only to the scripts, not to any
scientific claim; see § About this exploration and § Progress for what
remains unproved.

## About this exploration

A *cactus graph* is a connected graph in which every edge lies on at most
one cycle. A companion paper (Speyser, 2026) enumerates *strict m-gonal
cacti* - cacti in which every block is a cycle of a single fixed length
*m*, in the free (non-plane) setting - for *m* ≥ 5, deriving closed-form
and asymptotic results via split-decomposition. Its concluding remarks
note that the method extends "without difficulty" to a finite set Ω of
admissible cycle lengths - a cactus mixing, say, pentagons and hexagons at
different cut vertices - but this extension was never carried out, for any
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

## Repository structure

```
python/   the six Python scripts (solver, cross-checks)
pari/     the one PARI/GP script (verify_pari_omega.gp)
```

## Core scripts

- **`mgonal_cactus_series_omega.py`** - computes the rooted and unrooted
  enumeration series for strict cactus graphs admitting a finite set Ω of
  cycle lengths, generalizing `mgonal_cactus_series.py` from the original
  repository: the kernel *K_C* becomes a sum of one USEQ term per size in
  Ω (§5.1 of the original paper, generalized), using the same exact
  rational formal power series arithmetic (Python `Fraction`).
  Cross-checked against the pure *m* = 5, *m* = 6 series of [1] at
  degree 21 (see CHANGELOG).
- **`growth_rate_omega.py`** - estimates the exponential growth rate 1/ρ_Ω
  from the coefficients computed above, via the same *n*−3/2-corrected
  ratio test used in `mgonal_cactus_growth_rate.py`, and compares it
  against the already-published values 1/ρ_5, 1/ρ_6. Independent of
  `critical_point_omega.py` (no shared code); result: 1/ρ_Ω ≈ 1.865 by
  raw ratio test, versus 1.882 from the critical-point method — a small
  gap consistent with the known slow convergence of the raw ratio test,
  not a contradiction (see CHANGELOG).
- **`critical_point_omega.py`** - solves the critical-point system
  Φ(ρ,τ)=τ, Φ_y(ρ,τ)=1 numerically for Ω = {5, 6}. As predicted by hand
  before running it (the *m* = 6 branch of *K_C* reintroduces a term
  linear in *y* that blocks the substitution of Theorem 2 of the original
  paper), no closed form for τ_Ω was found among several natural
  candidates. More notably: ρ_Ω = 0.531336 is smaller than *both* ρ_5
  (0.604765) and ρ_6 (0.633235) — the growth rate exceeds both pure cases
  rather than interpolating between them (see CHANGELOG for verification
  details).

## Supplementary verification scripts

- **`exhaustive_iso_omega.py`** - builds strict cacti for Ω = {5, 6} with
  1 and 2 blocks directly as graphs (no functional equation involved),
  including configurations explicitly mixing a pentagon and a hexagon,
  and deduplicates by graph isomorphism (via `networkx`) — in the same
  spirit as `exhaustive_iso.py` in the original repository. Matches the
  solver's *k* = 1, *k* = 2 unrooted coefficients exactly (see CHANGELOG).
- **`split_tree_omega.py`** - a brute-force split-decomposition search
  (Definition 1 of the original paper), extending `split_tree_v2.py` to
  test whether Theorem 1's characterization still holds unchanged when
  two different cycle sizes can meet at the same cut vertex — a
  structural question the original paper never addresses, since it only
  ever treats a single fixed *m*. Result: it does, on the cases tested
  (a pentagon glued to a hexagon, plus the same three negative controls
  as the original script). Two bugs were found and fixed in this script
  during development (see CHANGELOG) — the result is evidence, not proof,
  of the wider claim.
- **`verify_dissymmetry_omega.py`** - verifies the unrooted series G(x)
  via the dissymmetry theorem for mixed Ω, where the re-rooted term T_Cm
  becomes a sum of dihedral cycle indices Z_D5 + Z_D6 rather than a
  single one — extending `verify_dissymmetry_all_m.py` to the mixed case.
  Uses `sympy.Rational` with hand-written convolutions rather than
  `fractions.Fraction`. Matches the main solver exactly on the first 8
  non-zero terms (see CHANGELOG for a documented dead end — a fully
  symbolic SymPy approach that was too slow and was abandoned).
- **`verify_pari_omega.gp`** - a second, independent solver in PARI/GP,
  using native truncated power series arithmetic, as a cross-check of
  `mgonal_cactus_series_omega.py` by a different code path (and language)
  entirely — extending `verify_pari.gp` to the mixed-Ω case. Matches both
  Python implementations exactly, and extends the check to x^40 (see
  CHANGELOG for three non-mathematical PARI/GP pitfalls encountered and
  documented along the way).

## Repository structure

```
python/    all Python scripts (six of the seven)
pari/      the PARI/GP script (verify_pari_omega.gp)
```

## Usage

```bash
python3 python/mgonal_cactus_series_omega.py --omega 5,6
python3 python/growth_rate_omega.py --omega 5,6
python3 python/critical_point_omega.py --omega 5,6
python3 python/exhaustive_iso_omega.py       # requires: pip install networkx
python3 python/split_tree_omega.py
python3 python/verify_dissymmetry_omega.py --omega 5,6   # requires: pip install sympy
gp -q pari/verify_pari_omega.gp < /dev/null               # requires: PARI/GP
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
- [x] Independent cross-check of ρ_{5,6} by the coefficient-ratio method
      — done, agrees with `critical_point_omega.py` within ~1%
      (`growth_rate_omega.py`)
- [x] Search for a closed form for τ_Ω, or evidence of a structural
      obstruction comparable to Proposition 1 of [1] — no closed form
      found among natural candidates tested (`critical_point_omega.py`)
- [x] Check whether Theorem 1's split-decomposition characterization
      holds unchanged for mixed cycle sizes — confirmed on a pentagon +
      hexagon pair and the standard negative controls, after fixing two
      bugs in the test script itself (`split_tree_omega.py`)
- [x] Independent verification of the unrooted series via the dissymmetry
      theorem — matches exactly on the first 8 non-zero terms
      (`verify_dissymmetry_omega.py`)
- [x] Second, independent solver in PARI/GP — matches both Python
      implementations exactly, extends the check to x^40
      (`verify_pari_omega.gp`)

This list will be updated as the work progresses; nothing above should be
taken as established until its box is checked.

## Acknowledgments

None yet. This section will be updated as this exploration develops and is
discussed with others.

## Citation

A citable archive of this repository is available via Zenodo. The latest
tagged version is v1.0
([10.5281/zenodo.21838871](https://doi.org/10.5281/zenodo.21838871)); a
v1.1 archive ("Complete toolkit" — all seven scripts in place, organized
into `python/` and `pari/`) is expected shortly, see CHANGELOG. No
written account of this work has been deposited yet — see § About this
exploration for the working paper planned once the computations are
further along.

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
