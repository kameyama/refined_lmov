# Refined Chern–Simons invariants and refined BPS invariants

English | [日本語](README.ja.md)

A research repository on refined Chern–Simons invariants and refined BPS (LMOV) invariants of knots. It holds the code, the tests, the computed data and the manuscripts for the themes in the table below. Themes are referred to by their directory names.

| Theme | Content | Entry point |
|---|---|---|
| framed_knot | Extends Kameyama–Nawata, "Refined large N duality for knots" (arXiv:1703.05408), to knots other than torus knots and to framed knots, and tests the positivity of the refined LMOV invariants | [progress.md](docs/framed_knot/progress.md), [REPORT.md](docs/framed_knot/reports/REPORT.md), [conventions.md](docs/framed_knot/conventions.md) |
| dgr_lmov | DGR-type and quadruply graded structures and refined LMOV invariants | [README.md](docs/dgr_lmov/README.md) |
| refined_kauffman | Large N duality for the refined Kauffman polynomial | [README.md](docs/refined_kauffman/README.md) |

The working rules are in [CLAUDE.md](CLAUDE.md) and .claude/rules/.

## Papers

| Theme | Title | Source |
|---|---|---|
| framed_knot | Refined large N duality for framed knots and positivity of refined LMOV invariants beyond torus knots | [docs/framed_knot/paper/en/](docs/framed_knot/paper/en/) ([PDF](docs/framed_knot/paper/en/main.pdf)) |
| refined_kauffman | Refined crosscap invariants and the large N duality for colored Kauffman invariants | [docs/refined_kauffman/paper/en/](docs/refined_kauffman/paper/en/) ([PDF](docs/refined_kauffman/paper/en/main.pdf)) |

## Language

This research was carried out in dialogue with an AI assistant (Claude, Anthropic), and the dialogue was conducted in Japanese. For that reason the working documents under docs/ (research plans, progress logs, conventions, reports), the working rules (CLAUDE.md, .claude/), the commit messages and many of the code comments are written in Japanese. The papers listed above and this README are in English.

This English README is the primary version. [README.ja.md](README.ja.md) is its Japanese counterpart, and the two are kept in step.

## Layout

The files of each theme are placed in a subdirectory named after the theme inside each top-level directory.

```
docs/framed_knot/     framed_knot: research plan, novelty report, conventions, progress, results, literature notes (literature/), paper (paper/)
docs/dgr_lmov/        dgr_lmov documents (start from README.md)
docs/refined_kauffman/  refined_kauffman documents and paper (start from README.md)
references/           paper PDFs and TeX sources (not tracked; sources are listed in references/README.md)
src/rcsknots/         computation code (a Python package shared by all themes)
src/rcsknots/dgr_lmov/  modules used only by dgr_lmov (dgr, lmovstruct, large_n_inverse)
src/rcsknots/refined_kauffman/  modules used only by refined_kauffman (koornwinder, daha_dn, daha_cn, flint_rf*, kauffman_*)
tests/framed_knot/    verification tests for framed_knot (pytest)
tests/dgr_lmov/       verification tests for dgr_lmov
tests/refined_kauffman/  verification tests for refined_kauffman
scripts/framed_knot/  one-off computation and verification scripts for framed_knot
scripts/dgr_lmov/     scripts for dgr_lmov (they also import scripts in framed_knot/)
scripts/refined_kauffman/  scripts for refined_kauffman (rk_*.py)
data/framed_knot/     computed results for framed_knot (JSON and human-readable tables)
data/dgr_lmov/        computed results for dgr_lmov (including lmovstruct/)
data/refined_kauffman/  computed results for refined_kauffman (pkl files in daha/ and others)
data/cache/           cache of intermediate results (not tracked; shared by all themes)
logs/framed_knot/     verification logs and records of run time and memory for framed_knot
logs/dgr_lmov/        logs for dgr_lmov
logs/refined_kauffman/  logs for refined_kauffman
notebooks/            scratch work
```

Scripts are run with `PYTHONPATH=src:scripts/framed_knot:scripts/dgr_lmov` (for refined_kauffman, `PYTHONPATH=src:scripts/refined_kauffman`). Imports between scripts are resolved through this path.

## Setup and tests

```bash
just setup
```

```bash
just fetch-refs
```

```bash
.venv/bin/pytest
```

## Environment

As checked on 2026-10-05.

| Tool | Status | Version | How to run |
|---|---|---|---|
| Python | available | 3.14.8 (Homebrew) | `.venv/bin/python` |
| sympy / numpy / mpmath | available | 1.14.0 / 2.5.3 / 1.3.0 | .venv |
| python-flint | available | 0.9.0 | .venv (`nmod`, `nmod_mat`, `nmod_poly`, `fmpz_mpoly`, etc.) |
| pytest | available | 9.1.1 | `.venv/bin/pytest` |
| khoca | available (sl(2) only; for N ≥ 3 it returns wrong values for the unknot, so it is not used) | 1.5 | .venv. Built against pari from the Sage conda environment: `CFLAGS=-I$S/include LDFLAGS="-L$S/lib -Wl,-rpath,$S/lib" pip install khoca` (`S` is the path of the conda environment `sage`). Requires the `py` package |
| SnapPy | available | 3.3.2 (knot_floer_homology 1.2.2) | .venv |
| scipy | available | 1.18.1 | .venv (`milp`, integer programming) |
| SageMath | available | 10.9 (conda-forge, environment `sage`, Python 3.14.7) | `DOT_SAGE=$TMPDIR/dotsage /opt/homebrew/Caskroom/miniforge/base/envs/sage/bin/sage file.sage` |
| Wolfram Engine | available (only outside the Claude Code sandbox) | 15.0.0 for Mac OS X ARM (64-bit) | `wolframscript -code ...` / `wolframscript -file x.wls` |

Minimal checks of SageMath (all of them worked; about 12 seconds in total):

- `SymmetricFunctions(QQ['q','t'].fraction_field()).macdonald().P()`: P[2,1] is displayed in the power-sum basis.
- `macdonald().Ht()`: H̃[2,1] = s[3] + (q+t) s[2,1] + qt s[1,1,1].
- Plethysm: `p(P[2]).plethysm(p[2])` is expanded in the P basis.
- `Link(BraidGroup(2)([1,1,1])).homfly_polynomial()`: L^-2 M^2 − 2L^-2 − L^-4 (the trefoil).

Wolfram Engine (activated by the author on 2026-10-05):

- Inside the sandbox it cannot read its configuration and licence and fails with "A WolframKernel location could not be determined". Run it outside the sandbox.
- The notebooks of arXiv:1703.05408 can be run by extracting the input cells with `NotebookImport[file, "Input" -> "HeldExpression"]`.
  Example: after evaluating the definition cells of the non-torus notebook, `FunctionExpand[torusknotsym[1,1,a,q,t]/torusknotsym[1,0,a,q,t]]` gives
  a(1 + qt − aq)/√(qt) (the reduced value for the trefoil). It also agrees with the formula ported to Python (`src/rcsknots/notebook_ref.py`).
