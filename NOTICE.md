# Licenses and third-party material

This file says which license applies to which part of the repository, and lists the material that is not covered by either license.

## What is covered

| Part | License |
|---|---|
| Source code: the Python package in `src/`, the scripts in `scripts/` (Python, Sage, Wolfram Language), the tests in `tests/`, and the scripts under `.claude/` | MIT, see [LICENSE](LICENSE) |
| Computed data in `data/`, logs in `logs/`, and the documents in `docs/`, `.claude/`, `README.md`, `README.ja.md` and `CLAUDE.md` | [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/) |

Copyright (c) 2026 Masaya Kameyama. For CC BY 4.0, attribute the material to Masaya Kameyama and link to this repository.

If you use the code or the data in published work, please cite it as described in [CITATION.cff](CITATION.cff) and in the README.

## What is not covered

### Manuscripts

The manuscripts of the papers, in `docs/framed_knot/paper/`, `docs/refined_kauffman/paper/` and `docs/refined_kauffman/tex/`, are not released under CC BY 4.0. All rights are reserved for now. Their license will be the one chosen for the arXiv versions.

### Material from other works

The files below reproduce or are derived from the work of other authors, or of the author together with co-authors. They are included so that the computations can be checked against the literature. Neither license above applies to them, and the rights remain with the original authors. Paths are relative to `data/framed_knot/reference/` unless they start with `src/`.

| File | Source | Terms stated by the source |
|---|---|---|
| `yui_kr/homology.json` | Computed results of yui-kr (https://github.com/taketo1024/yui-kr), the program of K. Nakagane and T. Sano, "Computations of HOMFLY homology", arXiv:2111.00388 | MIT, Copyright (c) 2024 Taketo Sano; see `yui_kr/LICENSE` and `yui_kr/SOURCE.md` |
| `yui_kr/braid.csv`, `yui_kr/homfly.csv` | Braid representatives and HOMFLY-PT polynomials from KnotInfo, obtained through the yui-kr repository | KnotInfo states no license. It asks to be cited as: C. Livingston and A. H. Moore, *KnotInfo: Table of Knot Invariants*, knotinfo.org, with the date of access |
| `gs_1112.0030/s2_thick.py` | Formulas transcribed from S. Gukov and M. Stošić, "Homological algebra of knots and BPS states", arXiv:1112.0030 | arXiv non-exclusive license |
| `stosic_2402.03066/quivers.json` | Quiver data transcribed from M. Stošić, "Generalized knots-quivers correspondence", arXiv:2402.03066 | arXiv non-exclusive license |
| `daha_examples.txt` | Explicit examples from I. Cherednik and I. Danilenko, "DAHA and iterated torus knots", arXiv:1408.4348 | arXiv non-exclusive license |
| `kn_tables.json`, `notebook_torus.txt`, `src/rcsknots/notebook_ref.py` | Tables and ancillary notebooks of M. Kameyama and S. Nawata, "Refined large N duality for knots", arXiv:1703.05408 | arXiv non-exclusive license |
| `kntz_package/values.txt`, `src/rcsknots/double_twist_refined.py` | Output of, and a Python port following, the ancillary package of M. Kameyama, S. Nawata, R. Tao and H. D. Zhang, "Cyclotomic expansions of HOMFLY-PT colored by rectangular Young diagrams", arXiv:1902.02275 | arXiv non-exclusive license |

Other files in the repository quote formulas and numerical values from the literature for comparison. The source is given in each file, and those formulas and values remain the work of their authors.

### Libraries

No third-party library is bundled. The code depends on sympy, numpy, mpmath and scipy (BSD), python-flint (MIT and LGPL-3.0-or-later) and pytest (MIT). Two scripts call khoca and SnapPy, and three scripts are run with SageMath; these are distributed under the GPL and are installed separately by the user.
