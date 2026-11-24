# DICEBOX-V

## I. Overview

`DICEBOX-V` (where "V" stands for *Versatile*) is an independent, multi-threaded OpenMP fork of `DICEBOX-1.0.1` for Monte Carlo simulations of nuclear γ-ray cascades within the extreme statistical model. It is developed and maintained independently of the original `DICEBOX-1.0.1` authors, who are unaffiliated with this release and neither maintain nor assume any responsibility for this modified code.

Preserving the cascade physics and full input compatibility of `DICEBOX-1.0.1`, `DICEBOX-V` runs one to two orders of magnitude faster (`12.8×–15.4×` single-thread acceleration and `210×–236×` cumulative pure-compute speedups on a 24-core server) and adds coupled temperature-dependent NLD and 2D PSF tables. For further details, refer to:

- The accompanying manuscript submitted to *Computer Physics Communications* (algorithms, OpenMP design, consistency fixes, validation, and benchmarks; see Section V).
- [DICEBOX_V_Supplement.pdf](doc/DICEBOX_V_Supplement.pdf) (summary of updates, installation, coupled 2D table input design, and example summary table).
- [dicebox_users_manual_1.0.1.pdf](doc/dicebox_users_manual_1.0.1.pdf) (standard `DICEBOX-1.0.1` input cards, analytical models, and output formats).

## II. Directory Structure

```text
DICEBOX-V/
├── LICENSE                            # GNU General Public License v3
├── Makefile                           # Build targets (all, v, origin, clean, help)
├── README.md                          # Quick-start guide
├── check_install.sh                   # Build + Example1 installation check in /tmp
├── doc/                               # User manuals (PDF)
│   ├── DICEBOX_V_Supplement.pdf       # Summary of updates, installation, 2D NLD/PSF input design
│   └── dicebox_users_manual_1.0.1.pdf # Standard DICEBOX-1.0.1 cards, analytical models, outputs
├── src/
│   └── dicebox_v.for                  # DICEBOX-V Fortran + OpenMP 5.0 source
├── v1.0.1/
│   └── dicebox_1.0.1.for              # Unmodified DICEBOX-1.0.1 serial reference source
└── examples/                          # Example1–Example6, Example7a, Example7b (+ reference_output/)
```

## III. Installation and Quick Start

**Requirements:** An OpenMP 5.0-capable `gfortran` compiler and GNU `make`.

**Build targets:** `make` / `make all` (`bin/dicebox_v` and `bin/dicebox_1.0.1`), `make v` (`bin/dicebox_v` only), `make origin` (`bin/dicebox_1.0.1` only), `make clean`, and `make help`.

```bash
make v
bash check_install.sh                 # Builds bin/dicebox_v if needed & verifies Example1 in /tmp

cd examples/Example1
export OMP_STACKSIZE=32M              # Recommended when the default OpenMP thread stack is small
OMP_NUM_THREADS=4 ../../bin/dicebox_v
```

Set `OMP_NUM_THREADS` to the available physical core count (`> 128` is clamped to `128`). Static common-block memory usage is `720.4 MB` plus `~0.35 MB` per thread on the stack.

## IV. Bundled Examples

Eight self-contained example decks (`Example1`–`Example6`, `Example7a`, and `Example7b`; all with `NEVENTS = 1000`) are bundled under `examples/`, covering analytical and tabulated NLD/PSF models, two-channel and resonance capture, and the coupled temperature-dependent 2D pair (`Example7a`/`Example7b`). Every example directory includes a complete set of precomputed sample outputs in `examples/<Deck>/reference_output/` (`DICE.PRO`, `PSF_GS.DAT`, `PSF_INI.DAT`, `EVENTS.*`, and `LEVELS.*` when `ISWLS = 1`). See [DICEBOX_V_Supplement.pdf](doc/DICEBOX_V_Supplement.pdf) for the full parameter summary table.

## V. License and Citations

`DICEBOX-V` is distributed under the GNU General Public License v3 (see [LICENSE](LICENSE)). The unmodified `DICEBOX-1.0.1` reference code in [dicebox_1.0.1.for](v1.0.1/dicebox_1.0.1.for) is included in accordance with GPL redistribution terms. Users publishing results obtained with `DICEBOX-V` are requested to cite:

1. N. N. Anh, T. Q. Viet, L. T. Phuc, P. B. Q. Hieu, N. X. Hai, and N. Q. Hung, "DICEBOX-V: A Versatile and High-Performance Multi-Threaded Extension of the DICEBOX Statistical Model Code," submitted to *Comput. Phys. Commun.* (2026).
2. F. Bečvář, "Simulation of γ cascades in complex nuclei with emphasis on assessment of uncertainties of cascade-related quantities," *Nucl. Instrum. Methods Phys. Res. A* **417** (1998) 434–449.
3. M. Krtička, *DICEBOX Version 1.0 User Manual*, Charles University, Prague (`doc/dicebox_users_manual_1.0.1.pdf`).
