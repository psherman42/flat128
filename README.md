# flat128

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23246819.svg)](https://doi.org/10.5281/zenodo.23246819)

**Depth or Breadth? Flat128 Arithmetic and the RV128 Opportunity**
Paul Sherman (RISC-V International) · Kurt Keville (MIT / Posit Foundation)

Flat128 is an optional wide fixed-point extension for RISC-V providing bitwise-reproducible accumulation across threads and NUMA nodes without per-operation rounding.

## Files
- `flatdot.c` — dot product reproducibility demonstration
- `Mult128.c` — 128×128-bit multiply from sixteen RV64I MUL instructions
- `paper/main.tex` — LaTeX source

## Results
Tested on SG2042 Pioneer (64-hart RISC-V, Debian sid, GCC 13.3.0).
Flat128 accumulation is bitwise identical across K∈{64..4096} and T∈{1,4,16,64}.
IEEE 754 double diverges up to 8.9×10⁻¹⁶ at K=4096.

## Cite
See CITATION.cff or [Zenodo](https://doi.org/10.5281/zenodo.23246819).
