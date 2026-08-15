# flat128
Partial-product multiply demonstration for flat128 — a proposed
128-bit fixed-point arithmetic extension for RISC-V.

Pedagogical C implementation of 128-bit partial-product
multiplication using four 32-bit limbs, base 2^32.
Each limb product maps to one RV64I MUL instruction.

Companion code to:
> "Explicit Wide Accumulation: A Fixed-Point Extension Proposal for RISC-V"
> Kurt L. Keville, Paul D. Sherman (2026)

## What this demonstrates

flat128 treats a 128-bit register as a uniform two's complement
significand with compiler-managed binary point. No bits spent on
exponent overhead. Multiplication yields a 256-bit exact intermediate
product with explicit, architecturally-specified narrowing.

This code models 128×128-bit multiplication as sixteen 32×32 partial
products, each mapping to a single RV64I MUL instruction. No exotic
hardware required. The two-pass accumulation (partial products first,
carry ripple second) maps directly to a carry-save adder tree followed
by a carry-propagate adder in hardware — the classic CSA+CPA structure
that SiFive's P870 and similar out-of-order cores already exploit.

## Build

```bash
gcc -o Mult128 Mult128.c
```

Tested: GCC on Windows (mingw64), Linux, RISC-V (SG2042 Milk-V Pioneer).
No dependencies beyond the C standard library.

## Experimental Validation

`flatdot.c` demonstrates bitwise reproducibility of flat128 dot product
accumulation vs IEEE 754 double on real RISC-V silicon.

**Platform:** Milk-V Pioneer, SG2042, 64 THead C920 harts, 4 NUMA nodes,
128GB DDR4, Fedora kernel 6.1.37

**Build:**
```bash
gcc -O2 -o flatdot flatdot.c -lpthread
```

**Run:**
```bash
./flatdot <K> <nthreads>   # K=elements, nthreads=1..64
```

**Results (1 thread vs 64 threads, all 4 NUMA nodes):**

| K    | flat128 hex               | IEEE 754 double delta |
|------|---------------------------|-----------------------|
| 64   | `002b4c60`                | 0                     |
| 256  | `0a33d180`                | ~2.7e+14              |
| 1024 | `00000002 83394600`       | ~7e+15                |
| 4096 | `000000a0 33851800`       | ~4.5e+17              |

flat128 hex is identical across all thread counts and NUMA nodes.
IEEE 754 double diverges monotonically with K.

## Usage
```
Mult128 [<multiplier> <multiplicand> [<binpoint>]]
```

No args runs the self-test.

## Money shots

### Self-test: (2^32 + 1)^2 = 2^64 + 2·2^32 + 1

```
Mult128
m[]= 1 1 0 0
n[]= 1 1 0 0
s[]= 1 2 1 0 0 0 0 0
```

Verify by inspection: (x+1)^2 = x^2 + 2x + 1 with x = 2^32.
Three non-zero limbs. Exact. No rounding.

### 2^32 × 2^32 = 2^64 — carry into limb 2

```
Mult128 4294967296 4294967296
s[]= 0 0 1 0 0 0 0 0
```

### Near-UINT64_MAX × 2 — carry propagation across limb boundary

```
Mult128 18446744073709551615 2
m[]= 4294967295 4294967295 0 0
n[]= 2 0 0 0
s[]= 4294967294 4294967295 1 0 0 0 0 0
value=18446744073709551614
```

(2^64 - 1) × 2 = 2^65 - 2. Limb 2 = 1 is the carry propagating
visibly across the limb boundary. value= shows low 64 bits only —
truncation expected and documented.

## The flat128 argument in one picture

```
float32:  [S][ exponent  8b ][      mantissa 23b      ]
float128: [S][  exponent 15b ][          mantissa 112b          ]
flat128:  [               significand 128b                      ]
                                    ^
                                    binary point (compiler-managed)
```

The exponent field purchases dynamic range at the cost of significand
bits. flat128 eliminates it. ~38.5 decimal digits of uniform precision
at the chosen binary-point position.

## ISA mapping (RV64I)

Each 32×32→64 limb product is one instruction:

```asm
MUL   t0, m_i, n_j    ; full 64-bit product of two 32-bit operands
SRLI  t1, t0,  32     ; high 32 bits = carry limb
                      ; low  32 bits of t0 = result limb
```

Four 64-bit operand pairs reduce to four MUL/MULHU pairs for a full
128×128→256 product. Base RV64I suffices today. Hardware single-cycle
extension is incremental.

## History

- 20260807 pds  initial cut, 128-bit extension, base 2^32 limbs, flat128 demo
- 20260808 pds  (prod_t) cast on multiply inner loop, Windows LLP64 fix
- 20260814 pds  flatdot.c: bitwise reproducibility test, Pioneer silicon results

## License

MIT
