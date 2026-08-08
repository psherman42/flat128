# flat128
Partial-product multiply demonstration for flat128 -- a proposed 128-bit fixed-point arithmetic extension for RISC-V.

Pedagogical C implementation of 128-bit partial-product
multiplication using four 32-bit limbs, base 2^32.
Each limb product maps to one RV64I MUL instruction.
