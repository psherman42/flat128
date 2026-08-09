/*
 * Mult128.c -- flat128 partial-product multiply demonstration
 *
 * Models 128-bit fixed-point multiplication as four 32-bit limbs,
 * base 2^32.  Each limb multiply maps directly to one MUL/MULHU
 * pair on RV64I.  The 256-bit exact product is assembled in p[][]
 * and summed into s[] with explicit carry propagation.
 *
 * Binary point (binpoint) is caller-managed -- no bits spent on
 * exponent overhead.  This is the flat128 accumulation argument
 * made executable.
 *
 * Companion to: "Depth or Breadth? RV128 and Next Generation
 * Arithmetic" -- K.L. Keville, P.D. Sherman
 *
 * 20260807 pds  initial cut, 128-bit by extension, base 2^32 limbs, flat128 demo
 * 20260808 pds  (prod_t) cast on multiply inner loop, Windows LLP64 fix
 */

#include <stdlib.h>  /* for atoi(), strtoull() */
#include <stdio.h>  /* for printf() */
#include <stdint.h>  /* for uint32_t, uint64_t, __uint128_t */
#include <inttypes.h>  /* for PRIu64 */
#include <time.h>  /* for clock() */

#define M_LIMBS 4   /* 128-bit multiplicand: 128 / 32 bits per limb */
#define N_LIMBS 4   /* 128-bit multiplier:   128 / 32 bits per limb */

#define BASE      0x100000000ULL  /* 2^32, one limb */

typedef uint32_t limb_t;    /* one 32-bit limb */
typedef uint64_t prod_t;    /* intermediate product, won't overflow */

int value_to_limbs(__uint128_t value, limb_t *limbs)
{
    int i;
    i = 0;
    while (value)
    {
        limbs[i++] = value % BASE;  // 1s place
        value /= BASE;  // 10s place
    }
    return i;
}

/*
 * On RV64I, a 32x32->64 limb product is one instruction:
 *
 *   MUL   t0, m_i, n_j    -- full 64-bit product of two 32-bit operands
 *   SRLI  t1, t0,  32     -- high 32 bits = carry limb
 *                         -- low  32 bits of t0 = result limb
 *
 * MULHU is not needed here: for 32-bit operands zero-extended into
 * 64-bit registers, the upper 64 bits of a 64x64 product are zero.
 *
 * This code uses 16 x 32x32 products (4 limbs x 4 limbs) for
 * transparency.  An RV64 hardware implementation can instead use
 * two 64-bit limbs, reducing to 4 x 64x64 products via MUL/MULHU
 * on full 64-bit registers -- the article's hardware claim.
 * Both decompositions are mathematically equivalent.
 */
void multiply(limb_t *m, limb_t *n, prod_t p[M_LIMBS+N_LIMBS][N_LIMBS])
{
    int i, mm;
    int j, nn;
    mm=M_LIMBS;
    nn=N_LIMBS;
    for (j=0; j<nn; j++)
    {
        for (i=0; i<mm; i++)
        {
            p[i+j][j] += ((prod_t)m[i] * (prod_t)n[j]);
            p[i+j+1][j] += p[i+j][j] / BASE;  /* 20260807 pds. was assignment, need accumulation */
            p[i+j][j] %= BASE;
        }
    }
}

void add(prod_t p[M_LIMBS+N_LIMBS][N_LIMBS], prod_t s[M_LIMBS+N_LIMBS])
{
    /*
     * Two-pass accumulation with full ripple carry.
     *
     * Naive interleaved carry:
     *     s[i+1] += s[i] / BASE;
     *     s[i]   %= BASE;
     * is fragile -- if s[i] accumulates N_LIMBS partial products
     * before the carry fires, a single carry step may not be enough.
     * Multiple partial products can sum beyond 2*BASE, requiring
     * carry to ripple more than one limb at a time.
     *
     * Solution: accumulate all partial products first, then ripple.
     * Pass 1 -- accumulate columns into s[], no carry yet.
     * Pass 2 -- ripple carry cleanly from limb 0 to limb M+N-1.
     *
     * RTL/circuit note: this is the classic carry-propagate adder
     * (CPA) stage following a carry-save adder (CSA) tree.
     * Pass 1 represents the carry-delayed accumulation that hardware
     * can implement with a CSA tree; Pass 2 represents the final CPA.
     * Efficient hardware resolves Pass 2 in O(log n) gate delays
     * via carry-lookahead or Han-Carlson rather than O(n) ripple.
     * SiFive's P870 out-of-order core already issues the independent
     * limb MUL/MULHU pairs in parallel without software intervention --
     * exactly the division of labor this two-pass model respects.
     *
     * Single-pass is phenomenologically impossible in software --
     * carry out of limb i depends on carry into limb i (serial
     * dependency chain).  Hardware escapes via carry-lookahead:
     * generate/propagate signals computed in parallel, O(log n)
     * gate delays.  This ripple is the honest software model.
     *
     * The generate/propagate signals G[i]/P[i] are locally independent
     * and an out-of-order engine could exploit this -- but the carry
     * dependency chain C[i] -> C[i+1] is a true data hazard that OOO
     * rename cannot dissolve.  Expressing carry-lookahead explicitly
     * in software requires manual assembly restructuring: non-portable,
     * compiler-hostile, fragile across microarchitectures, and no
     * faster than what the OOO engine already does with the ripple.
     * Let hardware do what hardware does best.  Software expresses the
     * algorithm honestly; microarchitecture optimizes the execution.
     * Conflating the two layers opens cans of worms -- barriers,
     * hazards, flush vulnerabilities -- for no gain.
     * Simplicity wins.  -- A. Shugart, paraphrased.
     */
    int i, mm;
    int j, nn;
    mm=M_LIMBS;
    nn=N_LIMBS;

    /* Pass 1 -- accumulate all partial products, no carry yet */
    for (i=0; i<(mm+nn); i++)
    {
        for (j=0; j<nn; j++)
        {
            s[i] += p[i][j];
        }
    }

    /* Pass 2 -- ripple carry cleanly limb 0 to limb M+N-1 */
    for (i=0; i<(mm+nn-1); i++)
    {
        s[i+1] += s[i] / BASE;
        s[i] %= BASE;
    }
}

/* demo: returns low 64 bits of 256-bit result */
uint64_t limbs_to_value(prod_t *limbs, int binpoint)
{
    uint64_t value;
    __uint128_t place;
    int k;
    place = 1;
    value = 0;
    for (k=0; k<(M_LIMBS+N_LIMBS); k++)
    {
        if (k >= binpoint)
        {
            value += (limbs[k] * place);
            place *= BASE;
        }
    }
    return value;
}

int main(int argc, char *argv[])
{
    int rc;
    int i, j, k;
    uint64_t M;   /* input value, 64-bit portable -- LP64 vs LLP64 okay */
    uint64_t N;   /* input value, 64-bit portable -- LP64 vs LLP64 okay */
    /* 4 limbs × 4 limbs = 8-limb (256-bit) exact product */
    limb_t m[M_LIMBS];    /* 32-bit limbs of M */
    limb_t n[N_LIMBS];    /* 32-bit limbs of N */
    uint64_t p[M_LIMBS+N_LIMBS][N_LIMBS];
    uint64_t s[M_LIMBS+N_LIMBS];
    uint64_t value;
    int binpoint;

    clock_t t0 = clock();

    for (i=0; i<M_LIMBS; i++)
    {
        m[i] = 0;
    }
    for (j=0; j<N_LIMBS; j++)
    {
        n[j] = 0;
    }
    for (j=0; j<N_LIMBS; j++)
    {
        for (i=0; i<(M_LIMBS+N_LIMBS); i++)
        {
            p[i][j] = 0;
        }
    }
    for (k=0; k<(M_LIMBS+N_LIMBS); k++)
    {
        s[k] = 0;
    }

    if (argc == 1)
    {
        M = 4294967297ULL;   /* 2^32 + 1, two limbs: [1, 1, 0, 0] */
        N = 4294967297ULL;   /* 2^32 + 1, two limbs: [1, 1, 0, 0] */
        binpoint = 0;
        /* expected: s[]= 1 2 1 0 0 0 0 0                         */
        /*           i.e. (2^32+1)^2 = 2^64 + 2*2^32 + 1          */
        /*           value= low 64 bits only, truncation expected */
    }
    /* verified 20260808:                                              */
    /* Mult128              --> s[]= 1 2 1 0 0 0 0 0                   */
    /* Mult128 4294967296^2 --> s[]= 0 0 1 0 0 0 0 0                   */
    /* Mult128 2^64-1 2     --> s[]= 4294967294 4294967295 1 0 0 0 0 0 */
    else if (argc == 3)
    {
        M = (uint64_t)strtoull( argv[1], NULL, 10 );
        N = (uint64_t)strtoull( argv[2], NULL, 10 );
        binpoint = 0;
    }
    else if (argc == 4)
    {
        M = (uint64_t)strtoull( argv[1], NULL, 10 );
        N = (uint64_t)strtoull( argv[2], NULL, 10 );
        binpoint = atoi( argv[3] );
    }
    else
    {
        printf("\n\nUsage:");
        printf("\n\tMult128 [<multiplier> <multiplicand> [<binpoint>]]");
        printf("\n");
        printf("\n\tNo args: runs self-test: (2^32+1)^2, expected s[]= 1 2 1 0 0 0 0 0");
        printf("\n");
        printf("\n\tExamples:");
        printf("\n\t  Mult128 4294967297 4294967297   -- same as no-arg self-test");
        printf("\n\t  Mult128 4294967296 4294967296   -- 2^32 * 2^32 = 2^64, one limb carry");
        printf("\n\t  Mult128 18446744073709551615 2  -- near ULONG_MAX * 2, exercises high limbs");
        printf("\n");
        printf("\n\tWarnings print if M or N exceeds M_LIMBS or N_LIMBS limbs (4 limbs = 128 bits)");
        printf("\n\tbinpoint: limb offset of binary point, compiler-managed in real flat128 ISA");
        printf("\n\tOverflow warning demo requires >128-bit input; use __uint128_t");
        printf("\n\tinput parsing not shown -- left as exercise (or future binpoint demo)");
        printf("\n");
        return -1;
    }

    rc = value_to_limbs( M, m );
    /* NOTE: overflow unreachable from argv with uint64_t input --
     *       full __uint128_t argv parsing left as reader exercise
     */
    if (rc > M_LIMBS)
        printf("WARNING: value %d exceeds M_LIMBS %d\n", rc, M_LIMBS);

    rc = value_to_limbs( N, n );
    if (rc > N_LIMBS)
        printf("WARNING: value %d exceeds N_LIMBS %d\n", rc, N_LIMBS);

    /////////////////////////////////////
    //
    printf("\nm[]=");
    for (i=0; i<M_LIMBS; i++)
    {
        printf(" %u", m[i]);
    }
    printf("\nn[]=");
    for (i=0; i<N_LIMBS; i++) {
        printf(" %u", n[i]);
    }
    //
    /////////////////////////////////////

    multiply( m, n, p );

    /////////////////////////////////////
    //
    printf("\np[][]=\n");
    for (j=0; j<N_LIMBS; j++)
    {
        for (i=0; i<M_LIMBS+N_LIMBS; i++)
        {
            printf(" %lu", p[i][j]);
        }
        printf("\n");
    }
    //
    /////////////////////////////////////

    add(p, s);

    /////////////////////////////////////
    //
    printf("\ns[]=\n");
    for (k=0; k<M_LIMBS+N_LIMBS; k++)
    {
        printf(" %lu", s[k]);
    }
    //
    /////////////////////////////////////

    value = limbs_to_value(s, binpoint);

    printf("\nvalue=%" PRIu64, value);

    clock_t t1 = clock();
    printf("\nelapsed: %.6f sec\n", (double)(t1-t0)/CLOCKS_PER_SEC);

    return 0;
}
