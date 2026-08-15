/*
 * flatdot.c -- bitwise reproducibility test for flat128 accumulation
 *
 * Demonstrates exact, reproducible dot product accumulation using
 * flat128 (128-bit two's complement fixed-point) vs IEEE 754 double.
 *
 * Usage: ./flatdot <K> <nthreads>
 *   K        : number of elements in dot product (try 64,256,1024,4096)
 *   nthreads : number of parallel threads (try 1,4,16,64)
 *
 * Build (Pioneer / Fedora):
 *   gcc -O2 -o flatdot flatdot.c -lpthread
 *
 * Companion code to:
 *   "Explicit Wide Accumulation: A Fixed-Point Extension Proposal for RISC-V"
 *   Keville, Sherman -- targeting SC26 / CoNGA'26
 *
 * Platform: SG2042 Milk-V Pioneer, Fedora userland
 *           64 C920 harts, 4 NUMA nodes, 128GB DDR4
 *
 * Paul Sherman / Kurt Keville -- August 2026
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>

/* ------------------------------------------------------------------ */
/* flat128: 128-bit two's complement significand, four 32-bit limbs    */
/* limb[0] = least significant, limb[3] = most significant            */
/* binary point convention: fixed at bit 64 (limb[2] boundary)        */
/* ------------------------------------------------------------------ */

#define N_LIMBS 4

typedef struct {
    uint32_t limb[N_LIMBS];
} flat128_t;

typedef struct {
    uint32_t limb[N_LIMBS * 2];   /* 256-bit exact intermediate product */
} flat256_t;

/* ------------------------------------------------------------------ */
/* flat128 addition -- exact when operands share scale                 */
/* two-pass: accumulate partial sums, then ripple carry                */
/* mirrors CSA+CPA hardware described in Mult128.c                     */
/* ------------------------------------------------------------------ */

static flat128_t
f128_add(flat128_t a, flat128_t b)
{
    flat128_t r;
    uint64_t carry = 0;
    int i;

    for (i = 0; i < N_LIMBS; i++) {
        uint64_t s = (uint64_t)a.limb[i] + (uint64_t)b.limb[i] + carry;
        r.limb[i] = (uint32_t)(s & 0xFFFFFFFFULL);
        carry = s >> 32;
    }
    /* carry out silently discarded -- overflow is explicit concern     */
    /* for accumulation depths used here (K<=4096) no overflow occurs   */
    return r;
}

/* ------------------------------------------------------------------ */
/* flat128 multiply -- 128x128 -> 256 exact intermediate               */
/* sixteen 32x32 partial products, each maps to one RV64I MUL          */
/* (see Mult128.c for full derivation and carry commentary)            */
/* narrow to flat128 by dropping low 128 bits (explicit, not hidden)   */
/* ------------------------------------------------------------------ */

static flat128_t
f128_mul_narrow(flat128_t a, flat128_t b)
{
    flat256_t p;
    uint64_t carry;
    int i, j;

    memset(&p, 0, sizeof(p));

    for (i = 0; i < N_LIMBS; i++) {
        carry = 0;
        for (j = 0; j < N_LIMBS; j++) {
            uint64_t prod = (uint64_t)a.limb[i] * (uint64_t)b.limb[j];
            uint64_t s    = (uint64_t)p.limb[i+j] + prod + carry;
            p.limb[i+j]   = (uint32_t)(s & 0xFFFFFFFFULL);
            carry          = s >> 32;
        }
        p.limb[i + N_LIMBS] += (uint32_t)carry;
    }

    /* explicit narrowing: take upper 128 bits of 256-bit product      */
    /* binary point shifts by 64 bits -- caller must account for scale */
    flat128_t r;
    for (i = 0; i < N_LIMBS; i++)
        r.limb[i] = p.limb[i + N_LIMBS];
    return r;
}

/* ------------------------------------------------------------------ */
/* deterministic input generation -- same seed = same values every run */
/* no randomness: reproducibility is the variable, not the inputs      */
/* ------------------------------------------------------------------ */

static flat128_t
gen_a(int idx)
{
    flat128_t v;
    /* simple deterministic pattern: spread index across limbs         */
    v.limb[0] = (uint32_t)(idx + 1);
    v.limb[1] = (uint32_t)(idx * 3 + 7);
    v.limb[2] = (uint32_t)(idx * 5 + 13);
    v.limb[3] = 0;
    return v;
}

static flat128_t
gen_b(int idx)
{
    flat128_t v;
    v.limb[0] = (uint32_t)(idx * 2 + 1);
    v.limb[1] = (uint32_t)(idx * 4 + 3);
    v.limb[2] = (uint32_t)(idx * 6 + 11);
    v.limb[3] = 0;
    return v;
}

/* ------------------------------------------------------------------ */
/* thread work: partial dot product over assigned slice of K elements  */
/* ------------------------------------------------------------------ */

typedef struct {
    int     tid;
    int     K;
    int     nthreads;
    flat128_t result_f128;
    double    result_f64;
} thread_arg_t;

static void *
dot_worker(void *arg)
{
    thread_arg_t *ta = (thread_arg_t *)arg;
    int start = (ta->tid * ta->K) / ta->nthreads;
    int end   = ((ta->tid + 1) * ta->K) / ta->nthreads;
    int i;

    flat128_t acc_f128;
    double    acc_f64 = 0.0;
    memset(&acc_f128, 0, sizeof(acc_f128));

    for (i = start; i < end; i++) {
        flat128_t a = gen_a(i);
        flat128_t b = gen_b(i);

        /* flat128: exact multiply then exact add                      */
        flat128_t prod = f128_mul_narrow(a, b);
        acc_f128 = f128_add(acc_f128, prod);

        /* IEEE 754 double: subject to reordering non-associativity    */
        double da = (double)a.limb[0] + (double)a.limb[1] * 4294967296.0;
        double db = (double)b.limb[0] + (double)b.limb[1] * 4294967296.0;
        acc_f64 += da * db;
    }

    ta->result_f128 = acc_f128;
    ta->result_f64  = acc_f64;
    return NULL;
}

/* ------------------------------------------------------------------ */
/* print flat128 as hex, most-significant limb first                   */
/* ------------------------------------------------------------------ */

static void
print_f128(const char *label, flat128_t v)
{
    printf("%s: %08x %08x %08x %08x\n",
        label,
        v.limb[3], v.limb[2], v.limb[1], v.limb[0]);
}

/* ------------------------------------------------------------------ */
/* main                                                                 */
/* ------------------------------------------------------------------ */

int
main(int argc, char *argv[])
{
    if (argc < 3) {
        fprintf(stderr, "usage: %s <K> <nthreads>\n", argv[0]);
        fprintf(stderr, "  K       : dot product length (64,256,1024,4096)\n");
        fprintf(stderr, "  nthreads: parallel threads  (1,4,16,64)\n");
        return 1;
    }

    int K        = (int)strtol(argv[1], NULL, 10);
    int nthreads = (int)strtol(argv[2], NULL, 10);

    if (K < 1 || nthreads < 1 || nthreads > 64) {
        fprintf(stderr, "K>=1, 1<=nthreads<=64\n");
        return 1;
    }

    printf("flatdot: K=%d nthreads=%d\n", K, nthreads);
    printf("Platform: SG2042 Pioneer / Fedora\n\n");

    pthread_t      *threads = malloc(nthreads * sizeof(pthread_t));
    thread_arg_t   *args    = malloc(nthreads * sizeof(thread_arg_t));
    int i;

    for (i = 0; i < nthreads; i++) {
        args[i].tid      = i;
        args[i].K        = K;
        args[i].nthreads = nthreads;
        memset(&args[i].result_f128, 0, sizeof(flat128_t));
        args[i].result_f64 = 0.0;
        pthread_create(&threads[i], NULL, dot_worker, &args[i]);
    }

    /* reduce partial results across threads                           */
    flat128_t total_f128;
    double    total_f64 = 0.0;
    memset(&total_f128, 0, sizeof(total_f128));

    for (i = 0; i < nthreads; i++) {
        pthread_join(threads[i], NULL);
        total_f128 = f128_add(total_f128, args[i].result_f128);
        total_f64 += args[i].result_f64;
    }

    /* results -- compare hex across runs and thread counts            */
    print_f128("flat128 result", total_f128);
    printf("double  result: %.20e\n", total_f64);
    printf("\nflat128 hex should be IDENTICAL across all runs,\n");
    printf("thread counts, and NUMA nodes on Pioneer.\n");
    printf("double  hex will diverge under thread reordering.\n");

    free(threads);
    free(args);
    return 0;
}
