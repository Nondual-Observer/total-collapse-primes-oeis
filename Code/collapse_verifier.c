// collapse_verifier.c -- exact C verifier for the total-collapse prime sequences (Cards 00-18).
// Family id = card number:
//   0 (x^2-2)^2 ->{1} (A401173)   1 (z+2)(z^3+2z-1)+z ->{-2}   2 V3^2 ->{0}   3 (x^2-2)^4 ->{1}
//   4 V5^2 ->{0,1}   5 V8^2 ->{1}   6 V9^2 ->{0}   7 V7^2 ->{0,1}   8 (z+2)P251(z)+z ->{-2}
//   9 V6^2 -> value 4 (unique cycle through 4)   10 x^2(x-2)^2 -> fixed points
//   11..17: K_n(u) = (u^3-(n^2+n+1)u)/(n(n+1)) -> {-1,0,1} for n = 2,4,1,5,6,7,8
// Rejection occurs only when an extraneous periodic cycle is certified via unbounded Brent cycle
// detection; acceptance occurs only after a complete functional-graph traversal of GF(p).
// Build:  clang -O3 -o collapse_verifier collapse_verifier.c -lpthread     (or gcc)
// Usage:
//   ./collapse_verifier scan <fam> <pmin> <pmax> <threads>   list all terms in [pmin,pmax]
//   ./collapse_verifier check <fam> <p>                      full traversal for one prime
//   ./collapse_verifier analyze <fam> <p>                    all cycles with lengths and basins
//   ./collapse_verifier card18 <nmin> <nmax> <pcap>          Card 18 values (search p <= pcap)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <stdatomic.h>
#include <time.h>

typedef uint64_t u64;
typedef unsigned __int128 u128;

typedef struct { u64 p; u64 c[16]; u64 invB; u64 A; int fam; } Ctx;

static inline u64 mulm(u64 a, u64 b, u64 p) {
    if (p < (1ULL << 32)) return (a * b) % p;
    return (u64)(((u128)a * b) % p);
}
static inline u64 addm(u64 a, u64 b, u64 p) { u64 s = a + b; return s >= p ? s - p : s; }
static inline u64 subm(u64 a, u64 b, u64 p) { return a >= b ? a - b : a + p - b; }

static u64 powm(u64 a, u64 e, u64 p) { u64 r = 1 % p; a %= p; while (e) { if (e & 1) r = mulm(r, a, p); a = mulm(a, a, p); e >>= 1; } return r; }

// K_n parameters per family
static const int KN[18] = {0,0,0,0,0,0,0,0,0,0,0, 2,4,1,5,6,7,8};

static void ctx_init(Ctx *C, int fam, u64 p, int n_override) {
    C->p = p; C->fam = fam;
    for (int i = 0; i < 16; i++) C->c[i] = (u64)i % p;
    int n = n_override ? n_override : KN[fam];
    if (n) {
        u64 B = ((u64)n * (n + 1)) % p, A = ((u64)n * n + n + 1) % p;
        C->A = A; C->invB = B ? powm(B, p - 2, p) : 0;
    }
}

static inline u64 V3(u64 x, const Ctx *C) { u64 p = C->p; u64 x2 = mulm(x, x, p); return mulm(x, subm(x2, C->c[3], p), p); }

static inline u64 fmap(u64 x, const Ctx *C) {
    u64 p = C->p; const u64 *c = C->c;
    switch (C->fam) {
    case 0: { u64 y = subm(mulm(x, x, p), c[2], p); return mulm(y, y, p); }
    case 1: { u64 z2 = mulm(x, x, p); u64 t = subm(mulm(addm(z2, c[2], p), x, p), c[1], p); return addm(mulm(addm(x, c[2], p), t, p), x, p); }
    case 2: { u64 v = V3(x, C); return mulm(v, v, p); }
    case 3: { u64 y = subm(mulm(x, x, p), c[2], p); y = mulm(y, y, p); return mulm(y, y, p); }
    case 4: { u64 x2 = mulm(x, x, p); u64 t = addm(subm(mulm(x2, x2, p), mulm(c[5], x2, p), p), c[5], p); u64 v = mulm(x, t, p); return mulm(v, v, p); }
    case 5: { u64 y = subm(mulm(x, x, p), c[2], p); y = subm(mulm(y, y, p), c[2], p); y = subm(mulm(y, y, p), c[2], p); return mulm(y, y, p); }
    case 6: { u64 v = V3(V3(x, C), C); return mulm(v, v, p); }
    case 7: { u64 x2 = mulm(x, x, p); u64 t = subm(x2, c[7], p); t = addm(mulm(t, x2, p), c[14], p); t = subm(mulm(t, x2, p), c[7], p); u64 v = mulm(x, t, p); return mulm(v, v, p); }
    case 8: { // (z+2)*P251(z) + z, P = z^7+5z^6+6z^5-2z^4-4z^3+2z^2+9z-1
        u64 t = addm(x, c[5], p);
        t = addm(mulm(t, x, p), c[6], p);
        t = subm(mulm(t, x, p), c[2], p);
        t = subm(mulm(t, x, p), c[4], p);
        t = addm(mulm(t, x, p), c[2], p);
        t = addm(mulm(t, x, p), c[9], p);
        t = subm(mulm(t, x, p), c[1], p);
        return addm(mulm(addm(x, c[2], p), t, p), x, p); }
    case 9: { u64 v = V3(x, C); v = subm(mulm(v, v, p), c[2], p); return mulm(v, v, p); }
    case 10: { u64 y = mulm(x, subm(x, c[2], p), p); return mulm(y, y, p); }
    default: { u64 t = mulm(x, subm(mulm(x, x, p), C->A, p), p); return mulm(t, C->invB, p); }
    }
}

// is cycle through point c (length lam) acceptable?  For fam 9: contains 4.
static int cycle_good(u64 cpt, u64 lam, const Ctx *C) {
    u64 p = C->p;
    switch (C->fam) {
    case 0: case 3: case 5: return lam == 1 && cpt == 1 % p;
    case 1: case 8: return lam == 1 && cpt == (p - 2 % p) % p;
    case 2: case 6: return lam == 1 && cpt == 0;
    case 4: case 7: return lam == 1 && (cpt == 0 || cpt == 1 % p);
    case 9: { u64 four = 4 % p, x = cpt; for (u64 i = 0; i < lam; i++) { if (x == four) return 1; x = fmap(x, C); } return 0; }
    case 10: return lam == 1;
    default: return (lam == 1 && cpt == 0) || (lam == 2 && (cpt == 1 % p || cpt == p - 1));
    }
}

static int admissible(int fam, u64 p) {
    if (fam >= 11) { int n = KN[fam]; if (p < 5) return 0; if (((u64)n * (n + 1)) % p == 0) return 0; }
    return 1;
}
// Proven necessary conditions (used only when USE_FILTERS): Dickson permutation barrier (Thm. 3.6, families 4 and 7).
static int USE_FILTERS = 1;
static int filtered_out(int fam, u64 p) {
    if (!USE_FILTERS) return 0;
    if ((fam == 4 || fam == 7) && p > 7) {
        u64 m = fam == 4 ? 5 : 7; u64 r = p % m;
        int perm = (r * r) % m != 1; // gcd(m, p^2-1) = 1
        if (perm && p % 4 == 3) return 1;
    }
    return 0;
}

// Brent cycle finding from x0: returns a point on the cycle and its length.
static u64 brent(u64 x0, const Ctx *C, u64 *lam_out) {
    u64 power = 1, lam = 1, tort = x0, hare = fmap(x0, C);
    while (tort != hare) {
        if (power == lam) { tort = hare; power <<= 1; lam = 0; }
        hare = fmap(hare, C); lam++;
    }
    *lam_out = lam; return hare;
}

// full traversal; returns 1 if collapse. mem must be zeroed array of p bytes.
static int full_traversal(const Ctx *C, uint8_t *st, u64 *ncyc_out) {
    u64 p = C->p, ncyc = 0;
    for (u64 x = 0; x < p; x++) {
        if (st[x]) continue;
        u64 c = x;
        while (st[c] == 0) { st[c] = 1; c = fmap(c, C); }
        if (st[c] == 1) {
            u64 lam = 1, y = fmap(c, C); while (y != c) { y = fmap(y, C); lam++; }
            ncyc++;
            if (!cycle_good(c, lam, C)) { *ncyc_out = ncyc; return 0; }
            if (C->fam == 9 && ncyc > 1) { *ncyc_out = ncyc; return 0; }
        }
        c = x; while (st[c] == 1) { st[c] = 2; c = fmap(c, C); }
    }
    *ncyc_out = ncyc; return 1;
}

static atomic_ullong phase2_count;
static const int NPROBE = 24;

static int test_prime(int fam, u64 p, uint8_t **buf, u64 *bufcap) {
    if (!admissible(fam, p)) return 0;
    if (filtered_out(fam, p)) return 0;
    Ctx C; ctx_init(&C, fam, p, 0);
    // phase 1: probes (rigorous rejection)
    u64 s = 0x9E3779B97F4A7C15ULL ^ p;
    for (int i = 0; i < NPROBE && (u64)i < p; i++) {
        u64 x0;
        if (i < 4) x0 = (u64)(i + 2) % p; else { s ^= s << 13; s ^= s >> 7; s ^= s << 17; x0 = s % p; }
        u64 lam, cpt = brent(x0, &C, &lam);
        if (!cycle_good(cpt, lam, &C)) return 0;
    }
    // phase 2: exact
    atomic_fetch_add(&phase2_count, 1);
    if (*bufcap < p) { free(*buf); *buf = malloc(p); *bufcap = p; }
    memset(*buf, 0, p);
    u64 nc; return full_traversal(&C, *buf, &nc);
}

// ---------------- scan driver ----------------
static u64 *primes; static u64 nprimes; static atomic_ullong next_idx; static int g_fam;
static u64 *found; static atomic_ullong nfound;

static void *worker(void *arg) {
    (void)arg; uint8_t *buf = NULL; u64 cap = 0;
    for (;;) {
        u64 i0 = atomic_fetch_add(&next_idx, 256);
        if (i0 >= nprimes) break;
        u64 i1 = i0 + 256 < nprimes ? i0 + 256 : nprimes;
        for (u64 i = i0; i < i1; i++) {
            if (test_prime(g_fam, primes[i], &buf, &cap)) { u64 k = atomic_fetch_add(&nfound, 1); found[k] = primes[i]; }
        }
        // free very large buffers between chunks to cap memory
        if (cap > 600000000ULL) { free(buf); buf = NULL; cap = 0; }
    }
    free(buf); return NULL;
}

static int cmpu(const void *a, const void *b) { u64 x = *(const u64 *)a, y = *(const u64 *)b; return x < y ? -1 : x > y; }

static void sieve_range(u64 pmin, u64 pmax) {
    u64 n = pmax + 1; uint8_t *comp = calloc(n / 2 + 1, 1); // odd numbers
    for (u64 i = 3; i * i <= pmax; i += 2) if (!comp[i / 2]) for (u64 j = i * i; j <= pmax; j += 2 * i) comp[j / 2] = 1;
    u64 cnt = 0; if (pmin <= 2 && pmax >= 2) cnt++;
    for (u64 i = 3; i <= pmax; i += 2) if (!comp[i / 2] && i >= pmin) cnt++;
    primes = malloc(sizeof(u64) * (cnt + 1)); nprimes = 0;
    if (pmin <= 2 && pmax >= 2) primes[nprimes++] = 2;
    for (u64 i = 3; i <= pmax; i += 2) if (!comp[i / 2] && i >= pmin) primes[nprimes++] = i;
    free(comp);
}

static int isprime_u64(u64 n) { if (n < 2) return 0; for (u64 d = 2; d * d <= n; d++) if (n % d == 0) return 0; return 1; }

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s <scan|check|analyze|card18> ...\n", argv[0]); return 1; }
    if (!strcmp(argv[1], "scan") && argc >= 6) {
        g_fam = atoi(argv[2]); u64 pmin = strtoull(argv[3], 0, 10), pmax = strtoull(argv[4], 0, 10); int nt = atoi(argv[5]);
        struct timespec t0, t1; clock_gettime(CLOCK_MONOTONIC, &t0);
        sieve_range(pmin, pmax); found = malloc(sizeof(u64) * 100000);
        pthread_t th[64]; for (int i = 0; i < nt; i++) pthread_create(&th[i], 0, worker, 0);
        for (int i = 0; i < nt; i++) pthread_join(th[i], 0);
        qsort(found, nfound, sizeof(u64), cmpu);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        printf("FAM %02d range [%llu,%llu] primes=%llu phase2=%llu time=%.1fs\nTERMS:", g_fam, (unsigned long long)pmin, (unsigned long long)pmax,
               (unsigned long long)nprimes, (unsigned long long)phase2_count, (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec));
        for (u64 i = 0; i < nfound; i++) printf(" %llu", (unsigned long long)found[i]);
        printf("\n"); return 0;
    }
    if (!strcmp(argv[1], "check") && argc >= 4) {
        int fam = atoi(argv[2]); u64 p = strtoull(argv[3], 0, 10);
        int pr = (p < 4000000000000ULL) ? isprime_u64(p) : -1;
        if (!admissible(fam, p)) { printf("FAM %02d p=%llu prime=%d NOT ADMISSIBLE\n", fam, (unsigned long long)p, pr); return 0; }
        Ctx C; ctx_init(&C, fam, p, 0);
        uint8_t *st = calloc(p, 1); if (!st) { printf("alloc fail\n"); return 1; }
        u64 nc; int ok = full_traversal(&C, st, &nc);
        printf("FAM %02d p=%llu prime=%d COLLAPSE=%s cycles_seen=%llu\n", fam, (unsigned long long)p, pr, ok ? "YES" : "NO", (unsigned long long)nc);
        return 0;
    }
    if (!strcmp(argv[1], "analyze") && argc >= 4) {
        int fam = atoi(argv[2]); u64 p = strtoull(argv[3], 0, 10);
        Ctx C; ctx_init(&C, fam, p, 0);
        uint32_t *lab = calloc(p, 4); const uint32_t ONP = 0xFFFFFFFFu; uint32_t ncyc = 0;
        u64 cap = 1024; u64 *clen = malloc(8 * cap), *cmin = malloc(8 * cap), *bas = calloc(cap, 8);
        for (u64 x = 0; x < p; x++) {
            if (lab[x]) continue;
            u64 c = x; while (lab[c] == 0) { lab[c] = ONP; c = fmap(c, &C); }
            if (lab[c] == ONP) {
                ncyc++; if (ncyc >= cap) { cap *= 2; clen = realloc(clen, 8 * cap); cmin = realloc(cmin, 8 * cap); bas = realloc(bas, 8 * cap); memset(bas + cap / 2, 0, 8 * cap / 2); }
                u64 lam = 0, mn = c, y = c; do { lab[y] = ncyc; if (y < mn) mn = y; y = fmap(y, &C); lam++; } while (y != c);
                clen[ncyc] = lam; cmin[ncyc] = mn;
            }
            uint32_t id = lab[c]; c = x; while (lab[c] == ONP) { lab[c] = id; c = fmap(c, &C); }
        }
        for (u64 x = 0; x < p; x++) bas[lab[x]]++;
        printf("FAM %02d p=%llu cycles=%u\n", fam, (unsigned long long)p, ncyc);
        for (uint32_t i = 1; i <= ncyc; i++) printf("  cycle %u: len=%llu min=%llu basin=%llu (%.6f%%)\n", i, (unsigned long long)clen[i], (unsigned long long)cmin[i], (unsigned long long)bas[i], 100.0 * bas[i] / p);
        return 0;
    }
    if (!strcmp(argv[1], "card18") && argc >= 5) {
        int nmin = atoi(argv[2]), nmax = atoi(argv[3]); u64 pcap = strtoull(argv[4], 0, 10);
        sieve_range(2, pcap); uint8_t *buf = NULL; u64 bcap = 0;
        for (int n = nmin; n <= nmax; n++) {
            u64 ans = 0;
            for (u64 i = 0; i < nprimes; i++) {
                u64 p = primes[i]; if (p < 5 || p <= (u64)n + 1) continue;
                Ctx C; ctx_init(&C, 11, p, n);
                int bad = 0; u64 s = p * 2654435761ULL;
                for (int k = 0; k < NPROBE && !bad; k++) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; u64 lam, cpt = brent(k < 3 ? (u64)(k + 2) % p : s % p, &C, &lam); if (!cycle_good(cpt, lam, &C)) bad = 1; }
                if (bad) continue;
                if (bcap < p) { free(buf); buf = malloc(p); bcap = p; }
                memset(buf, 0, p); u64 nc;
                if (full_traversal(&C, buf, &nc)) { ans = p; break; }
            }
            printf("a(%d) = %llu%s\n", n, (unsigned long long)ans, ans ? "" : " (none below cap)");
            fflush(stdout);
        }
        return 0;
    }
    return 1;
}
