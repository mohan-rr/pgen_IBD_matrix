#ifndef PGL_BASE_H
#define PGL_BASE_H

#define _POSIX_C_SOURCE 200809L

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__GNUC__) || defined(__clang__)
#define unlikely(expr) __builtin_expect(!!(expr), 0)
#define likely(expr) __builtin_expect(!!(expr), 1)
#else
#define unlikely(expr) (expr)
#define likely(expr) (expr)
#endif

typedef enum {
  kPglRetSuccess = 0,
  kPglRetSkipped = 1,
  kPglRetNomem = 2,
  kPglRetOpenFail = 3,
  kPglRetReadFail = 4,
  kPglRetWriteFail = 5,
  kPglRetMalformedInput = 6,
  kPglRetInconsistentInput = 7,
  kPglRetInvalidCmdline = 8,
  kPglRetThreadCreateFail = 9,
  kPglRetNetworkFail = 10,
  kPglRetVarRecordTooLarge = 11,
  kPglRetUnsupportedInstructions = 12,
  kPglRetDegenerateData = 13,
  kPglRetDecompressFail = 14,
  kPglRetRewindFail = 15,
  kPglRetGpuFail = 16,
  kPglRetInternalUse1 = 17,
  kPglRetSampleMajorBed = 32,
  kPglRetInternalError = 60,
  kPglRetNotYetSupported = 63,
  kPglRetLongLine = 126
} PglErr;

#define kBitsPerWord ((uint32_t)(8u * sizeof(uintptr_t)))
#define kBitsPerWordD2 (kBitsPerWord / 2u)

/* 01 repeated, and 10 repeated. Genotype words use the low bit of each pair. */
#define kMask5555 (~(uintptr_t)0 / 3u)
#define kMaskAAAA (kMask5555 << 1)

#define kPglDifflistGroupSize 64u
#define kPglVblockSize 65536u

static inline uint32_t BitCtToWordCt(uint32_t bit_ct) {
  return (bit_ct + kBitsPerWord - 1u) / kBitsPerWord;
}

static inline uint32_t NypCtToWordCt(uint32_t nyp_ct) {
  return (nyp_ct * 2u + kBitsPerWord - 1u) / kBitsPerWord;
}

static inline uint32_t NypCtToByteCt(uint32_t nyp_ct) {
  return (nyp_ct * 2u + 7u) / 8u;
}

static inline uint32_t DivUpU32(uint32_t numer, uint32_t denom) {
  return (numer + denom - 1u) / denom;
}

static inline uint32_t BytesToRepresentNzU32(uint32_t uii) {
  if (uii == 0u) {
    return 1u;
  }
  return 1u + ((31u - (uint32_t)__builtin_clz(uii)) / 8u);
}

static inline uint32_t PopcountW(uintptr_t ww) {
#if UINTPTR_MAX > 0xffffffffu
  return (uint32_t)__builtin_popcountll((unsigned long long)ww);
#else
  return (uint32_t)__builtin_popcount((unsigned)ww);
#endif
}

static inline uint32_t GetBit(const uintptr_t* bitarr, uint32_t idx) {
  return (uint32_t)((bitarr[idx / kBitsPerWord] >> (idx % kBitsPerWord)) & 1u);
}

static inline void SetBit(uintptr_t* bitarr, uint32_t idx) {
  bitarr[idx / kBitsPerWord] |= (uintptr_t)1 << (idx % kBitsPerWord);
}

static inline uint32_t GetGeno(const uintptr_t* genovec, uint32_t sample_idx) {
  uint32_t shift = 2u * (sample_idx % kBitsPerWordD2);
  return (uint32_t)((genovec[sample_idx / kBitsPerWordD2] >> shift) & 3u);
}

static inline void SetGeno(uintptr_t* genovec, uint32_t sample_idx, uint32_t geno) {
  uint32_t word = sample_idx / kBitsPerWordD2;
  uint32_t shift = 2u * (sample_idx % kBitsPerWordD2);
  uintptr_t mask = (uintptr_t)3 << shift;
  genovec[word] = (genovec[word] & ~mask) | ((uintptr_t)(geno & 3u) << shift);
}

static inline void ZeroTrailingNyps(uint32_t nyp_ct, uintptr_t* bitarr) {
  uint32_t bit_ct = nyp_ct * 2u;
  uint32_t rem = bit_ct % kBitsPerWord;
  if (rem) {
    bitarr[bit_ct / kBitsPerWord] &= (((uintptr_t)1 << rem) - (uintptr_t)1);
  }
}

/* Swap hom-ref (0) with hom-alt (2). Het and missing stay put. */
static inline uintptr_t InvertGenoWordUnsafe(uintptr_t geno_word) {
  return geno_word ^ ((~(geno_word << 1)) & kMaskAAAA);
}

static inline void GenovecInvertUnsafe(uint32_t sample_ct, uintptr_t* genovec) {
  uint32_t nword = NypCtToWordCt(sample_ct);
  uint32_t word_idx;
  for (word_idx = 0; word_idx < nword; word_idx++) {
    genovec[word_idx] = InvertGenoWordUnsafe(genovec[word_idx]);
  }
  ZeroTrailingNyps(sample_ct, genovec);
}

static inline void GenovecFill(uint32_t sample_ct, uint32_t geno, uintptr_t* genovec) {
  uint32_t nword = NypCtToWordCt(sample_ct);
  uintptr_t ww = (uintptr_t)(geno & 3u) * kMask5555;
  uint32_t word_idx;
  for (word_idx = 0; word_idx < nword; word_idx++) {
    genovec[word_idx] = ww;
  }
  ZeroTrailingNyps(sample_ct, genovec);
}

typedef struct Arena {
  unsigned char* base;
  uintptr_t cap;
  uintptr_t off;
} Arena;

static inline PglErr ArenaInit(Arena* arena, uintptr_t cap) {
  arena->base = 0;
  arena->cap = 0;
  arena->off = 0;
  if (cap < 64u) {
    cap = 64u;
  }
  arena->base = (unsigned char*)malloc(cap);
  if (unlikely(!arena->base)) {
    return kPglRetNomem;
  }
  arena->cap = cap;
  return kPglRetSuccess;
}

static inline void ArenaDestroy(Arena* arena) {
  free(arena->base);
  arena->base = 0;
  arena->cap = 0;
  arena->off = 0;
}

static inline uintptr_t ArenaMark(const Arena* arena) {
  return arena->off;
}

static inline void ArenaReset(Arena* arena, uintptr_t mark) {
  arena->off = mark;
}

static inline void* ArenaAlloc(Arena* arena, uintptr_t nbytes) {
  uintptr_t align = 64u;
  uintptr_t off = (arena->off + (align - 1u)) & ~(align - 1u);
  if (nbytes > arena->cap || off > arena->cap - nbytes) {
    return 0;
  }
  void* ptr = arena->base + off;
  arena->off = off + nbytes;
  return ptr;
}

static inline void* ArenaAllocZero(Arena* arena, uintptr_t nbytes) {
  void* ptr = ArenaAlloc(arena, nbytes);
  if (ptr) {
    memset(ptr, 0, nbytes);
  }
  return ptr;
}

static inline char* ArenaStrdup(Arena* arena, const char* src) {
  uintptr_t n = strlen(src) + 1u;
  char* dst = (char*)ArenaAlloc(arena, n);
  if (!dst) {
    return 0;
  }
  memcpy(dst, src, n);
  return dst;
}

static inline uint32_t PtrRoom(const unsigned char* ptr, const unsigned char* end, uint32_t nbytes) {
  if (ptr > end) {
    return 0;
  }
  return (uintptr_t)(end - ptr) >= (uintptr_t)nbytes;
}

#ifdef __cplusplus
}
#endif

#endif
