#ifndef PGEN_READ_H
#define PGEN_READ_H

#include "pgen_ldmeta.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Reader for variable-width phased biallelic hardcalls (storage mode 0x10):
 * header, 2^16-variant block index, vrtype, raw, sparse (1), LD (2/3),
 * simple difflist (4/6/7), and the hardcall phase track. Dosages are not
 * decoded. Genotype codes are the plink2 2-bit codes: 0 hom-ref, 1 het,
 * 2 hom-alt, 3 missing.
 */
typedef struct PgenFile {
  FILE* fp;
  uint32_t variant_ct;
  uint32_t sample_ct;
  unsigned char* vrtypes;
  uint32_t* rec_lens;
  uint64_t* rec_fpos;
  unsigned char* recbuf;
  uint32_t recbuf_cap;
  uintptr_t* anchor_genovec;
  uint32_t anchor_vidx;
  uint32_t anchor_nword;
} PgenFile;

void PgenClose(PgenFile* pg);

PglErr PgenOpen(PgenFile* pg, const char* path);

static inline uint32_t PgenVariantCt(const PgenFile* pg) {
  return pg->variant_ct;
}

static inline uint32_t PgenSampleCt(const PgenFile* pg) {
  return pg->sample_ct;
}

static inline uint32_t PgenVrtypeAt(const PgenFile* pg, uint32_t vidx) {
  return pg->vrtypes[vidx];
}

static inline uint32_t PgenRecLen(const PgenFile* pg, uint32_t vidx) {
  return pg->rec_lens[vidx];
}

/* LD records only. Reads the leading difflist varint and nothing else. */
PglErr PgenLdDiffCount(PgenFile* pg, uint32_t vidx, uint32_t* count);

/*
 * Decode hardcalls and the phase track.
 * phasepresent bit i is set when sample i is a phased het.
 * phaseinfo bit i is set when that het is swapped (1|0). Clear means 0|1.
 * Both bit arrays are sized BitCtToWordCt(sample_ct) and are fully written.
 */
PglErr PgenGetP(PgenFile* pg, uint32_t vidx,
                uintptr_t* genovec, uintptr_t* phasepresent, uintptr_t* phaseinfo,
                uint32_t* phasepresent_ct);

#ifdef __cplusplus
}
#endif

#endif
