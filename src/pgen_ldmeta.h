#ifndef PGEN_LDMETA_H
#define PGEN_LDMETA_H

#include "pgl_base.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Main-track compression is the low 3 bits of vrtype. Phase is bit 4, a
 * separate aux track. Do not use a 0x02/0x04/0x06 table.
 *
 *   0 raw 2-bit
 *   1 one-bit + difflist (sparse)
 *   2 LD difflist against the latest uncompressed anchor
 *   3 inverted LD difflist (patch, then swap hom-ref/hom-alt)
 *   4/6/7 simple difflist, base genotype 0/2/3
 *   5 all homozygous-ref (pgenlib; the written spec reserved this code)
 */

#define PGEN_VR_RAW 0u
#define PGEN_VR_SPARSE 1u
#define PGEN_VR_LD 2u
#define PGEN_VR_LD_INV 3u
#define PGEN_VR_DIFF0 4u
#define PGEN_VR_ALLREF 5u
#define PGEN_VR_DIFF2 6u
#define PGEN_VR_DIFF3 7u
#define PGEN_VR_PHASE 0x10u
#define PGEN_VR_MULTI 0x08u

static inline uint32_t PgenVrtypeMain(uint32_t vrtype) {
  return vrtype & 7u;
}

/* Types 2 and 3: bit 1 set, bit 2 clear. Matches pgenlib VrtypeLdCompressed. */
static inline uint32_t PgenVrtypeLd(uint32_t vrtype) {
  return (vrtype & 6u) == 2u;
}

static inline uint32_t PgenVrtypePhase(uint32_t vrtype) {
  return (vrtype & PGEN_VR_PHASE) != 0u;
}

static inline uint32_t PgenVrtypeMulti(uint32_t vrtype) {
  return (vrtype & PGEN_VR_MULTI) != 0u;
}

/*
 * Difflist length is the leading base-128 varint, not record_bytes/3.
 * Returns 1 on success. *consumed is the varint width in bytes.
 * Values that do not fit in 31 bits are rejected.
 */
static inline uint32_t PgenReadDifflistCount(const unsigned char* rec,
                                             const unsigned char* rec_end,
                                             uint32_t* count,
                                             uint32_t* consumed) {
  const unsigned char* ptr = rec;
  uint32_t vint = 0;
  uint32_t shift = 0;
  if (ptr >= rec_end) {
    return 0;
  }
  while (1) {
    uint32_t byte_val;
    if (ptr >= rec_end) {
      return 0;
    }
    byte_val = *ptr++;
    if (shift == 28u) {
      if (byte_val & 0xf0u) {
        return 0;
      }
      vint |= byte_val << 28;
      break;
    }
    vint |= (byte_val & 127u) << shift;
    if ((byte_val & 128u) == 0u) {
      break;
    }
    shift += 7u;
  }
  if (vint > 0x7fffffffu) {
    return 0;
  }
  *count = vint;
  *consumed = (uint32_t)(ptr - rec);
  return 1;
}

#ifdef __cplusplus
}
#endif

#endif
