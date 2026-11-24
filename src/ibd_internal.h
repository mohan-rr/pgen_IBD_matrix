#ifndef IBD_INTERNAL_H
#define IBD_INTERNAL_H

#include "pgen_read.h"

#ifdef __cplusplus
extern "C" {
#endif

#define kIbdWindowCm 3.0
#define kIbdHalfCm 1.5
/* Merge only overlaps and sub-step gaps so isolated one-site hits do not chain. */
#define kIbdStitchGapCm 0.05
#define kIbdRawMinCm 0.5

enum {
  kIbdRoleDrop = 0,
  kIbdRoleExtend = 1,
  kIbdRoleSeeded = 2,
  kIbdRoleLdSkip = 3
};

typedef struct IbdConfig {
  const char* pgen_path;
  const char* pvar_path;
  const char* psam_path;
  const char* map_path;
  const char* out_prefix;
  double min_length_cm;
  double maf_seed;
  uint32_t mac_extend;
  uint32_t perm_ct;
  uint32_t band_ct;
  int32_t ld_skip; /* negative: default sample_ct/100 */
} IbdConfig;

typedef struct IbdTrace {
  uint8_t* roles; /* variant_ct entries; caller frees */
  uint32_t variant_ct;
  uint32_t sample_ct;
  uint32_t ld_skip_ct;
  uint32_t segment_ct;
  uint32_t pair_ct;
} IbdTrace;

typedef struct SampleInfo {
  char* fid;
  char* iid;
} SampleInfo;

typedef struct ExtSite {
  uint32_t variant_idx;
  uint32_t bp;
  double cm;
  char* chrom;
  char* id;
  uintptr_t* hap_alt; /* bit per haplotype, length BitCtToWordCt(hap_ct) */
} ExtSite;

typedef struct HapPair {
  uint32_t lo;
  uint32_t hi;
} HapPair;

typedef struct RawSeg {
  uint32_t hap_lo;
  uint32_t hap_hi;
  uint32_t var_lo;
  uint32_t var_hi;
  uint32_t bp_lo;
  uint32_t bp_hi;
  double cm_lo;
  double cm_hi;
  double sim;
  char* chrom;
  char* id_lo;
  char* id_hi;
} RawSeg;

void IbdConfigDefaults(IbdConfig* cfg);

const char* PglErrName(PglErr reterr);

/* Keyed 32-bit hash of the variant index. Never returns UINT32_MAX. */
uint32_t IbdPermIndex(uint32_t variant_idx, uint32_t perm_idx);

PglErr IbdExtendWindow(const ExtSite* sites, uint32_t ext_ct,
                       const uintptr_t* hap_bits, uint32_t hap_word_ct, uint32_t hap_ct,
                       const HapPair* pairs, uint32_t pair_ct,
                       double win_cm0, double win_cm1,
                       RawSeg** segs, uint32_t* seg_ct, uint32_t* seg_cap);

PglErr IbdStitchAndWrite(const SampleInfo* samples, uint32_t sample_ct,
                         RawSeg* segs, uint32_t seg_ct,
                         double min_length_cm, const char* out_prefix,
                         uint32_t* out_seg_ct, uint32_t* out_pair_ct);

PglErr IbdRun(const IbdConfig* cfg, IbdTrace* trace);

#ifdef __cplusplus
}
#endif

#endif
