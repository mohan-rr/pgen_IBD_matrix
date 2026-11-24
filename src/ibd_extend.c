#include "ibd_internal.h"

static uint32_t HapBit(const uintptr_t* hap_bits, uint32_t hap_word_ct, uint32_t hap, uint32_t ext_i) {
  return GetBit(hap_bits + (uintptr_t)hap * hap_word_ct, ext_i);
}

static uint32_t SiteMismatch(const uintptr_t* hap_bits, uint32_t hap_word_ct,
                             uint32_t hap_a, uint32_t hap_b, uint32_t ext_i) {
  return HapBit(hap_bits, hap_word_ct, hap_a, ext_i) ^ HapBit(hap_bits, hap_word_ct, hap_b, ext_i);
}

/* XOR + popcount over extend-site bits [lo, hi). */
static uint32_t MismatchRange(const uintptr_t* hap_bits, uint32_t hap_word_ct,
                              uint32_t hap_a, uint32_t hap_b, uint32_t lo, uint32_t hi) {
  const uintptr_t* aa = hap_bits + (uintptr_t)hap_a * hap_word_ct;
  const uintptr_t* bb = hap_bits + (uintptr_t)hap_b * hap_word_ct;
  uint32_t err = 0;
  uint32_t idx = lo;
  if (lo >= hi) {
    return 0;
  }
  while (idx < hi && (idx % kBitsPerWord) != 0u) {
    err += (uint32_t)(((aa[idx / kBitsPerWord] ^ bb[idx / kBitsPerWord]) >> (idx % kBitsPerWord)) & 1u);
    idx++;
  }
  while (idx + kBitsPerWord <= hi) {
    err += PopcountW(aa[idx / kBitsPerWord] ^ bb[idx / kBitsPerWord]);
    idx += kBitsPerWord;
  }
  while (idx < hi) {
    err += (uint32_t)(((aa[idx / kBitsPerWord] ^ bb[idx / kBitsPerWord]) >> (idx % kBitsPerWord)) & 1u);
    idx++;
  }
  return err;
}

/* About one mismatch per 200 extend sites. */
static uint32_t WithinBudget(uint32_t err, uint32_t len) {
  if (len == 0u) {
    return 1u;
  }
  return err <= len / 200u;
}

static void TrimToBudget(const uintptr_t* hap_bits, uint32_t hap_word_ct,
                         uint32_t hap_a, uint32_t hap_b,
                         uint32_t w0, uint32_t w1, uint32_t* out_lo, uint32_t* out_hi) {
  uint32_t lo = w0;
  uint32_t hi = w1;
  uint32_t err;
  uint32_t guard = 0;
  if (w0 >= w1) {
    *out_lo = w0;
    *out_hi = w0;
    return;
  }
  err = MismatchRange(hap_bits, hap_word_ct, hap_a, hap_b, lo, hi);
  while (lo < hi && !WithinBudget(err, hi - lo)) {
    uint32_t prev_lo = lo;
    uint32_t prev_hi = hi;
    if (++guard > (w1 - w0) + 2u) {
      lo = hi;
      break;
    }
    if (SiteMismatch(hap_bits, hap_word_ct, hap_a, hap_b, lo)) {
      err--;
      lo++;
      continue;
    }
    if (SiteMismatch(hap_bits, hap_word_ct, hap_a, hap_b, hi - 1u)) {
      err--;
      hi--;
      continue;
    }
    {
      uint32_t left_end = lo;
      uint32_t right_start = hi;
      uint32_t left_len;
      uint32_t right_len;
      while (left_end < hi && !SiteMismatch(hap_bits, hap_word_ct, hap_a, hap_b, left_end)) {
        left_end++;
      }
      while (right_start > lo && !SiteMismatch(hap_bits, hap_word_ct, hap_a, hap_b, right_start - 1u)) {
        right_start--;
      }
      left_len = left_end - lo;
      right_len = hi - right_start;
      if (left_len >= right_len) {
        hi = left_end;
      } else {
        lo = right_start;
      }
      err = MismatchRange(hap_bits, hap_word_ct, hap_a, hap_b, lo, hi);
    }
    if (lo == prev_lo && hi == prev_hi) {
      lo = hi;
      break;
    }
  }
  *out_lo = lo;
  *out_hi = hi;
}

static void ExtendRight(const uintptr_t* hap_bits, uint32_t hap_word_ct,
                        uint32_t hap_a, uint32_t hap_b, uint32_t ext_ct,
                        uint32_t* lo, uint32_t* hi, uint32_t* err) {
  while (*hi < ext_ct) {
    if ((*hi % kBitsPerWord) == 0u && *hi + kBitsPerWord <= ext_ct) {
      uint32_t add = MismatchRange(hap_bits, hap_word_ct, hap_a, hap_b, *hi, *hi + kBitsPerWord);
      if (WithinBudget(*err + add, *hi + kBitsPerWord - *lo)) {
        *err += add;
        *hi += kBitsPerWord;
        continue;
      }
    }
    {
      uint32_t bit = SiteMismatch(hap_bits, hap_word_ct, hap_a, hap_b, *hi);
      if (!WithinBudget(*err + bit, *hi + 1u - *lo)) {
        break;
      }
      *err += bit;
      (*hi)++;
    }
  }
}

static void ExtendLeft(const uintptr_t* hap_bits, uint32_t hap_word_ct,
                       uint32_t hap_a, uint32_t hap_b,
                       uint32_t* lo, uint32_t* hi, uint32_t* err) {
  while (*lo > 0u) {
    uint32_t bit = SiteMismatch(hap_bits, hap_word_ct, hap_a, hap_b, *lo - 1u);
    if (!WithinBudget(*err + bit, *hi - (*lo - 1u))) {
      break;
    }
    *err += bit;
    (*lo)--;
  }
}

static PglErr SegPush(RawSeg** segs, uint32_t* seg_ct, uint32_t* seg_cap, const RawSeg* seg) {
  if (*seg_ct == *seg_cap) {
    uint32_t cap = *seg_cap ? *seg_cap * 2u : 64u;
    RawSeg* fresh = (RawSeg*)realloc(*segs, (size_t)cap * sizeof(RawSeg));
    if (unlikely(!fresh)) {
      return kPglRetNomem;
    }
    *segs = fresh;
    *seg_cap = cap;
  }
  (*segs)[(*seg_ct)++] = *seg;
  return kPglRetSuccess;
}

PglErr IbdExtendWindow(const ExtSite* sites, uint32_t ext_ct,
                       const uintptr_t* hap_bits, uint32_t hap_word_ct, uint32_t hap_ct,
                       const HapPair* pairs, uint32_t pair_ct,
                       double win_cm0, double win_cm1,
                       RawSeg** segs, uint32_t* seg_ct, uint32_t* seg_cap) {
  uint32_t w0 = 0;
  uint32_t w1 = 0;
  uint32_t pair_idx;
  (void)hap_ct;
  while (w0 < ext_ct && sites[w0].cm < win_cm0) {
    w0++;
  }
  w1 = w0;
  while (w1 < ext_ct && sites[w1].cm < win_cm1) {
    w1++;
  }
  for (pair_idx = 0; pair_idx < pair_ct; pair_idx++) {
    uint32_t hap_a = pairs[pair_idx].lo;
    uint32_t hap_b = pairs[pair_idx].hi;
    uint32_t lo = 0;
    uint32_t hi = 0;
    uint32_t err;
    uint32_t nsite;
    RawSeg seg;
    TrimToBudget(hap_bits, hap_word_ct, hap_a, hap_b, w0, w1, &lo, &hi);
    if (lo >= hi) {
      continue;
    }
    err = MismatchRange(hap_bits, hap_word_ct, hap_a, hap_b, lo, hi);
    ExtendRight(hap_bits, hap_word_ct, hap_a, hap_b, ext_ct, &lo, &hi, &err);
    ExtendLeft(hap_bits, hap_word_ct, hap_a, hap_b, &lo, &hi, &err);
    if (hi <= lo) {
      continue;
    }
    nsite = hi - lo;
    if (sites[hi - 1u].cm - sites[lo].cm < kIbdRawMinCm) {
      continue;
    }
    memset(&seg, 0, sizeof(seg));
    seg.hap_lo = hap_a;
    seg.hap_hi = hap_b;
    seg.var_lo = sites[lo].variant_idx;
    seg.var_hi = sites[hi - 1u].variant_idx;
    seg.bp_lo = sites[lo].bp;
    seg.bp_hi = sites[hi - 1u].bp;
    seg.cm_lo = sites[lo].cm;
    seg.cm_hi = sites[hi - 1u].cm;
    seg.sim = 1.0 - ((double)err / (double)nsite);
    seg.chrom = sites[lo].chrom;
    seg.id_lo = sites[lo].id;
    seg.id_hi = sites[hi - 1u].id;
    if (unlikely(SegPush(segs, seg_ct, seg_cap, &seg))) {
      return kPglRetNomem;
    }
  }
  return kPglRetSuccess;
}

static int CmpRawSeg(const void* va, const void* vb) {
  const RawSeg* sa = (const RawSeg*)va;
  const RawSeg* sb = (const RawSeg*)vb;
  int cmp;
  if (sa->hap_lo != sb->hap_lo) {
    return sa->hap_lo < sb->hap_lo ? -1 : 1;
  }
  if (sa->hap_hi != sb->hap_hi) {
    return sa->hap_hi < sb->hap_hi ? -1 : 1;
  }
  cmp = strcmp(sa->chrom, sb->chrom);
  if (cmp) {
    return cmp;
  }
  if (sa->cm_lo < sb->cm_lo) {
    return -1;
  }
  if (sa->cm_lo > sb->cm_lo) {
    return 1;
  }
  if (sa->cm_hi < sb->cm_hi) {
    return -1;
  }
  if (sa->cm_hi > sb->cm_hi) {
    return 1;
  }
  return 0;
}

static const char* HapName(uint32_t hap) {
  return (hap & 1u) ? "SAMPLE_1" : "SAMPLE_0";
}

static uint32_t StitchSegs(RawSeg* segs, uint32_t seg_ct, double min_length_cm) {
  uint32_t kept = 0;
  uint32_t idx;
  if (seg_ct == 0u) {
    return 0;
  }
  qsort(segs, seg_ct, sizeof(RawSeg), CmpRawSeg);
  for (idx = 0; idx < seg_ct; idx++) {
    RawSeg* cur = &segs[idx];
    if (kept > 0u) {
      RawSeg* prev = &segs[kept - 1u];
      if (prev->hap_lo == cur->hap_lo && prev->hap_hi == cur->hap_hi &&
          strcmp(prev->chrom, cur->chrom) == 0 &&
          cur->cm_lo <= prev->cm_hi + kIbdStitchGapCm) {
        if (cur->cm_hi > prev->cm_hi) {
          prev->cm_hi = cur->cm_hi;
          prev->var_hi = cur->var_hi;
          prev->bp_hi = cur->bp_hi;
          prev->id_hi = cur->id_hi;
        }
        if (cur->cm_lo < prev->cm_lo) {
          prev->cm_lo = cur->cm_lo;
          prev->var_lo = cur->var_lo;
          prev->bp_lo = cur->bp_lo;
          prev->id_lo = cur->id_lo;
        }
        if (cur->sim > prev->sim) {
          prev->sim = cur->sim;
        }
        continue;
      }
    }
    segs[kept++] = *cur;
  }
  {
    uint32_t write = 0;
    for (idx = 0; idx < kept; idx++) {
      if (segs[idx].cm_hi - segs[idx].cm_lo + 1e-12 >= min_length_cm) {
        segs[write++] = segs[idx];
      }
    }
    return write;
  }
}

typedef struct PairSum {
  uint32_t hap_lo;
  uint32_t hap_hi;
  double shared_cm;
  uint32_t nseg;
} PairSum;

PglErr IbdStitchAndWrite(const SampleInfo* samples, uint32_t sample_ct,
                         RawSeg* segs, uint32_t seg_ct,
                         double min_length_cm, const char* out_prefix,
                         uint32_t* out_seg_ct, uint32_t* out_pair_ct) {
  char path[4096];
  FILE* seg_fp = 0;
  FILE* mat_fp = 0;
  PairSum* pairs = 0;
  uint32_t pair_ct = 0;
  uint32_t pair_cap = 0;
  uint32_t idx;
  PglErr reterr = kPglRetSuccess;
  int nprint;
  (void)sample_ct;
  seg_ct = StitchSegs(segs, seg_ct, min_length_cm);
  nprint = snprintf(path, sizeof(path), "%s.segments", out_prefix);
  if (nprint < 0 || (size_t)nprint >= sizeof(path)) {
    return kPglRetInvalidCmdline;
  }
  seg_fp = fopen(path, "wb");
  if (unlikely(!seg_fp)) {
    reterr = kPglRetOpenFail;
    goto cleanup;
  }
  if (fprintf(seg_fp, "#FID1 IID1 HAP1 FID2 IID2 HAP2 CHR START_BP END_BP START_ID END_ID CM SIMILARITY\n") < 0) {
    reterr = kPglRetWriteFail;
    goto cleanup;
  }
  for (idx = 0; idx < seg_ct; idx++) {
    const RawSeg* seg = &segs[idx];
    uint32_t s0 = seg->hap_lo / 2u;
    uint32_t s1 = seg->hap_hi / 2u;
    if (fprintf(seg_fp, "%s %s %s %s %s %s %s %u %u %s %s %.6f %.6f\n",
                samples[s0].fid, samples[s0].iid, HapName(seg->hap_lo),
                samples[s1].fid, samples[s1].iid, HapName(seg->hap_hi),
                seg->chrom, seg->bp_lo, seg->bp_hi, seg->id_lo, seg->id_hi,
                seg->cm_hi - seg->cm_lo, seg->sim) < 0) {
      reterr = kPglRetWriteFail;
      goto cleanup;
    }
    if (pair_ct == 0u || pairs[pair_ct - 1u].hap_lo != seg->hap_lo || pairs[pair_ct - 1u].hap_hi != seg->hap_hi) {
      if (pair_ct == pair_cap) {
        uint32_t cap = pair_cap ? pair_cap * 2u : 32u;
        PairSum* fresh = (PairSum*)realloc(pairs, (size_t)cap * sizeof(PairSum));
        if (unlikely(!fresh)) {
          reterr = kPglRetNomem;
          goto cleanup;
        }
        pairs = fresh;
        pair_cap = cap;
      }
      pairs[pair_ct].hap_lo = seg->hap_lo;
      pairs[pair_ct].hap_hi = seg->hap_hi;
      pairs[pair_ct].shared_cm = 0;
      pairs[pair_ct].nseg = 0;
      pair_ct++;
    }
    pairs[pair_ct - 1u].shared_cm += seg->cm_hi - seg->cm_lo;
    pairs[pair_ct - 1u].nseg++;
  }
  nprint = snprintf(path, sizeof(path), "%s.ibd", out_prefix);
  if (nprint < 0 || (size_t)nprint >= sizeof(path)) {
    reterr = kPglRetInvalidCmdline;
    goto cleanup;
  }
  mat_fp = fopen(path, "wb");
  if (unlikely(!mat_fp)) {
    reterr = kPglRetOpenFail;
    goto cleanup;
  }
  if (fprintf(mat_fp, "#SAMPLE1 HAP1 SAMPLE2 HAP2 SHARED_CM N_SEGMENTS\n") < 0) {
    reterr = kPglRetWriteFail;
    goto cleanup;
  }
  for (idx = 0; idx < pair_ct; idx++) {
    uint32_t s0 = pairs[idx].hap_lo / 2u;
    uint32_t s1 = pairs[idx].hap_hi / 2u;
    if (fprintf(mat_fp, "%s %s %s %s %.6f %u\n",
                samples[s0].iid, HapName(pairs[idx].hap_lo),
                samples[s1].iid, HapName(pairs[idx].hap_hi),
                pairs[idx].shared_cm, pairs[idx].nseg) < 0) {
      reterr = kPglRetWriteFail;
      goto cleanup;
    }
  }
  if (out_seg_ct) {
    *out_seg_ct = seg_ct;
  }
  if (out_pair_ct) {
    *out_pair_ct = pair_ct;
  }
  cleanup:
  if (seg_fp) {
    if (fclose(seg_fp) != 0 && !reterr) {
      reterr = kPglRetWriteFail;
    }
  }
  if (mat_fp) {
    if (fclose(mat_fp) != 0 && !reterr) {
      reterr = kPglRetWriteFail;
    }
  }
  free(pairs);
  return reterr;
}
