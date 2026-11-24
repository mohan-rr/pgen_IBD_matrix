#include "ibd_internal.h"

void IbdConfigDefaults(IbdConfig* cfg) {
  memset(cfg, 0, sizeof(*cfg));
  cfg->min_length_cm = 2.0;
  cfg->maf_seed = 0.01;
  cfg->mac_extend = 2u;
  cfg->perm_ct = 40u;
  cfg->band_ct = 10u;
  cfg->ld_skip = -1;
}

const char* PglErrName(PglErr reterr) {
  switch (reterr) {
    case kPglRetSuccess: return "success";
    case kPglRetNomem: return "out of memory";
    case kPglRetOpenFail: return "open failed";
    case kPglRetReadFail: return "read failed";
    case kPglRetWriteFail: return "write failed";
    case kPglRetMalformedInput: return "malformed input";
    case kPglRetInconsistentInput: return "inconsistent input";
    case kPglRetInvalidCmdline: return "invalid command line";
    case kPglRetNotYetSupported: return "not yet supported";
    case kPglRetVarRecordTooLarge: return "variant record too large";
    case kPglRetLongLine: return "line too long";
    case kPglRetDegenerateData: return "degenerate data";
    default: return "error";
  }
}

uint32_t IbdPermIndex(uint32_t variant_idx, uint32_t perm_idx) {
  uint32_t mix = (variant_idx + 1u) * 0x9E3779B9u;
  mix ^= (perm_idx + 1u) * 0x85EBCA6Bu;
  mix ^= mix >> 16;
  mix *= 0x7FEB352Du;
  mix ^= mix >> 15;
  mix *= 0x846CA68Bu;
  mix ^= mix >> 16;
  if (mix == 0xffffffffu) {
    mix = 0;
  }
  return mix;
}

static uint32_t HalfOf(double cm, double origin) {
  double x = (cm - origin) / kIbdHalfCm;
  double snapped;
  if (x < 0.0) {
    x = 0.0;
  }
  snapped = (double)(uint32_t)(x + 0.5);
  if (x < snapped && snapped - x < 1e-8) {
    x = snapped;
  }
  return (uint32_t)x;
}

static int ReadLine(FILE* fp, char* buf, uint32_t cap, PglErr* reterr) {
  size_t len;
  if (!fgets(buf, (int)cap, fp)) {
    if (feof(fp)) {
      return 0;
    }
    *reterr = kPglRetReadFail;
    return -1;
  }
  len = strlen(buf);
  if (len && buf[len - 1u] == '\n') {
    buf[--len] = 0;
    if (len && buf[len - 1u] == '\r') {
      buf[--len] = 0;
    }
    return 1;
  }
  if (!feof(fp)) {
    *reterr = kPglRetLongLine;
    return -1;
  }
  return 1;
}

static uint32_t SplitWs(char* line, char** tok, uint32_t max_tok) {
  uint32_t ntok = 0;
  char* ptr = line;
  while (*ptr && ntok < max_tok) {
    while (*ptr == ' ' || *ptr == '\t' || *ptr == '\r') {
      ptr++;
    }
    if (*ptr == 0) {
      break;
    }
    tok[ntok++] = ptr;
    while (*ptr && *ptr != ' ' && *ptr != '\t' && *ptr != '\r') {
      ptr++;
    }
    if (*ptr) {
      *ptr++ = 0;
    }
  }
  return ntok;
}

static uint32_t ParseU32(const char* text, uint32_t* out) {
  char* end = 0;
  unsigned long value;
  if (!text || !text[0]) {
    return 0;
  }
  value = strtoul(text, &end, 10);
  if (end == text || *end != 0 || value > 0xfffffffful) {
    return 0;
  }
  *out = (uint32_t)value;
  return 1;
}

static uint32_t ParseF64(const char* text, double* out) {
  char* end = 0;
  double value;
  if (!text || !text[0]) {
    return 0;
  }
  value = strtod(text, &end);
  if (end == text || *end != 0) {
    return 0;
  }
  *out = value;
  return 1;
}

static uint32_t IsAcgtBase(char base) {
  if (base >= 'a' && base <= 'z') {
    base = (char)(base - 32);
  }
  return base == 'A' || base == 'C' || base == 'G' || base == 'T';
}

static uint32_t IsBiallelicSnp(const char* ref, const char* alt) {
  if (!ref || !alt || !ref[0] || ref[1] || !alt[0] || alt[1]) {
    return 0;
  }
  if (strchr(alt, ',') || alt[0] == '<' || ref[0] == '<') {
    return 0;
  }
  return IsAcgtBase(ref[0]) && IsAcgtBase(alt[0]);
}

typedef struct VarInfo {
  char* chrom;
  char* id;
  uint32_t bp;
  double cm;
  uint8_t vrtype;
  uint8_t is_snp;
  uint8_t role;
  uintptr_t* hap_alt;
} VarInfo;

typedef struct MapPoint {
  char* chrom;
  uint32_t bp;
  double cm;
} MapPoint;

typedef struct ChromMap {
  char* chrom;
  uint32_t start;
  uint32_t count;
} ChromMap;

static int CmpMapPoint(const void* va, const void* vb) {
  const MapPoint* pa = (const MapPoint*)va;
  const MapPoint* pb = (const MapPoint*)vb;
  int cmp = strcmp(pa->chrom, pb->chrom);
  if (cmp) {
    return cmp;
  }
  if (pa->bp < pb->bp) {
    return -1;
  }
  if (pa->bp > pb->bp) {
    return 1;
  }
  return 0;
}

static double MapCm(const MapPoint* pts, uint32_t count, uint32_t bp) {
  uint32_t lo;
  uint32_t hi;
  double frac;
  if (count == 0u) {
    return 0.0;
  }
  if (bp <= pts[0].bp) {
    return pts[0].cm;
  }
  if (bp >= pts[count - 1u].bp) {
    return pts[count - 1u].cm;
  }
  lo = 0;
  hi = count - 1u;
  while (lo + 1u < hi) {
    uint32_t mid = lo + (hi - lo) / 2u;
    if (pts[mid].bp <= bp) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  if (pts[hi].bp == pts[lo].bp) {
    return pts[lo].cm;
  }
  frac = (double)(bp - pts[lo].bp) / (double)(pts[hi].bp - pts[lo].bp);
  return pts[lo].cm + frac * (pts[hi].cm - pts[lo].cm);
}

static const ChromMap* FindChrom(const ChromMap* maps, uint32_t map_ct, const char* chrom) {
  uint32_t idx;
  for (idx = 0; idx < map_ct; idx++) {
    if (strcmp(maps[idx].chrom, chrom) == 0) {
      return &maps[idx];
    }
  }
  return 0;
}

static PglErr LoadMap(Arena* arena, const char* path, MapPoint** out_pts, uint32_t* out_ct,
                      ChromMap** out_maps, uint32_t* out_map_ct) {
  FILE* fp = fopen(path, "rb");
  char line[8192];
  MapPoint* pts = 0;
  uint32_t count = 0;
  uint32_t cap = 0;
  ChromMap* maps = 0;
  uint32_t map_ct = 0;
  PglErr reterr = kPglRetSuccess;
  int got;
  if (unlikely(!fp)) {
    fprintf(stderr, "pgen_ibd: cannot open map %s\n", path);
    return kPglRetOpenFail;
  }
  while ((got = ReadLine(fp, line, sizeof(line), &reterr)) != 0) {
    char* tok[8];
    uint32_t ntok;
    uint32_t bp = 0;
    double cm = 0.0;
    MapPoint point;
    if (got < 0) {
      goto cleanup;
    }
    if (line[0] == 0 || line[0] == '#') {
      continue;
    }
    ntok = SplitWs(line, tok, 8);
    if (ntok < 4u || !ParseF64(tok[2], &cm) || !ParseU32(tok[3], &bp)) {
      continue;
    }
    if (count == cap) {
      uint32_t ncap = cap ? cap * 2u : 256u;
      MapPoint* fresh = (MapPoint*)realloc(pts, (size_t)ncap * sizeof(MapPoint));
      if (unlikely(!fresh)) {
        reterr = kPglRetNomem;
        goto cleanup;
      }
      pts = fresh;
      cap = ncap;
    }
    point.chrom = ArenaStrdup(arena, tok[0]);
    point.bp = bp;
    point.cm = cm;
    if (unlikely(!point.chrom)) {
      reterr = kPglRetNomem;
      goto cleanup;
    }
    pts[count++] = point;
  }
  if (count == 0u) {
    fprintf(stderr, "pgen_ibd: map %s has no PLINK chrom/id/cm/bp rows\n", path);
    reterr = kPglRetMalformedInput;
    goto cleanup;
  }
  qsort(pts, count, sizeof(MapPoint), CmpMapPoint);
  maps = (ChromMap*)malloc(sizeof(ChromMap));
  if (unlikely(!maps)) {
    reterr = kPglRetNomem;
    goto cleanup;
  }
  maps[0].chrom = pts[0].chrom;
  maps[0].start = 0;
  maps[0].count = 1;
  map_ct = 1;
  {
    uint32_t idx;
    for (idx = 1; idx < count; idx++) {
      if (strcmp(pts[idx].chrom, maps[map_ct - 1u].chrom) == 0) {
        maps[map_ct - 1u].count++;
      } else {
        ChromMap* fresh = (ChromMap*)realloc(maps, (size_t)(map_ct + 1u) * sizeof(ChromMap));
        if (unlikely(!fresh)) {
          reterr = kPglRetNomem;
          goto cleanup;
        }
        maps = fresh;
        maps[map_ct].chrom = pts[idx].chrom;
        maps[map_ct].start = idx;
        maps[map_ct].count = 1;
        map_ct++;
      }
    }
  }
  *out_pts = pts;
  *out_ct = count;
  *out_maps = maps;
  *out_map_ct = map_ct;
  pts = 0;
  maps = 0;
  cleanup:
  if (fp) {
    fclose(fp);
  }
  free(pts);
  free(maps);
  return reterr;
}

static PglErr LoadPsam(Arena* arena, const char* path, uint32_t expect_ct,
                       SampleInfo** out_samples) {
  FILE* fp = fopen(path, "rb");
  char line[8192];
  SampleInfo* samples;
  uint32_t count = 0;
  int fid_col = 0;
  int iid_col = 1;
  int have_header = 0;
  PglErr reterr = kPglRetSuccess;
  int got;
  if (unlikely(!fp)) {
    fprintf(stderr, "pgen_ibd: cannot open psam %s\n", path);
    return kPglRetOpenFail;
  }
  samples = (SampleInfo*)ArenaAllocZero(arena, (uintptr_t)expect_ct * sizeof(SampleInfo));
  if (unlikely(!samples)) {
    reterr = kPglRetNomem;
    goto cleanup;
  }
  while ((got = ReadLine(fp, line, sizeof(line), &reterr)) != 0) {
    char* tok[16];
    uint32_t ntok;
    char* fid;
    char* iid;
    if (got < 0) {
      goto cleanup;
    }
    if (line[0] == 0) {
      continue;
    }
    if (line[0] == '#') {
      if (!have_header && (strncmp(line, "#FID", 4) == 0 || strncmp(line, "#IID", 4) == 0)) {
        uint32_t col;
        ntok = SplitWs(line, tok, 16);
        fid_col = -1;
        iid_col = -1;
        for (col = 0; col < ntok; col++) {
          if (strcmp(tok[col], "#FID") == 0 || strcmp(tok[col], "FID") == 0) {
            fid_col = (int)col;
          } else if (strcmp(tok[col], "#IID") == 0 || strcmp(tok[col], "IID") == 0) {
            iid_col = (int)col;
          }
        }
        if (iid_col < 0) {
          reterr = kPglRetMalformedInput;
          goto cleanup;
        }
        have_header = 1;
      }
      continue;
    }
    ntok = SplitWs(line, tok, 16);
    if (fid_col >= 0 && (uint32_t)fid_col >= ntok) {
      reterr = kPglRetMalformedInput;
      goto cleanup;
    }
    if ((uint32_t)iid_col >= ntok) {
      reterr = kPglRetMalformedInput;
      goto cleanup;
    }
    fid = fid_col >= 0 ? tok[fid_col] : (char*)"0";
    iid = tok[iid_col];
    if (count >= expect_ct) {
      fprintf(stderr, "pgen_ibd: psam has more samples than the pgen\n");
      reterr = kPglRetInconsistentInput;
      goto cleanup;
    }
    samples[count].fid = ArenaStrdup(arena, fid);
    samples[count].iid = ArenaStrdup(arena, iid);
    if (unlikely(!samples[count].fid || !samples[count].iid)) {
      reterr = kPglRetNomem;
      goto cleanup;
    }
    count++;
  }
  if (count != expect_ct) {
    fprintf(stderr, "pgen_ibd: psam sample count %u != pgen %u\n", count, expect_ct);
    reterr = kPglRetInconsistentInput;
    goto cleanup;
  }
  *out_samples = samples;
  cleanup:
  if (fp) {
    fclose(fp);
  }
  return reterr;
}

static int ColOf(char** tok, uint32_t ntok, const char* name) {
  uint32_t col;
  for (col = 0; col < ntok; col++) {
    if (strcmp(tok[col], name) == 0) {
      return (int)col;
    }
  }
  return -1;
}

static PglErr LoadPvar(Arena* arena, const char* path, const PgenFile* pg,
                       const MapPoint* pts, const ChromMap* maps, uint32_t map_ct,
                       VarInfo** out_vars) {
  FILE* fp = fopen(path, "rb");
  char line[8192];
  VarInfo* vars;
  uint32_t count = 0;
  int chrom_col = 0;
  int pos_col = 1;
  int id_col = 2;
  int ref_col = 3;
  int alt_col = 4;
  int have_header = 0;
  PglErr reterr = kPglRetSuccess;
  int got;
  if (unlikely(!fp)) {
    fprintf(stderr, "pgen_ibd: cannot open pvar %s\n", path);
    return kPglRetOpenFail;
  }
  vars = (VarInfo*)ArenaAllocZero(arena, (uintptr_t)pg->variant_ct * sizeof(VarInfo));
  if (unlikely(!vars)) {
    reterr = kPglRetNomem;
    goto cleanup;
  }
  while ((got = ReadLine(fp, line, sizeof(line), &reterr)) != 0) {
    char* tok[32];
    uint32_t ntok;
    uint32_t bp = 0;
    const ChromMap* cmap;
    if (got < 0) {
      goto cleanup;
    }
    if (line[0] == 0) {
      continue;
    }
    if (line[0] == '#') {
      if (!have_header && strncmp(line, "#CHROM", 6) == 0) {
        uint32_t col;
        ntok = SplitWs(line, tok, 32);
        for (col = 0; col < ntok; col++) {
          if (strcmp(tok[col], "FORMAT") == 0) {
            ntok = col;
            break;
          }
        }
        chrom_col = ColOf(tok, ntok, "#CHROM");
        pos_col = ColOf(tok, ntok, "POS");
        id_col = ColOf(tok, ntok, "ID");
        ref_col = ColOf(tok, ntok, "REF");
        alt_col = ColOf(tok, ntok, "ALT");
        if (chrom_col < 0 || pos_col < 0 || id_col < 0 || ref_col < 0 || alt_col < 0) {
          fprintf(stderr, "pgen_ibd: pvar header needs CHROM POS ID REF ALT\n");
          reterr = kPglRetMalformedInput;
          goto cleanup;
        }
        have_header = 1;
      }
      continue;
    }
    ntok = SplitWs(line, tok, 32);
    if (!have_header) {
      if (ntok >= 6u) {
        chrom_col = 0;
        id_col = 1;
        pos_col = 3;
        alt_col = 4;
        ref_col = 5;
      } else if (ntok == 5u) {
        chrom_col = 0;
        id_col = 1;
        pos_col = 2;
        alt_col = 3;
        ref_col = 4;
      } else {
        reterr = kPglRetMalformedInput;
        goto cleanup;
      }
      have_header = 1;
    }
    if ((uint32_t)chrom_col >= ntok || (uint32_t)pos_col >= ntok || (uint32_t)id_col >= ntok ||
        (uint32_t)ref_col >= ntok || (uint32_t)alt_col >= ntok || !ParseU32(tok[pos_col], &bp)) {
      reterr = kPglRetMalformedInput;
      goto cleanup;
    }
    if (count >= pg->variant_ct) {
      fprintf(stderr, "pgen_ibd: pvar has more variants than the pgen\n");
      reterr = kPglRetInconsistentInput;
      goto cleanup;
    }
    vars[count].chrom = ArenaStrdup(arena, tok[chrom_col]);
    vars[count].id = ArenaStrdup(arena, tok[id_col]);
    if (unlikely(!vars[count].chrom || !vars[count].id)) {
      reterr = kPglRetNomem;
      goto cleanup;
    }
    vars[count].bp = bp;
    vars[count].vrtype = pg->vrtypes[count];
    vars[count].is_snp = (uint8_t)IsBiallelicSnp(tok[ref_col], tok[alt_col]);
    cmap = FindChrom(maps, map_ct, vars[count].chrom);
    if (!cmap) {
      fprintf(stderr, "pgen_ibd: chromosome %s is missing from the map\n", vars[count].chrom);
      reterr = kPglRetInconsistentInput;
      goto cleanup;
    }
    vars[count].cm = MapCm(pts + cmap->start, cmap->count, bp);
    count++;
  }
  if (count != pg->variant_ct) {
    fprintf(stderr, "pgen_ibd: pvar variant count %u != pgen %u\n", count, pg->variant_ct);
    reterr = kPglRetInconsistentInput;
    goto cleanup;
  }
  *out_vars = vars;
  cleanup:
  if (fp) {
    fclose(fp);
  }
  return reterr;
}

static PglErr BuildMasks(PgenFile* pg, VarInfo* vars, uint32_t mac_extend, double maf_seed,
                         uintptr_t* hap_slab, uint32_t hap_word, uintptr_t* genovec,
                         uintptr_t* phasepresent, uintptr_t* phaseinfo) {
  uint32_t vidx;
  uint32_t sample_ct = pg->sample_ct;
  uint32_t hap_ct = sample_ct * 2u;
  for (vidx = 0; vidx < pg->variant_ct; vidx++) {
    uint32_t phase_ct = 0;
    uint32_t sample_idx;
    uint32_t ref_ac = 0;
    uint32_t alt_ac = 0;
    uint32_t mac;
    uint32_t denom;
    uint32_t major_alt;
    double maf;
    PglErr reterr;
    vars[vidx].hap_alt = hap_slab + (uintptr_t)vidx * hap_word;
    vars[vidx].role = kIbdRoleDrop;
    if (!vars[vidx].is_snp || PgenVrtypeMulti(vars[vidx].vrtype)) {
      continue;
    }
    reterr = PgenGetP(pg, vidx, genovec, phasepresent, phaseinfo, &phase_ct);
    if (unlikely(reterr)) {
      fprintf(stderr, "pgen_ibd: failed to decode %s (%s)\n", vars[vidx].id, PglErrName(reterr));
      return reterr;
    }
    for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
      uint32_t geno = GetGeno(genovec, sample_idx);
      if (geno == 0u) {
        ref_ac += 2u;
      } else if (geno == 1u) {
        ref_ac++;
        alt_ac++;
      } else if (geno == 2u) {
        alt_ac += 2u;
      }
    }
    denom = ref_ac + alt_ac;
    if (denom == 0u) {
      continue;
    }
    mac = ref_ac < alt_ac ? ref_ac : alt_ac;
    if (mac < mac_extend) {
      continue;
    }
    maf = (double)mac / (double)denom;
    major_alt = alt_ac > ref_ac;
    for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
      uint32_t geno = GetGeno(genovec, sample_idx);
      uint32_t h0;
      uint32_t h1;
      if (geno == 3u) {
        h0 = major_alt;
        h1 = major_alt;
      } else if (geno == 0u) {
        h0 = 0;
        h1 = 0;
      } else if (geno == 2u) {
        h0 = 1;
        h1 = 1;
      } else {
        if (!GetBit(phasepresent, sample_idx)) {
          fprintf(stderr, "pgen_ibd: unphased het at %s sample %u\n", vars[vidx].id, sample_idx);
          return kPglRetInconsistentInput;
        }
        if (GetBit(phaseinfo, sample_idx)) {
          h0 = 1;
          h1 = 0;
        } else {
          h0 = 0;
          h1 = 1;
        }
      }
      if (h0) {
        SetBit(vars[vidx].hap_alt, sample_idx * 2u);
      }
      if (h1) {
        SetBit(vars[vidx].hap_alt, sample_idx * 2u + 1u);
      }
    }
    (void)hap_ct;
    vars[vidx].role = (maf >= maf_seed) ? (uint8_t)kIbdRoleSeeded : (uint8_t)kIbdRoleExtend;
  }
  return kPglRetSuccess;
}

typedef struct LshKey {
  uint32_t h0;
  uint32_t h1;
  uint32_t h2;
  uint32_t h3;
  uint32_t hap;
} LshKey;

static uint32_t KeyByte(const LshKey* key, uint32_t pass) {
  uint32_t word_i = 3u - pass / 4u;
  uint32_t byte_i = pass % 4u;
  uint32_t word = key->h3;
  if (word_i == 0u) {
    word = key->h0;
  } else if (word_i == 1u) {
    word = key->h1;
  } else if (word_i == 2u) {
    word = key->h2;
  }
  return (word >> (8u * byte_i)) & 255u;
}

static void RadixSortLsh(LshKey* keys, LshKey* tmp, uint32_t nkey) {
  uint32_t pass;
  LshKey* src = keys;
  LshKey* dst = tmp;
  if (nkey < 2u) {
    return;
  }
  for (pass = 0; pass < 16u; pass++) {
    uint32_t count[256];
    uint32_t sum = 0;
    uint32_t bucket;
    uint32_t idx;
    memset(count, 0, sizeof(count));
    for (idx = 0; idx < nkey; idx++) {
      count[KeyByte(&src[idx], pass)]++;
    }
    for (bucket = 0; bucket < 256u; bucket++) {
      uint32_t c = count[bucket];
      count[bucket] = sum;
      sum += c;
    }
    for (idx = 0; idx < nkey; idx++) {
      uint32_t b = KeyByte(&src[idx], pass);
      dst[count[b]++] = src[idx];
    }
    {
      LshKey* swap = src;
      src = dst;
      dst = swap;
    }
  }
}

static int CmpPair(const void* va, const void* vb) {
  const HapPair* pa = (const HapPair*)va;
  const HapPair* pb = (const HapPair*)vb;
  if (pa->lo != pb->lo) {
    return pa->lo < pb->lo ? -1 : 1;
  }
  if (pa->hi != pb->hi) {
    return pa->hi < pb->hi ? -1 : 1;
  }
  return 0;
}

static uint32_t SigHash(const uint32_t* buf_a, const uint32_t* buf_b, uint32_t hap, uint32_t perm, uint32_t perm_ct) {
  uint32_t left = buf_a[(uintptr_t)hap * perm_ct + perm];
  uint32_t right = buf_b[(uintptr_t)hap * perm_ct + perm];
  return left < right ? left : right;
}

static void UpdateMinhash(uint32_t* partial, uint32_t hap_ct, uint32_t perm_ct,
                          uint32_t variant_idx, const uintptr_t* hap_alt) {
  uint32_t hap;
  for (hap = 0; hap < hap_ct; hap++) {
    uint32_t perm;
    uint32_t* row;
    if (!GetBit(hap_alt, hap)) {
      continue;
    }
    row = partial + (uintptr_t)hap * perm_ct;
    for (perm = 0; perm < perm_ct; perm++) {
      uint32_t pix = IbdPermIndex(variant_idx, perm);
      if (pix < row[perm]) {
        row[perm] = pix;
      }
    }
  }
}

static void FillU32(uint32_t* dst, uint32_t n, uint32_t value) {
  uint32_t idx;
  for (idx = 0; idx < n; idx++) {
    dst[idx] = value;
  }
}

static PglErr EmitWindow(Arena* arena, const uint32_t* buf_a, const uint32_t* buf_b,
                         uint32_t hap_ct, uint32_t perm_ct, uint32_t band_ct,
                         const ExtSite* sites, uint32_t ext_ct,
                         const uintptr_t* hap_bits, uint32_t hap_word_ct,
                         double win_cm0, double win_cm1,
                         RawSeg** segs, uint32_t* seg_ct, uint32_t* seg_cap) {
  uint32_t hashes_per_band = perm_ct / band_ct;
  uint32_t band;
  uintptr_t mark = ArenaMark(arena);
  PglErr reterr = kPglRetSuccess;
  if (hashes_per_band == 0u || hashes_per_band > 4u) {
    ArenaReset(arena, mark);
    return kPglRetInvalidCmdline;
  }
  HapPair* pairs = 0;
  uint32_t pair_ct = 0;
  uint32_t pair_cap = 0;
  for (band = 0; band < band_ct; band++) {
    LshKey* keys;
    LshKey* tmp;
    uint32_t nkey = 0;
    uint32_t hap;
    uint32_t cursor;
    uint32_t base = band * hashes_per_band;
    keys = (LshKey*)ArenaAlloc(arena, (uintptr_t)hap_ct * sizeof(LshKey));
    tmp = (LshKey*)ArenaAlloc(arena, (uintptr_t)hap_ct * sizeof(LshKey));
    if (unlikely(!keys || !tmp)) {
      reterr = kPglRetNomem;
      break;
    }
    for (hap = 0; hap < hap_ct; hap++) {
      uint32_t hashes[4];
      uint32_t hash_i;
      uint32_t all_max = 1;
      LshKey key;
      for (hash_i = 0; hash_i < hashes_per_band; hash_i++) {
        hashes[hash_i] = SigHash(buf_a, buf_b, hap, base + hash_i, perm_ct);
        if (hashes[hash_i] != UINT32_MAX) {
          all_max = 0;
        }
      }
      for (; hash_i < 4u; hash_i++) {
        hashes[hash_i] = 0;
      }
      if (all_max) {
        continue;
      }
      key.h0 = hashes[0];
      key.h1 = hashes_per_band > 1u ? hashes[1] : 0u;
      key.h2 = hashes_per_band > 2u ? hashes[2] : 0u;
      key.h3 = hashes_per_band > 3u ? hashes[3] : 0u;
      key.hap = hap;
      keys[nkey++] = key;
    }
    RadixSortLsh(keys, tmp, nkey);
    cursor = 0;
    while (cursor < nkey) {
      uint32_t run = cursor + 1u;
      while (run < nkey && keys[run].h0 == keys[cursor].h0 && keys[run].h1 == keys[cursor].h1 &&
             keys[run].h2 == keys[cursor].h2 && keys[run].h3 == keys[cursor].h3) {
        run++;
      }
      if (run - cursor >= 2u) {
        uint32_t left;
        for (left = cursor; left + 1u < run; left++) {
          uint32_t right;
          for (right = left + 1u; right < run; right++) {
            HapPair pair;
            pair.lo = keys[left].hap;
            pair.hi = keys[right].hap;
            if (pair.lo > pair.hi) {
              uint32_t swap = pair.lo;
              pair.lo = pair.hi;
              pair.hi = swap;
            }
            if (pair.lo == pair.hi) {
              continue;
            }
            if (pair_ct == pair_cap) {
              uint32_t ncap = pair_cap ? pair_cap * 2u : 32u;
              HapPair* fresh = (HapPair*)realloc(pairs, (size_t)ncap * sizeof(HapPair));
              if (unlikely(!fresh)) {
                reterr = kPglRetNomem;
                break;
              }
              pairs = fresh;
              pair_cap = ncap;
            }
            pairs[pair_ct++] = pair;
          }
        }
      }
      cursor = run;
    }
    if (reterr) {
      break;
    }
  }
  if (!reterr && pair_ct > 1u) {
    uint32_t write = 0;
    uint32_t idx;
    qsort(pairs, pair_ct, sizeof(HapPair), CmpPair);
    for (idx = 0; idx < pair_ct; idx++) {
      if (write == 0u || pairs[idx].lo != pairs[write - 1u].lo || pairs[idx].hi != pairs[write - 1u].hi) {
        pairs[write++] = pairs[idx];
      }
    }
    pair_ct = write;
  }
  if (!reterr && pair_ct) {
    reterr = IbdExtendWindow(sites, ext_ct, hap_bits, hap_word_ct, hap_ct, pairs, pair_ct,
                             win_cm0, win_cm1, segs, seg_ct, seg_cap);
  }
  free(pairs);
  ArenaReset(arena, mark);
  return reterr;
}

static PglErr RunChrom(PgenFile* pg, Arena* arena, VarInfo* vars, uint32_t c0, uint32_t c1,
                       const IbdConfig* cfg, uint32_t ld_lim, uint32_t* ld_skip_ct,
                       RawSeg** segs, uint32_t* seg_ct, uint32_t* seg_cap) {
  uint32_t sample_ct = pg->sample_ct;
  uint32_t hap_ct = sample_ct * 2u;
  uint32_t perm_ct = cfg->perm_ct;
  uint32_t vidx;
  uint32_t ext_ct = 0;
  uint32_t ext_i;
  uint32_t hap_word_ct;
  uint32_t partial_n;
  uint32_t max_half;
  uint32_t half;
  uint32_t anchor_seed_half = UINT32_MAX;
  double origin;
  ExtSite* sites;
  uintptr_t* hap_bits;
  uint32_t* buf[2];
  uintptr_t mark;
  PglErr reterr = kPglRetSuccess;
  if (c0 >= c1) {
    return kPglRetSuccess;
  }
  for (vidx = c0 + 1u; vidx < c1; vidx++) {
    if (vars[vidx].cm + 1e-9 < vars[vidx - 1u].cm) {
      fprintf(stderr, "pgen_ibd: genetic positions decrease at %s\n", vars[vidx].id);
      return kPglRetInconsistentInput;
    }
  }
  for (vidx = c0; vidx < c1; vidx++) {
    if (vars[vidx].role != kIbdRoleDrop) {
      ext_ct++;
    }
  }
  mark = ArenaMark(arena);
  hap_word_ct = ext_ct ? BitCtToWordCt(ext_ct) : 1u;
  sites = (ExtSite*)ArenaAllocZero(arena, ((uintptr_t)ext_ct + 1u) * sizeof(ExtSite));
  hap_bits = (uintptr_t*)ArenaAllocZero(arena, (uintptr_t)hap_ct * hap_word_ct * sizeof(uintptr_t));
  partial_n = hap_ct * perm_ct;
  buf[0] = (uint32_t*)ArenaAlloc(arena, (uintptr_t)partial_n * sizeof(uint32_t));
  buf[1] = (uint32_t*)ArenaAlloc(arena, (uintptr_t)partial_n * sizeof(uint32_t));
  if (unlikely(!sites || !hap_bits || !buf[0] || !buf[1])) {
    ArenaReset(arena, mark);
    return kPglRetNomem;
  }
  ext_i = 0;
  for (vidx = c0; vidx < c1; vidx++) {
    uint32_t hap;
    if (vars[vidx].role == kIbdRoleDrop) {
      continue;
    }
    sites[ext_i].variant_idx = vidx;
    sites[ext_i].bp = vars[vidx].bp;
    sites[ext_i].cm = vars[vidx].cm;
    sites[ext_i].chrom = vars[vidx].chrom;
    sites[ext_i].id = vars[vidx].id;
    sites[ext_i].hap_alt = vars[vidx].hap_alt;
    for (hap = 0; hap < hap_ct; hap++) {
      if (GetBit(vars[vidx].hap_alt, hap)) {
        SetBit(hap_bits + (uintptr_t)hap * hap_word_ct, ext_i);
      }
    }
    ext_i++;
  }
  origin = vars[c0].cm;
  max_half = HalfOf(vars[c1 - 1u].cm, origin);
  FillU32(buf[0], partial_n, UINT32_MAX);
  FillU32(buf[1], partial_n, UINT32_MAX);
  vidx = c0;
  for (half = 0; half <= max_half; half++) {
    uint32_t* dest = buf[half & 1u];
    FillU32(dest, partial_n, UINT32_MAX);
    while (vidx < c1 && HalfOf(vars[vidx].cm, origin) == half) {
      uint32_t is_ld = PgenVrtypeLd(vars[vidx].vrtype);
      if (!is_ld) {
        anchor_seed_half = (vars[vidx].role == kIbdRoleSeeded) ? half : UINT32_MAX;
      }
      if (vars[vidx].role == kIbdRoleSeeded) {
        if (is_ld) {
          uint32_t diff_ct = 0;
          reterr = PgenLdDiffCount(pg, vidx, &diff_ct);
          if (unlikely(reterr)) {
            fprintf(stderr, "pgen_ibd: difflist count failed at %s\n", vars[vidx].id);
            ArenaReset(arena, mark);
            return reterr;
          }
          if (diff_ct <= ld_lim && anchor_seed_half != UINT32_MAX &&
              (anchor_seed_half == half || anchor_seed_half + 1u == half)) {
            vars[vidx].role = kIbdRoleLdSkip;
            (*ld_skip_ct)++;
            vidx++;
            continue;
          }
        }
        UpdateMinhash(dest, hap_ct, perm_ct, vidx, vars[vidx].hap_alt);
      }
      vidx++;
    }
    if (half >= 1u && ext_ct) {
      double win0 = origin + (double)(half - 1u) * kIbdHalfCm;
      double win1 = origin + (double)(half + 1u) * kIbdHalfCm;
      reterr = EmitWindow(arena, buf[0], buf[1], hap_ct, perm_ct, cfg->band_ct,
                          sites, ext_ct, hap_bits, hap_word_ct, win0, win1,
                          segs, seg_ct, seg_cap);
      if (unlikely(reterr)) {
        ArenaReset(arena, mark);
        return reterr;
      }
    }
  }
  if (max_half == 0u && ext_ct) {
    reterr = EmitWindow(arena, buf[0], buf[1], hap_ct, perm_ct, cfg->band_ct,
                        sites, ext_ct, hap_bits, hap_word_ct,
                        origin, origin + kIbdHalfCm, segs, seg_ct, seg_cap);
  }
  if (vidx != c1 && !reterr) {
    fprintf(stderr, "pgen_ibd: unsorted genetic map near %s\n", vars[vidx].id);
    reterr = kPglRetInconsistentInput;
  }
  ArenaReset(arena, mark);
  return reterr;
}

PglErr IbdRun(const IbdConfig* cfg, IbdTrace* trace) {
  Arena arena;
  PgenFile pg;
  MapPoint* pts = 0;
  ChromMap* maps = 0;
  uint32_t map_ct = 0;
  uint32_t point_ct = 0;
  VarInfo* vars = 0;
  SampleInfo* samples = 0;
  uintptr_t* hap_slab = 0;
  uintptr_t* genovec = 0;
  uintptr_t* phasepresent = 0;
  uintptr_t* phaseinfo = 0;
  uint8_t* roles = 0;
  RawSeg* segs = 0;
  uint32_t seg_ct = 0;
  uint32_t seg_cap = 0;
  uint32_t ld_skip_ct = 0;
  uint32_t out_seg_ct = 0;
  uint32_t out_pair_ct = 0;
  uint32_t hap_word;
  uint32_t ld_lim;
  uint32_t vidx;
  uint32_t c0;
  uint32_t phased = 0;
  PglErr reterr = kPglRetSuccess;
  memset(&arena, 0, sizeof(arena));
  memset(&pg, 0, sizeof(pg));
  pg.anchor_vidx = UINT32_MAX;
  if (!cfg->pgen_path || !cfg->pvar_path || !cfg->psam_path || !cfg->map_path || !cfg->out_prefix) {
    return kPglRetInvalidCmdline;
  }
  if (cfg->perm_ct == 0u || cfg->band_ct == 0u || (cfg->perm_ct % cfg->band_ct) != 0u ||
      cfg->perm_ct / cfg->band_ct > 4u || cfg->min_length_cm <= 0.0 || cfg->mac_extend == 0u) {
    fprintf(stderr, "pgen_ibd: perm-ct must be a multiple of band-ct, with at most 4 hashes per band\n");
    return kPglRetInvalidCmdline;
  }
  reterr = PgenOpen(&pg, cfg->pgen_path);
  if (unlikely(reterr)) {
    fprintf(stderr, "pgen_ibd: cannot read %s (%s)\n", cfg->pgen_path, PglErrName(reterr));
    goto cleanup;
  }
  {
    uintptr_t cap = ((uintptr_t)32 << 20);
    cap += (uintptr_t)pg.variant_ct * 512u;
    cap += (uintptr_t)pg.variant_ct * (uintptr_t)pg.sample_ct;
    reterr = ArenaInit(&arena, cap);
    if (unlikely(reterr)) {
      goto cleanup;
    }
  }
  for (vidx = 0; vidx < pg.variant_ct; vidx++) {
    if (PgenVrtypePhase(pg.vrtypes[vidx])) {
      phased = 1;
      break;
    }
  }
  if (!phased) {
    fprintf(stderr, "pgen_ibd: pgen has no hardcall phase track; unphased files are rejected\n");
    reterr = kPglRetInconsistentInput;
    goto cleanup;
  }
  reterr = LoadMap(&arena, cfg->map_path, &pts, &point_ct, &maps, &map_ct);
  if (unlikely(reterr)) {
    goto cleanup;
  }
  reterr = LoadPsam(&arena, cfg->psam_path, pg.sample_ct, &samples);
  if (unlikely(reterr)) {
    goto cleanup;
  }
  reterr = LoadPvar(&arena, cfg->pvar_path, &pg, pts, maps, map_ct, &vars);
  if (unlikely(reterr)) {
    goto cleanup;
  }
  hap_word = BitCtToWordCt(pg.sample_ct * 2u);
  hap_slab = (uintptr_t*)ArenaAllocZero(&arena, (uintptr_t)pg.variant_ct * hap_word * sizeof(uintptr_t));
  genovec = (uintptr_t*)calloc(NypCtToWordCt(pg.sample_ct), sizeof(uintptr_t));
  phasepresent = (uintptr_t*)calloc(BitCtToWordCt(pg.sample_ct), sizeof(uintptr_t));
  phaseinfo = (uintptr_t*)calloc(BitCtToWordCt(pg.sample_ct), sizeof(uintptr_t));
  roles = (uint8_t*)malloc(pg.variant_ct);
  if (unlikely(!hap_slab || !genovec || !phasepresent || !phaseinfo || !roles)) {
    reterr = kPglRetNomem;
    goto cleanup;
  }
  reterr = BuildMasks(&pg, vars, cfg->mac_extend, cfg->maf_seed, hap_slab, hap_word,
                      genovec, phasepresent, phaseinfo);
  if (unlikely(reterr)) {
    goto cleanup;
  }
  ld_lim = cfg->ld_skip >= 0 ? (uint32_t)cfg->ld_skip : (pg.sample_ct / 100u);
  if (ld_lim < 1u) {
    ld_lim = 1u;
  }
  for (c0 = 0; c0 < pg.variant_ct; ) {
    uint32_t c1 = c0 + 1u;
    while (c1 < pg.variant_ct && strcmp(vars[c1].chrom, vars[c0].chrom) == 0) {
      c1++;
    }
    reterr = RunChrom(&pg, &arena, vars, c0, c1, cfg, ld_lim, &ld_skip_ct, &segs, &seg_ct, &seg_cap);
    if (unlikely(reterr)) {
      goto cleanup;
    }
    c0 = c1;
  }
  reterr = IbdStitchAndWrite(samples, pg.sample_ct, segs, seg_ct, cfg->min_length_cm, cfg->out_prefix,
                             &out_seg_ct, &out_pair_ct);
  if (unlikely(reterr)) {
    goto cleanup;
  }
  printf("pgen_ibd: wrote %u segment(s), %u pair(s), skipped %u tiny LD diff(s)\n",
         out_seg_ct, out_pair_ct, ld_skip_ct);
  if (trace) {
    for (vidx = 0; vidx < pg.variant_ct; vidx++) {
      roles[vidx] = vars[vidx].role;
    }
    trace->roles = roles;
    roles = 0;
    trace->variant_ct = pg.variant_ct;
    trace->sample_ct = pg.sample_ct;
    trace->ld_skip_ct = ld_skip_ct;
    trace->segment_ct = out_seg_ct;
    trace->pair_ct = out_pair_ct;
  }
  cleanup:
  free(roles);
  free(genovec);
  free(phasepresent);
  free(phaseinfo);
  free(segs);
  free(pts);
  free(maps);
  PgenClose(&pg);
  ArenaDestroy(&arena);
  return reterr;
}
