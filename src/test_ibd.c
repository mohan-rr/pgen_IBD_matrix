#include "ibd_internal.h"

#include <stdio.h>

static int g_fails = 0;

static void Fail(const char* msg) {
  fprintf(stderr, "FAIL: %s\n", msg);
  g_fails++;
}

static uint32_t WritePgen(const char* path, uint32_t sample_ct, const unsigned char* vrtypes,
                          const unsigned char* const* recs, const uint32_t* lens, uint32_t variant_ct) {
  FILE* fp;
  unsigned char prefix[12];
  uint64_t offset;
  uint32_t vidx;
  uint32_t i;
  offset = 12u + 8u + (uint64_t)variant_ct + (uint64_t)variant_ct;
  fp = fopen(path, "wb");
  if (!fp) {
    return 0;
  }
  prefix[0] = 0x6c;
  prefix[1] = 0x1b;
  prefix[2] = 0x10;
  for (i = 0; i < 4u; i++) {
    prefix[3u + i] = (unsigned char)((variant_ct >> (8u * i)) & 255u);
    prefix[7u + i] = (unsigned char)((sample_ct >> (8u * i)) & 255u);
  }
  prefix[11] = 0x04;
  if (fwrite(prefix, 1, 12, fp) != 12) {
    fclose(fp);
    return 0;
  }
  for (i = 0; i < 8u; i++) {
    unsigned char b = (unsigned char)((offset >> (8u * i)) & 255u);
    if (fwrite(&b, 1, 1, fp) != 1) {
      fclose(fp);
      return 0;
    }
  }
  if (fwrite(vrtypes, 1, variant_ct, fp) != variant_ct) {
    fclose(fp);
    return 0;
  }
  for (vidx = 0; vidx < variant_ct; vidx++) {
    unsigned char b = (unsigned char)lens[vidx];
    if (lens[vidx] > 255u || fwrite(&b, 1, 1, fp) != 1) {
      fclose(fp);
      return 0;
    }
  }
  for (vidx = 0; vidx < variant_ct; vidx++) {
    if (fwrite(recs[vidx], 1, lens[vidx], fp) != lens[vidx]) {
      fclose(fp);
      return 0;
    }
  }
  if (fclose(fp) != 0) {
    return 0;
  }
  return 1;
}

static void ExpectGeno(const char* name, const uintptr_t* genovec, const uintptr_t* phasepresent,
                       const uintptr_t* phaseinfo, const uint8_t* geno, uint32_t present_mask, uint32_t info_mask) {
  uint32_t sample_idx;
  char buf[128];
  for (sample_idx = 0; sample_idx < 4u; sample_idx++) {
    if (GetGeno(genovec, sample_idx) != geno[sample_idx]) {
      snprintf(buf, sizeof(buf), "%s geno sample %u got %u want %u", name, sample_idx,
               GetGeno(genovec, sample_idx), geno[sample_idx]);
      Fail(buf);
    }
    if (GetBit(phasepresent, sample_idx) != ((present_mask >> sample_idx) & 1u)) {
      snprintf(buf, sizeof(buf), "%s phasepresent sample %u", name, sample_idx);
      Fail(buf);
    }
    if (GetBit(phaseinfo, sample_idx) != ((info_mask >> sample_idx) & 1u)) {
      snprintf(buf, sizeof(buf), "%s phaseinfo sample %u", name, sample_idx);
      Fail(buf);
    }
  }
}

static void TestVarint(void) {
  unsigned char buf[2];
  uint32_t count = 0;
  uint32_t consumed = 0;
  buf[0] = 0x80;
  buf[1] = 0x01;
  if (!PgenReadDifflistCount(buf, buf + 2, &count, &consumed) || count != 128u || consumed != 2u) {
    Fail("multi-byte difflist varint");
  }
}

static void TestHandcrafted(const char* dir) {
  /* Anchor raw [0,0,2,1], phase on sample 3 unswapped.
     Type-2 and type-3 both decode to [2,0,0,1]. Record length 5 is not count*3. */
  static const unsigned char rec0[] = {0x60, 0x00};
  static const unsigned char rec1[] = {0x02, 0x00, 0x02, 0x02, 0x00};
  static const unsigned char rec2[] = {0x01, 0x01, 0x02, 0x00};
  static const unsigned char rec3[] = {0x01, 0x03, 0x00, 0x06};
  static const unsigned char rec4[] = {0x01, 0x03, 0x00};
  static const unsigned char vrtypes[] = {0x10, 0x12, 0x13, 0x11, 0x10};
  static const uint32_t lens[] = {2, 5, 4, 4, 3};
  const unsigned char* recs[] = {rec0, rec1, rec2, rec3, rec4};
  char path[512];
  PgenFile pg;
  uintptr_t genovec[2];
  uintptr_t phasepresent[2];
  uintptr_t phaseinfo[2];
  uintptr_t genovec_b[2];
  uintptr_t present_b[2];
  uintptr_t info_b[2];
  uint32_t phase_ct = 0;
  uint32_t count = 0;
  uint8_t anchor_g[4] = {0, 0, 2, 1};
  uint8_t target_g[4] = {2, 0, 0, 1};
  uint8_t sparse_g[4] = {1, 1, 0, 0};
  uint8_t explicit_g[4] = {1, 0, 0, 0};
  PglErr reterr;
  snprintf(path, sizeof(path), "%s/unit.pgen", dir);
  if (!WritePgen(path, 4, vrtypes, recs, lens, 5)) {
    Fail("write unit pgen");
    return;
  }
  memset(&pg, 0, sizeof(pg));
  memset(genovec, 0, sizeof(genovec));
  memset(phasepresent, 0, sizeof(phasepresent));
  memset(phaseinfo, 0, sizeof(phaseinfo));
  memset(genovec_b, 0, sizeof(genovec_b));
  memset(present_b, 0, sizeof(present_b));
  memset(info_b, 0, sizeof(info_b));
  reterr = PgenOpen(&pg, path);
  if (reterr) {
    Fail("open unit pgen");
    return;
  }
  reterr = PgenLdDiffCount(&pg, 1, &count);
  if (reterr || count != 2u || PgenRecLen(&pg, 1) != 5u) {
    Fail("type-2 difflist count is not record_bytes/3");
  }
  reterr = PgenGetP(&pg, 0, genovec, phasepresent, phaseinfo, &phase_ct);
  if (reterr) {
    Fail("get anchor");
  } else {
    ExpectGeno("anchor", genovec, phasepresent, phaseinfo, anchor_g, 1u << 3, 0);
  }
  reterr = PgenGetP(&pg, 1, genovec, phasepresent, phaseinfo, &phase_ct);
  if (reterr) {
    Fail("get type2");
  } else {
    ExpectGeno("type2", genovec, phasepresent, phaseinfo, target_g, 1u << 3, 0);
  }
  reterr = PgenGetP(&pg, 2, genovec_b, present_b, info_b, &phase_ct);
  if (reterr) {
    Fail("get type3");
  } else {
    ExpectGeno("type3", genovec_b, present_b, info_b, target_g, 1u << 3, 0);
    if (memcmp(genovec, genovec_b, sizeof(genovec)) != 0 ||
        memcmp(phasepresent, present_b, sizeof(phasepresent)) != 0 ||
        memcmp(phaseinfo, info_b, sizeof(info_b)) != 0) {
      Fail("type2 and type3 haplotype bits differ");
    }
  }
  reterr = PgenGetP(&pg, 3, genovec, phasepresent, phaseinfo, &phase_ct);
  if (reterr) {
    Fail("get sparse");
  } else {
    ExpectGeno("sparse", genovec, phasepresent, phaseinfo, sparse_g, 0x3u, 0x3u);
  }
  reterr = PgenGetP(&pg, 4, genovec, phasepresent, phaseinfo, &phase_ct);
  if (reterr) {
    Fail("get explicit phase");
  } else {
    ExpectGeno("explicit-phase", genovec, phasepresent, phaseinfo, explicit_g, 0x1u, 0);
  }
  PgenClose(&pg);
}

static void HapsToCall(const uint8_t h0[4], const uint8_t h1[4], uint8_t geno[4], uint8_t swapped[4]) {
  uint32_t i;
  for (i = 0; i < 4u; i++) {
    if (h0[i] == 0u && h1[i] == 0u) {
      geno[i] = 0;
      swapped[i] = 0;
    } else if (h0[i] == 1u && h1[i] == 1u) {
      geno[i] = 2;
      swapped[i] = 0;
    } else if (h0[i] == 0u && h1[i] == 1u) {
      geno[i] = 1;
      swapped[i] = 0;
    } else {
      geno[i] = 1;
      swapped[i] = 1;
    }
  }
}

static unsigned char PhaseByte(const uint8_t geno[4], const uint8_t swapped[4]) {
  unsigned char phase = 0;
  uint32_t bit = 1;
  uint32_t i;
  for (i = 0; i < 4u; i++) {
    if (geno[i] == 1u) {
      if (swapped[i]) {
        phase = (unsigned char)(phase | (unsigned char)(1u << bit));
      }
      bit++;
    }
  }
  return phase;
}

static unsigned char GenoByte(const uint8_t geno[4]) {
  unsigned char packed = 0;
  uint32_t i;
  for (i = 0; i < 4u; i++) {
    packed = (unsigned char)(packed | (unsigned char)(geno[i] << (2u * i)));
  }
  return packed;
}

static uint8_t InvertGenoByte(uint8_t geno) {
  if (geno == 0u) {
    return 2u;
  }
  if (geno == 2u) {
    return 0u;
  }
  return geno;
}

static uint32_t AppendDifflist(unsigned char* dst, const uint32_t* ids, const uint8_t* genos, uint32_t difflist_len) {
  uint32_t nwrite = 0;
  uint32_t i;
  uint32_t geno_bytes;
  dst[nwrite++] = (unsigned char)difflist_len;
  if (difflist_len == 0u) {
    return nwrite;
  }
  dst[nwrite++] = (unsigned char)ids[0];
  geno_bytes = (difflist_len + 3u) / 4u;
  memset(dst + nwrite, 0, geno_bytes);
  for (i = 0; i < difflist_len; i++) {
    dst[nwrite + i / 4u] = (unsigned char)(dst[nwrite + i / 4u] | (unsigned char)((genos[i] & 3u) << (2u * (i % 4u))));
  }
  nwrite += geno_bytes;
  for (i = 1; i < difflist_len; i++) {
    dst[nwrite++] = (unsigned char)(ids[i] - ids[i - 1u]);
  }
  return nwrite;
}

static void Pattern(uint32_t vidx, uint8_t h0[4], uint8_t h1[4]) {
  uint32_t s;
  for (s = 0; s < 4u; s++) {
    h0[s] = 0;
    h1[s] = 0;
  }
  if (vidx == 7u) {
    h0[0] = 1;
    return;
  }
  if (vidx == 12u) {
    h0[0] = 1;
    h0[1] = 1;
    return;
  }
  if (vidx == 8u || vidx == 4u || vidx == 6u || vidx == 10u || vidx == 0u || vidx == 2u || vidx == 14u || vidx == 16u) {
    if (vidx == 0u || vidx == 2u || vidx == 14u || vidx == 16u) {
      h0[0] = 1;
      h1[2] = 1;
      return;
    }
    h0[0] = 1;
    h0[1] = 1;
    h1[2] = 1;
    return;
  }
  h0[0] = 1;
  if (vidx == 1u || vidx == 3u || vidx == 13u || vidx == 15u) {
    h1[3] = 1;
    return;
  }
  h0[1] = 1;
  h1[3] = 1;
}

static uint32_t EncodeVariant(uint32_t vidx, uint8_t geno[][4], unsigned char* dst, unsigned char* vrtype) {
  uint8_t h0[4];
  uint8_t h1[4];
  uint8_t swapped[4];
  unsigned char phase;
  Pattern(vidx, h0, h1);
  HapsToCall(h0, h1, geno[vidx], swapped);
  phase = PhaseByte(geno[vidx], swapped);
  if (vidx == 10u) {
    unsigned char bits = 0;
    uint32_t s;
    for (s = 0; s < 4u; s++) {
      if (geno[vidx][s] == 1u) {
        bits = (unsigned char)(bits | (unsigned char)(1u << s));
      }
    }
    dst[0] = 0x01;
    dst[1] = bits;
    dst[2] = 0x00;
    dst[3] = phase;
    *vrtype = 0x11;
    return 4;
  }
  if (vidx == 9u || vidx == 12u) {
    uint32_t anchor = (vidx == 9u) ? 8u : 11u;
    uint32_t ids[4];
    uint8_t values[4];
    uint32_t nlist = 0;
    uint32_t s;
    uint32_t nwrite;
    for (s = 0; s < 4u; s++) {
      uint8_t stored = geno[vidx][s];
      if (vidx == 9u) {
        stored = InvertGenoByte(stored);
      }
      if (stored != geno[anchor][s]) {
        ids[nlist] = s;
        values[nlist] = stored;
        nlist++;
      }
    }
    nwrite = AppendDifflist(dst, ids, values, nlist);
    dst[nwrite++] = phase;
    *vrtype = (vidx == 9u) ? 0x13 : 0x12;
    return nwrite;
  }
  dst[0] = GenoByte(geno[vidx]);
  dst[1] = phase;
  *vrtype = 0x10;
  return 2;
}

static uint32_t WriteTextFiles(const char* dir) {
  char path[512];
  FILE* fp;
  uint32_t vidx;
  snprintf(path, sizeof(path), "%s/test.psam", dir);
  fp = fopen(path, "wb");
  if (!fp) {
    return 0;
  }
  fputs("#FID IID\nF A\nF B\nF C\nF D\n", fp);
  fclose(fp);
  snprintf(path, sizeof(path), "%s/test.pvar", dir);
  fp = fopen(path, "wb");
  if (!fp) {
    return 0;
  }
  fputs("#CHROM POS ID REF ALT\n", fp);
  for (vidx = 0; vidx < 17u; vidx++) {
    const char* alt = (vidx == 8u) ? "ACGTACGT" : "C";
    fprintf(fp, "1 %u v%u A %s\n", vidx * 10000u, vidx, alt);
  }
  fclose(fp);
  snprintf(path, sizeof(path), "%s/test.map", dir);
  fp = fopen(path, "wb");
  if (!fp) {
    return 0;
  }
  for (vidx = 0; vidx < 17u; vidx++) {
    fprintf(fp, "1 v%u %.1f %u\n", vidx, (double)vidx * 0.5, vidx * 10000u);
  }
  fclose(fp);
  return 1;
}

static uint32_t CountDataLines(const char* path, char lines[][256], uint32_t max_lines) {
  FILE* fp = fopen(path, "rb");
  char buf[256];
  uint32_t n = 0;
  if (!fp) {
    return UINT32_MAX;
  }
  while (fgets(buf, sizeof(buf), fp)) {
    size_t len = strlen(buf);
    if (len && buf[len - 1u] == '\n') {
      buf[--len] = 0;
    }
    if (buf[0] == 0 || buf[0] == '#') {
      continue;
    }
    if (n < max_lines) {
      snprintf(lines[n], 256, "%s", buf);
    }
    n++;
  }
  fclose(fp);
  return n;
}

static void TestIntegration(const char* dir, const char* bin) {
  enum { kN = 17 };
  unsigned char recs_store[kN][16];
  unsigned char vrtypes[kN];
  uint32_t lens[kN];
  const unsigned char* rec_ptrs[kN];
  uint8_t geno[kN][4];
  char path[640];
  char out_prefix[512];
  char cli_prefix[512];
  IbdConfig cfg;
  IbdTrace trace;
  PgenFile pg;
  uintptr_t genovec[2];
  uintptr_t phasepresent[2];
  uintptr_t phaseinfo[2];
  uint32_t phase_ct = 0;
  uint32_t vidx;
  uint32_t diff_ct = 0;
  char lines[8][256];
  uint32_t nline;
  PglErr reterr;
  uint8_t in2[4] = {1, 1, 0, 1};
  uint8_t in1[4] = {1, 1, 1, 0};
  uint8_t endg[4] = {1, 1, 0, 0};
  for (vidx = 0; vidx < kN; vidx++) {
    lens[vidx] = EncodeVariant(vidx, geno, recs_store[vidx], &vrtypes[vidx]);
    rec_ptrs[vidx] = recs_store[vidx];
  }
  if (!WriteTextFiles(dir)) {
    Fail("write text inputs");
    return;
  }
  snprintf(path, sizeof(path), "%s/test.pgen", dir);
  if (!WritePgen(path, 4, vrtypes, rec_ptrs, lens, kN)) {
    Fail("write integration pgen");
    return;
  }
  memset(&pg, 0, sizeof(pg));
  if (PgenOpen(&pg, path)) {
    Fail("reopen integration pgen");
    return;
  }
  if (PgenGetP(&pg, 9, genovec, phasepresent, phaseinfo, &phase_ct)) {
    Fail("decode type3 integration");
  } else {
    ExpectGeno("int-type3", genovec, phasepresent, phaseinfo, in2, 0xbu, 0x3u);
  }
  if (PgenGetP(&pg, 10, genovec, phasepresent, phaseinfo, &phase_ct)) {
    Fail("decode sparse integration");
  } else {
    ExpectGeno("int-sparse", genovec, phasepresent, phaseinfo, in1, 0x7u, 0x3u);
  }
  if (PgenGetP(&pg, 12, genovec, phasepresent, phaseinfo, &phase_ct)) {
    Fail("decode boundary ld");
  } else {
    ExpectGeno("int-ld", genovec, phasepresent, phaseinfo, endg, 0x3u, 0x3u);
  }
  if (PgenLdDiffCount(&pg, 12, &diff_ct) || diff_ct != 1u) {
    Fail("boundary ld difflist count");
  }
  PgenClose(&pg);

  IbdConfigDefaults(&cfg);
  snprintf(out_prefix, sizeof(out_prefix), "%s/out", dir);
  cfg.pgen_path = path;
  snprintf(path, sizeof(path), "%s/test.pvar", dir);
  cfg.pvar_path = path;
  {
    static char pgen_path[512];
    static char pvar_path[512];
    static char psam_path[512];
    static char map_path[512];
    snprintf(pgen_path, sizeof(pgen_path), "%s/test.pgen", dir);
    snprintf(pvar_path, sizeof(pvar_path), "%s/test.pvar", dir);
    snprintf(psam_path, sizeof(psam_path), "%s/test.psam", dir);
    snprintf(map_path, sizeof(map_path), "%s/test.map", dir);
    cfg.pgen_path = pgen_path;
    cfg.pvar_path = pvar_path;
    cfg.psam_path = psam_path;
    cfg.map_path = map_path;
    cfg.out_prefix = out_prefix;
    memset(&trace, 0, sizeof(trace));
    reterr = IbdRun(&cfg, &trace);
    if (reterr) {
      fprintf(stderr, "IbdRun failed: %s\n", PglErrName(reterr));
      Fail("IbdRun");
      return;
    }
  }
  if (trace.ld_skip_ct != 1u || !trace.roles || trace.roles[12] != kIbdRoleLdSkip) {
    Fail("tiny LD diff was not skipped exactly once");
  }
  if (!trace.roles || trace.roles[7] != kIbdRoleDrop || trace.roles[8] != kIbdRoleDrop) {
    Fail("singleton or long indel entered a mask");
  }
  if (!trace.roles || trace.roles[4] != kIbdRoleSeeded || trace.roles[9] != kIbdRoleSeeded) {
    Fail("expected seed sites were not seeded");
  }
  if (trace.segment_ct != 1u || trace.pair_ct != 1u) {
    fprintf(stderr, "segments %u pairs %u\n", trace.segment_ct, trace.pair_ct);
    Fail("expected one tract and one pair");
  }
  snprintf(path, sizeof(path), "%s.ibd", out_prefix);
  nline = CountDataLines(path, lines, 8);
  if (nline != 1u) {
    Fail("matrix should contain one pair");
  } else {
    char s1[64], h1[64], s2[64], h2[64];
    double shared = 0;
    unsigned nseg = 0;
    if (sscanf(lines[0], "%63s %63s %63s %63s %lf %u", s1, h1, s2, h2, &shared, &nseg) != 6 ||
        strcmp(s1, "A") != 0 || strcmp(h1, "SAMPLE_0") != 0 ||
        strcmp(s2, "B") != 0 || strcmp(h2, "SAMPLE_0") != 0 || nseg != 1u ||
        shared < 3.999 || shared > 4.001) {
      fprintf(stderr, "matrix line: %s\n", lines[0]);
      Fail("sparse matrix row");
    }
  }
  snprintf(path, sizeof(path), "%s.segments", out_prefix);
  nline = CountDataLines(path, lines, 8);
  if (nline != 1u) {
    Fail("segment table should contain one row");
  } else {
    char f1[32], i1[32], h1[32], f2[32], i2[32], h2[32], chrom[32], id0[32], id1[32];
    unsigned bp0 = 0;
    unsigned bp1 = 0;
    double cm = 0;
    double sim = 0;
    if (sscanf(lines[0], "%31s %31s %31s %31s %31s %31s %31s %u %u %31s %31s %lf %lf",
               f1, i1, h1, f2, i2, h2, chrom, &bp0, &bp1, id0, id1, &cm, &sim) != 13 ||
        strcmp(i1, "A") != 0 || strcmp(i2, "B") != 0 || strcmp(h1, "SAMPLE_0") != 0 ||
        strcmp(id0, "v4") != 0 || strcmp(id1, "v12") != 0 || bp0 != 40000u || bp1 != 120000u ||
        cm < 3.999 || cm > 4.001 || sim < 0.999) {
      fprintf(stderr, "segment line: %s\n", lines[0]);
      Fail("segment boundary");
    }
  }
  free(trace.roles);
  if (bin && bin[0]) {
    char cmd[1024];
    snprintf(cli_prefix, sizeof(cli_prefix), "%s/cli", dir);
    snprintf(cmd, sizeof(cmd),
             "%s --pgen %s/test.pgen --pvar %s/test.pvar --psam %s/test.psam --map %s/test.map --out %s",
             bin, dir, dir, dir, dir, cli_prefix);
    if (system(cmd) != 0) {
      Fail("CLI run");
    } else {
      snprintf(path, sizeof(path), "%s.ibd", cli_prefix);
      nline = CountDataLines(path, lines, 8);
      if (nline != 1u) {
        Fail("CLI matrix");
      }
    }
  }
}

int main(int argc, char** argv) {
  const char* dir = "/tmp/pgen_ibd_test";
  const char* bin = argc > 1 ? argv[1] : "";
  if (system("rm -rf /tmp/pgen_ibd_test && mkdir -p /tmp/pgen_ibd_test") != 0) {
    fprintf(stderr, "cannot create test dir\n");
    return 1;
  }
  TestVarint();
  TestHandcrafted(dir);
  TestIntegration(dir, bin);
  if (g_fails) {
    fprintf(stderr, "%d failure(s)\n", g_fails);
    return 1;
  }
  printf("all tests passed\n");
  return 0;
}
