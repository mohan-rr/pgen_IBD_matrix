#include "pgen_read.h"

#include <inttypes.h>

static uint32_t LoadLeU32(const unsigned char* ptr, uint32_t nbytes) {
  uint32_t value = 0;
  uint32_t byte_idx;
  for (byte_idx = 0; byte_idx < nbytes; byte_idx++) {
    value |= (uint32_t)ptr[byte_idx] << (8u * byte_idx);
  }
  return value;
}

static uint64_t LoadLeU64(const unsigned char* ptr) {
  uint64_t value = 0;
  uint32_t byte_idx;
  for (byte_idx = 0; byte_idx < 8u; byte_idx++) {
    value |= (uint64_t)ptr[byte_idx] << (8u * byte_idx);
  }
  return value;
}

static PglErr ReadExact(FILE* fp, void* dst, uint64_t nbytes) {
  unsigned char* ptr = (unsigned char*)dst;
  while (nbytes) {
    size_t chunk = nbytes > (1u << 20) ? (size_t)(1u << 20) : (size_t)nbytes;
    size_t got = fread(ptr, 1, chunk, fp);
    if (got != chunk) {
      return kPglRetReadFail;
    }
    ptr += got;
    nbytes -= got;
  }
  return kPglRetSuccess;
}

static PglErr SeekTo(FILE* fp, uint64_t pos) {
  if (fseeko(fp, (off_t)pos, SEEK_SET) != 0) {
    return kPglRetReadFail;
  }
  return kPglRetSuccess;
}

void PgenClose(PgenFile* pg) {
  if (!pg) {
    return;
  }
  if (pg->fp) {
    fclose(pg->fp);
    pg->fp = 0;
  }
  free(pg->vrtypes);
  free(pg->rec_lens);
  free(pg->rec_fpos);
  free(pg->recbuf);
  free(pg->anchor_genovec);
  pg->vrtypes = 0;
  pg->rec_lens = 0;
  pg->rec_fpos = 0;
  pg->recbuf = 0;
  pg->anchor_genovec = 0;
  pg->recbuf_cap = 0;
  pg->variant_ct = 0;
  pg->sample_ct = 0;
  pg->anchor_vidx = UINT32_MAX;
  pg->anchor_nword = 0;
}

static PglErr EnsureRecbuf(PgenFile* pg, uint32_t nbytes) {
  unsigned char* fresh;
  if (nbytes <= pg->recbuf_cap) {
    return kPglRetSuccess;
  }
  fresh = (unsigned char*)realloc(pg->recbuf, nbytes);
  if (unlikely(!fresh)) {
    return kPglRetNomem;
  }
  pg->recbuf = fresh;
  pg->recbuf_cap = nbytes;
  return kPglRetSuccess;
}

static PglErr LoadRecord(PgenFile* pg, uint32_t vidx, const unsigned char** start, const unsigned char** end) {
  PglErr reterr;
  uint32_t rec_len = pg->rec_lens[vidx];
  if (rec_len > 0x7fffffffu) {
    return kPglRetVarRecordTooLarge;
  }
  reterr = EnsureRecbuf(pg, rec_len);
  if (unlikely(reterr)) {
    return reterr;
  }
  reterr = SeekTo(pg->fp, pg->rec_fpos[vidx]);
  if (unlikely(reterr)) {
    return reterr;
  }
  reterr = ReadExact(pg->fp, pg->recbuf, rec_len);
  if (unlikely(reterr)) {
    return reterr;
  }
  *start = pg->recbuf;
  *end = pg->recbuf + rec_len;
  return kPglRetSuccess;
}

static PglErr ApplyDifflist(const unsigned char* end, const unsigned char** pp,
                            uint32_t sample_ct, uintptr_t* genovec) {
  uint32_t difflist_len = 0;
  uint32_t consumed = 0;
  uint32_t group_ct;
  uint32_t id_bytes;
  uint32_t index_bytes;
  uint32_t geno_bytes;
  uint32_t filled = 0;
  uint32_t group_idx;
  const unsigned char* leaders;
  const unsigned char* genos;
  if (!PgenReadDifflistCount(*pp, end, &difflist_len, &consumed)) {
    return kPglRetMalformedInput;
  }
  *pp += consumed;
  if (difflist_len == 0u) {
    return kPglRetSuccess;
  }
  if (difflist_len > sample_ct) {
    return kPglRetMalformedInput;
  }
  group_ct = DivUpU32(difflist_len, kPglDifflistGroupSize);
  id_bytes = BytesToRepresentNzU32(sample_ct);
  index_bytes = group_ct * (id_bytes + 1u) - 1u;
  if (!PtrRoom(*pp, end, index_bytes)) {
    return kPglRetMalformedInput;
  }
  leaders = *pp;
  *pp += index_bytes;
  geno_bytes = DivUpU32(difflist_len, 4u);
  if (!PtrRoom(*pp, end, geno_bytes)) {
    return kPglRetMalformedInput;
  }
  genos = *pp;
  *pp += geno_bytes;
  for (group_idx = 0; group_idx < group_ct; group_idx++) {
    uint32_t sample_id = LoadLeU32(leaders + (uintptr_t)group_idx * id_bytes, id_bytes);
    uint32_t in_group = difflist_len - filled;
    uint32_t j;
    if (in_group > kPglDifflistGroupSize) {
      in_group = kPglDifflistGroupSize;
    }
    for (j = 0; j < in_group; j++) {
      uint32_t geno;
      if (j > 0u) {
        uint32_t delta = 0;
        uint32_t delta_bytes = 0;
        if (!PgenReadDifflistCount(*pp, end, &delta, &delta_bytes) || delta == 0u) {
          return kPglRetMalformedInput;
        }
        *pp += delta_bytes;
        if (sample_id > sample_ct - 1u - delta) {
          return kPglRetMalformedInput;
        }
        sample_id += delta;
      }
      if (sample_id >= sample_ct) {
        return kPglRetMalformedInput;
      }
      geno = (genos[filled / 4u] >> (2u * (filled % 4u))) & 3u;
      SetGeno(genovec, sample_id, geno);
      filled++;
    }
  }
  if (filled != difflist_len) {
    return kPglRetMalformedInput;
  }
  return kPglRetSuccess;
}

static PglErr DecodeOnebit(const unsigned char* end, const unsigned char** pp,
                           uint32_t sample_ct, uintptr_t* genovec) {
  uint32_t prefix_bytes = (sample_ct + 15u) / 8u;
  uint32_t code;
  uint32_t base_geno;
  uint32_t delta;
  uint32_t sample_idx;
  const unsigned char* bits;
  if (!PtrRoom(*pp, end, prefix_bytes)) {
    return kPglRetMalformedInput;
  }
  code = (*pp)[0];
  base_geno = code >> 2;
  delta = code & 3u;
  if (base_geno > 3u || base_geno + delta > 3u) {
    return kPglRetMalformedInput;
  }
  bits = *pp + 1;
  *pp += prefix_bytes;
  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
    uint32_t bit = (bits[sample_idx / 8u] >> (sample_idx % 8u)) & 1u;
    SetGeno(genovec, sample_idx, base_geno + (bit ? delta : 0u));
  }
  ZeroTrailingNyps(sample_ct, genovec);
  return ApplyDifflist(end, pp, sample_ct, genovec);
}

static PglErr DecodeNonLdMain(uint32_t vrtype, const unsigned char* end, const unsigned char** pp,
                              uint32_t sample_ct, uintptr_t* genovec) {
  uint32_t main_type = PgenVrtypeMain(vrtype);
  uint32_t nword = NypCtToWordCt(sample_ct);
  memset(genovec, 0, (size_t)nword * sizeof(uintptr_t));
  if (main_type == PGEN_VR_RAW) {
    uint32_t nbytes = NypCtToByteCt(sample_ct);
    if (!PtrRoom(*pp, end, nbytes)) {
      return kPglRetMalformedInput;
    }
    memcpy(genovec, *pp, nbytes);
    *pp += nbytes;
    ZeroTrailingNyps(sample_ct, genovec);
    return kPglRetSuccess;
  }
  if (main_type == PGEN_VR_SPARSE) {
    return DecodeOnebit(end, pp, sample_ct, genovec);
  }
  if (main_type == PGEN_VR_ALLREF) {
    return kPglRetSuccess;
  }
  if (main_type == PGEN_VR_DIFF0 || main_type == PGEN_VR_DIFF2 || main_type == PGEN_VR_DIFF3) {
    GenovecFill(sample_ct, main_type & 3u, genovec);
    return ApplyDifflist(end, pp, sample_ct, genovec);
  }
  return kPglRetMalformedInput;
}

static uint32_t PopcountBytes(const unsigned char* ptr, uint32_t nbytes) {
  uint32_t sum = 0;
  uint32_t byte_idx;
  for (byte_idx = 0; byte_idx < nbytes; byte_idx++) {
    sum += (uint32_t)__builtin_popcount((unsigned)ptr[byte_idx]);
  }
  return sum;
}

static PglErr DecodePhase(uint32_t vrtype, const unsigned char* end, const unsigned char** pp,
                          uint32_t sample_ct, const uintptr_t* genovec,
                          uintptr_t* phasepresent, uintptr_t* phaseinfo,
                          uint32_t* phasepresent_ct) {
  uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  uint32_t het_ct = 0;
  uint32_t sample_idx;
  memset(phasepresent, 0, (size_t)sample_ctl * sizeof(uintptr_t));
  memset(phaseinfo, 0, (size_t)sample_ctl * sizeof(uintptr_t));
  *phasepresent_ct = 0;
  if (!PgenVrtypePhase(vrtype)) {
    return kPglRetSuccess;
  }
  for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
    if (GetGeno(genovec, sample_idx) == 1u) {
      het_ct++;
    }
  }
  if (het_ct == 0u) {
    return kPglRetMalformedInput;
  }
  if (!PtrRoom(*pp, end, 1u)) {
    return kPglRetMalformedInput;
  }
  if (((*pp)[0] & 1u) == 0u) {
    uint32_t nbytes = 1u + het_ct / 8u;
    uint32_t bit_i = 1u;
    if (!PtrRoom(*pp, end, nbytes)) {
      return kPglRetMalformedInput;
    }
    for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
      if (GetGeno(genovec, sample_idx) == 1u) {
        uint32_t swapped = ((*pp)[bit_i / 8u] >> (bit_i % 8u)) & 1u;
        SetBit(phasepresent, sample_idx);
        if (swapped) {
          SetBit(phaseinfo, sample_idx);
        }
        bit_i++;
      }
    }
    *pp += nbytes;
    *phasepresent_ct = het_ct;
    return kPglRetSuccess;
  }
  {
    uint32_t first_bytes = 1u + het_ct / 8u;
    uint32_t raw_present;
    uint32_t info_bytes;
    uint32_t bit_i = 1u;
    uint32_t info_bit = 0;
    const unsigned char* info;
    if (!PtrRoom(*pp, end, first_bytes)) {
      return kPglRetMalformedInput;
    }
    raw_present = PopcountBytes(*pp, first_bytes);
    if (raw_present == 0u) {
      return kPglRetMalformedInput;
    }
    raw_present -= 1u;
    if (raw_present == 0u) {
      return kPglRetMalformedInput;
    }
    info = *pp + first_bytes;
    info_bytes = DivUpU32(raw_present, 8u);
    if (!PtrRoom(info, end, info_bytes)) {
      return kPglRetMalformedInput;
    }
    for (sample_idx = 0; sample_idx < sample_ct; sample_idx++) {
      if (GetGeno(genovec, sample_idx) == 1u) {
        uint32_t present = ((*pp)[bit_i / 8u] >> (bit_i % 8u)) & 1u;
        bit_i++;
        if (present) {
          uint32_t swapped = (info[info_bit / 8u] >> (info_bit % 8u)) & 1u;
          SetBit(phasepresent, sample_idx);
          if (swapped) {
            SetBit(phaseinfo, sample_idx);
          }
          info_bit++;
        }
      }
    }
    *pp = info + info_bytes;
    *phasepresent_ct = raw_present;
    return kPglRetSuccess;
  }
}

static PglErr FindLdAnchor(const PgenFile* pg, uint32_t vidx, uint32_t* anchor_vidx) {
  uint32_t block_start = vidx & ~(kPglVblockSize - 1u);
  uint32_t anchor = vidx;
  if (!PgenVrtypeLd(pg->vrtypes[vidx])) {
    return kPglRetMalformedInput;
  }
  while (anchor > block_start && PgenVrtypeLd(pg->vrtypes[anchor])) {
    anchor--;
  }
  if (PgenVrtypeLd(pg->vrtypes[anchor])) {
    return kPglRetMalformedInput;
  }
  *anchor_vidx = anchor;
  return kPglRetSuccess;
}

static PglErr DecodeAnchor(PgenFile* pg, uint32_t anchor_vidx) {
  const unsigned char* start;
  const unsigned char* end;
  const unsigned char* ptr;
  PglErr reterr;
  uint32_t vrtype;
  if (pg->anchor_vidx == anchor_vidx) {
    return kPglRetSuccess;
  }
  vrtype = pg->vrtypes[anchor_vidx];
  if (PgenVrtypeLd(vrtype) || PgenVrtypeMulti(vrtype)) {
    return kPglRetNotYetSupported;
  }
  reterr = LoadRecord(pg, anchor_vidx, &start, &end);
  if (unlikely(reterr)) {
    return reterr;
  }
  ptr = start;
  reterr = DecodeNonLdMain(vrtype, end, &ptr, pg->sample_ct, pg->anchor_genovec);
  if (unlikely(reterr)) {
    return reterr;
  }
  pg->anchor_vidx = anchor_vidx;
  return kPglRetSuccess;
}

PglErr PgenOpen(PgenFile* pg, const char* path) {
  unsigned char prefix[12];
  unsigned char* block_off_buf = 0;
  PglErr reterr = kPglRetSuccess;
  uint32_t variant_ct;
  uint32_t sample_ct;
  uint32_t control;
  uint32_t low;
  uint32_t vrec_len_byte_ct;
  uint32_t vrtype_8bit;
  uint32_t allele_byte_ct;
  uint32_t nonref_mode;
  uint32_t block_ct;
  uint32_t block_idx;
  uint32_t vidx;
  uint64_t header_pos;
  memset(pg, 0, sizeof(*pg));
  pg->anchor_vidx = UINT32_MAX;
  pg->fp = fopen(path, "rb");
  if (unlikely(!pg->fp)) {
    reterr = kPglRetOpenFail;
    goto cleanup;
  }
  reterr = ReadExact(pg->fp, prefix, 12);
  if (unlikely(reterr)) {
    goto cleanup;
  }
  if (prefix[0] != 0x6cu || prefix[1] != 0x1bu) {
    reterr = kPglRetMalformedInput;
    goto cleanup;
  }
  if (prefix[2] != 0x10u) {
    reterr = kPglRetNotYetSupported;
    goto cleanup;
  }
  variant_ct = LoadLeU32(prefix + 3, 4);
  sample_ct = LoadLeU32(prefix + 7, 4);
  control = prefix[11];
  if (variant_ct == 0u || sample_ct == 0u) {
    reterr = kPglRetMalformedInput;
    goto cleanup;
  }
  low = control & 15u;
  if (low >= 8u) {
    reterr = kPglRetNotYetSupported;
    goto cleanup;
  }
  vrec_len_byte_ct = (low & 3u) + 1u;
  vrtype_8bit = (low >> 2) & 1u;
  allele_byte_ct = (control >> 4) & 3u;
  nonref_mode = (control >> 6) & 3u;
  block_ct = DivUpU32(variant_ct, kPglVblockSize);
  pg->variant_ct = variant_ct;
  pg->sample_ct = sample_ct;
  pg->vrtypes = (unsigned char*)malloc(variant_ct);
  pg->rec_lens = (uint32_t*)malloc((size_t)variant_ct * sizeof(uint32_t));
  pg->rec_fpos = (uint64_t*)malloc((size_t)variant_ct * sizeof(uint64_t));
  pg->anchor_nword = NypCtToWordCt(sample_ct);
  pg->anchor_genovec = (uintptr_t*)calloc(pg->anchor_nword, sizeof(uintptr_t));
  block_off_buf = (unsigned char*)malloc((size_t)block_ct * 8u);
  if (unlikely(!pg->vrtypes || !pg->rec_lens || !pg->rec_fpos || !pg->anchor_genovec || !block_off_buf)) {
    reterr = kPglRetNomem;
    goto cleanup;
  }
  reterr = ReadExact(pg->fp, block_off_buf, (uint64_t)block_ct * 8u);
  if (unlikely(reterr)) {
    goto cleanup;
  }
  header_pos = 12u + (uint64_t)block_ct * 8u;
  vidx = 0;
  for (block_idx = 0; block_idx < block_ct; block_idx++) {
    uint32_t block_n = variant_ct - vidx;
    uint32_t vr_bytes;
    uint32_t len_bytes;
    uint32_t local;
    unsigned char* vr_buf;
    unsigned char* len_buf;
    uint64_t block_pos;
    uint64_t cursor;
    if (block_n > kPglVblockSize) {
      block_n = kPglVblockSize;
    }
    vr_bytes = vrtype_8bit ? block_n : DivUpU32(block_n, 2u);
    len_bytes = block_n * vrec_len_byte_ct;
    vr_buf = (unsigned char*)malloc(vr_bytes ? vr_bytes : 1u);
    len_buf = (unsigned char*)malloc(len_bytes ? len_bytes : 1u);
    if (unlikely(!vr_buf || !len_buf)) {
      free(vr_buf);
      free(len_buf);
      reterr = kPglRetNomem;
      goto cleanup;
    }
    reterr = SeekTo(pg->fp, header_pos);
    if (unlikely(!reterr)) {
      reterr = ReadExact(pg->fp, vr_buf, vr_bytes);
    }
    if (unlikely(!reterr)) {
      reterr = ReadExact(pg->fp, len_buf, len_bytes);
    }
    if (unlikely(reterr)) {
      free(vr_buf);
      free(len_buf);
      goto cleanup;
    }
    header_pos += vr_bytes;
    header_pos += len_bytes;
    if (allele_byte_ct) {
      header_pos += (uint64_t)block_n * allele_byte_ct;
    }
    if (nonref_mode == 3u) {
      header_pos += DivUpU32(block_n, 8u);
    }
    for (local = 0; local < block_n; local++) {
      if (vrtype_8bit) {
        pg->vrtypes[vidx + local] = vr_buf[local];
      } else {
        unsigned char packed = vr_buf[local / 2u];
        pg->vrtypes[vidx + local] = (local & 1u) ? (unsigned char)(packed >> 4) : (unsigned char)(packed & 15u);
      }
      pg->rec_lens[vidx + local] = LoadLeU32(len_buf + (uintptr_t)local * vrec_len_byte_ct, vrec_len_byte_ct);
      if (pg->rec_lens[vidx + local] == 0u) {
        free(vr_buf);
        free(len_buf);
        reterr = kPglRetMalformedInput;
        goto cleanup;
      }
    }
    free(vr_buf);
    free(len_buf);
    block_pos = LoadLeU64(block_off_buf + (uintptr_t)block_idx * 8u);
    cursor = block_pos;
    for (local = 0; local < block_n; local++) {
      pg->rec_fpos[vidx + local] = cursor;
      cursor += pg->rec_lens[vidx + local];
    }
    vidx += block_n;
  }
  cleanup:
  free(block_off_buf);
  if (reterr) {
    PgenClose(pg);
  }
  return reterr;
}

PglErr PgenLdDiffCount(PgenFile* pg, uint32_t vidx, uint32_t* count) {
  unsigned char buf[8];
  uint32_t rec_len;
  uint32_t nread;
  uint32_t consumed = 0;
  PglErr reterr;
  if (vidx >= pg->variant_ct || !PgenVrtypeLd(pg->vrtypes[vidx])) {
    return kPglRetMalformedInput;
  }
  rec_len = pg->rec_lens[vidx];
  nread = rec_len < 8u ? rec_len : 8u;
  if (nread == 0u) {
    return kPglRetMalformedInput;
  }
  reterr = SeekTo(pg->fp, pg->rec_fpos[vidx]);
  if (unlikely(reterr)) {
    return reterr;
  }
  reterr = ReadExact(pg->fp, buf, nread);
  if (unlikely(reterr)) {
    return reterr;
  }
  if (!PgenReadDifflistCount(buf, buf + nread, count, &consumed)) {
    return kPglRetMalformedInput;
  }
  return kPglRetSuccess;
}

PglErr PgenGetP(PgenFile* pg, uint32_t vidx,
                uintptr_t* genovec, uintptr_t* phasepresent, uintptr_t* phaseinfo,
                uint32_t* phasepresent_ct) {
  const unsigned char* start;
  const unsigned char* end;
  const unsigned char* ptr;
  PglErr reterr;
  uint32_t vrtype;
  if (vidx >= pg->variant_ct) {
    return kPglRetMalformedInput;
  }
  vrtype = pg->vrtypes[vidx];
  if (PgenVrtypeMulti(vrtype)) {
    return kPglRetNotYetSupported;
  }
  if (PgenVrtypeLd(vrtype)) {
    uint32_t anchor_vidx = 0;
    reterr = FindLdAnchor(pg, vidx, &anchor_vidx);
    if (unlikely(reterr)) {
      return reterr;
    }
    reterr = DecodeAnchor(pg, anchor_vidx);
    if (unlikely(reterr)) {
      return reterr;
    }
    reterr = LoadRecord(pg, vidx, &start, &end);
    if (unlikely(reterr)) {
      return reterr;
    }
    memcpy(genovec, pg->anchor_genovec, (size_t)pg->anchor_nword * sizeof(uintptr_t));
    ptr = start;
    reterr = ApplyDifflist(end, &ptr, pg->sample_ct, genovec);
    if (unlikely(reterr)) {
      return reterr;
    }
    if (PgenVrtypeMain(vrtype) == PGEN_VR_LD_INV) {
      GenovecInvertUnsafe(pg->sample_ct, genovec);
    }
    return DecodePhase(vrtype, end, &ptr, pg->sample_ct, genovec, phasepresent, phaseinfo, phasepresent_ct);
  }
  reterr = LoadRecord(pg, vidx, &start, &end);
  if (unlikely(reterr)) {
    return reterr;
  }
  ptr = start;
  reterr = DecodeNonLdMain(vrtype, end, &ptr, pg->sample_ct, genovec);
  if (unlikely(reterr)) {
    return reterr;
  }
  memcpy(pg->anchor_genovec, genovec, (size_t)pg->anchor_nword * sizeof(uintptr_t));
  pg->anchor_vidx = vidx;
  return DecodePhase(vrtype, end, &ptr, pg->sample_ct, genovec, phasepresent, phaseinfo, phasepresent_ct);
}
