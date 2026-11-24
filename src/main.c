#include "ibd_internal.h"

static void PrintUsage(FILE* fp) {
  fprintf(fp,
          "Usage: pgen_ibd --pgen FILE --pvar FILE --psam FILE --map FILE --out PREFIX\n"
          "           [--min-length CM] [--maf-seed F] [--mac-extend N]\n"
          "           [--perm-ct N] [--band-ct N] [--ld-skip N]\n"
          "\n"
          "Find IBD from a phased biallelic pgen with iLASH-style MinHash/LSH.\n"
          "Windows are 3.0 cM, stepped by 1.5 cM. Tiny LD diffs are skipped while\n"
          "seeding and visited again during exact extension.\n"
          "\n"
          "Defaults: --min-length 2 --maf-seed 0.01 --mac-extend 2\n"
          "          --perm-ct 40 --band-ct 10 --ld-skip sample_ct/100\n"
          "The skip threshold used is max(1, --ld-skip).\n"
          "\n"
          "Writes PREFIX.ibd (sparse matrix) and PREFIX.segments.\n"
          "The map is PLINK format: chromosome variant_id cM bp.\n");
}

static const char* TakeVal(int* idx, int argc, char** argv) {
  const char* arg = argv[*idx];
  const char* eq = strchr(arg, '=');
  if (eq) {
    return eq + 1;
  }
  if (*idx + 1 >= argc) {
    return 0;
  }
  *idx += 1;
  return argv[*idx];
}

static int ParseDoubleArg(const char* text, double* out) {
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

static int FlagIs(const char* arg, const char* name) {
  size_t nlen = strlen(name);
  if (strncmp(arg, name, nlen) != 0) {
    return 0;
  }
  return arg[nlen] == 0 || arg[nlen] == '=';
}

int main(int argc, char** argv) {
  IbdConfig cfg;
  int idx;
  PglErr reterr;
  IbdConfigDefaults(&cfg);
  if (argc <= 1) {
    PrintUsage(stderr);
    return kPglRetInvalidCmdline;
  }
  for (idx = 1; idx < argc; idx++) {
    const char* val = 0;
    if (strcmp(argv[idx], "--help") == 0 || strcmp(argv[idx], "-h") == 0) {
      PrintUsage(stdout);
      return kPglRetSuccess;
    }
    if (FlagIs(argv[idx], "--pgen")) {
      val = TakeVal(&idx, argc, argv);
      cfg.pgen_path = val;
    } else if (FlagIs(argv[idx], "--pvar")) {
      val = TakeVal(&idx, argc, argv);
      cfg.pvar_path = val;
    } else if (FlagIs(argv[idx], "--psam")) {
      val = TakeVal(&idx, argc, argv);
      cfg.psam_path = val;
    } else if (FlagIs(argv[idx], "--map")) {
      val = TakeVal(&idx, argc, argv);
      cfg.map_path = val;
    } else if (FlagIs(argv[idx], "--out")) {
      val = TakeVal(&idx, argc, argv);
      cfg.out_prefix = val;
    } else if (FlagIs(argv[idx], "--min-length")) {
      val = TakeVal(&idx, argc, argv);
      if (!val || !ParseDoubleArg(val, &cfg.min_length_cm) || cfg.min_length_cm <= 0.0) {
        fprintf(stderr, "pgen_ibd: --min-length needs a positive number\n");
        return kPglRetInvalidCmdline;
      }
    } else if (FlagIs(argv[idx], "--maf-seed")) {
      val = TakeVal(&idx, argc, argv);
      if (!val || !ParseDoubleArg(val, &cfg.maf_seed) || cfg.maf_seed < 0.0 || cfg.maf_seed > 0.5) {
        fprintf(stderr, "pgen_ibd: --maf-seed needs a number in [0, 0.5]\n");
        return kPglRetInvalidCmdline;
      }
    } else if (FlagIs(argv[idx], "--mac-extend")) {
      val = TakeVal(&idx, argc, argv);
      if (!val) {
        fprintf(stderr, "pgen_ibd: --mac-extend needs an integer\n");
        return kPglRetInvalidCmdline;
      }
      cfg.mac_extend = (uint32_t)strtoul(val, 0, 10);
    } else if (FlagIs(argv[idx], "--perm-ct")) {
      val = TakeVal(&idx, argc, argv);
      if (!val) {
        fprintf(stderr, "pgen_ibd: --perm-ct needs an integer\n");
        return kPglRetInvalidCmdline;
      }
      cfg.perm_ct = (uint32_t)strtoul(val, 0, 10);
    } else if (FlagIs(argv[idx], "--band-ct")) {
      val = TakeVal(&idx, argc, argv);
      if (!val) {
        fprintf(stderr, "pgen_ibd: --band-ct needs an integer\n");
        return kPglRetInvalidCmdline;
      }
      cfg.band_ct = (uint32_t)strtoul(val, 0, 10);
    } else if (FlagIs(argv[idx], "--ld-skip")) {
      val = TakeVal(&idx, argc, argv);
      if (!val) {
        fprintf(stderr, "pgen_ibd: --ld-skip needs an integer\n");
        return kPglRetInvalidCmdline;
      }
      cfg.ld_skip = (int32_t)strtol(val, 0, 10);
    } else {
      fprintf(stderr, "pgen_ibd: unknown argument %s\n", argv[idx]);
      PrintUsage(stderr);
      return kPglRetInvalidCmdline;
    }
    if (!val) {
      fprintf(stderr, "pgen_ibd: missing value for %s\n", argv[idx]);
      return kPglRetInvalidCmdline;
    }
  }
  reterr = IbdRun(&cfg, 0);
  return (int)reterr;
}
