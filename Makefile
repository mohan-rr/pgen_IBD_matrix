CC ?= gcc
WARN = -std=c99 -O2 -Wall -Wextra -Wshadow -Wformat-security
CPPFLAGS = -Isrc

SRCS = src/main.c src/ibd_lsh.c src/ibd_extend.c src/pgen_read.c
LIBSRCS = src/ibd_lsh.c src/ibd_extend.c src/pgen_read.c

.PHONY: all test clean

all: pgen_ibd

pgen_ibd: $(SRCS) src/pgl_base.h src/pgen_ldmeta.h src/pgen_read.h src/ibd_internal.h
	$(CC) $(WARN) $(CPPFLAGS) -o $@ $(SRCS)

test_ibd: src/test_ibd.c $(LIBSRCS) src/pgl_base.h src/pgen_ldmeta.h src/pgen_read.h src/ibd_internal.h
	$(CC) $(WARN) $(CPPFLAGS) -o $@ src/test_ibd.c $(LIBSRCS)

test: pgen_ibd test_ibd
	./test_ibd ./pgen_ibd

clean:
	rm -f pgen_ibd test_ibd
