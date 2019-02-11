CC = gcc
CFLAGS = -I.
DEPS = bacuov.c gfp.c gfpcircm.h gfpm.c Keccak-readable-and-compact-c89.c bacuov.h gfpcirc.c gfpe.c gfpm.h test_algebra.c csprng.c gfpcirc.h gfpe.h gfpx.c test_bacuov.c csprng.h gfpcircm.c gfp.h gfpx.h gfpem.c gfpem.h
ALGEBRA = gfp.c gfpcircm.h gfpm.c gfpcirc.c gfpe.c gfpm.h gfpcirc.h gfpe.h gfpx.c gfpcircm.c gfp.h gfpx.h gfpem.c gfpem.h
KECCAK = Keccak-readable-and-compact-c89.c

algebra: $(DEPS)
	$(CC) -o test_algebra test_algebra.c $(ALGEBRA) $(KECCAK) -DGF_PRIME_MODULUS=7 -DEXTENSION_DEGREE=5 -DDEFINING_POLYNOMIAL=0x0102010501 -DDEGREE_OF_CIRCULANCY=11

bacuov: $(DEPS)
	$(CC) -o test_bacuov test_bacuov.c $(ALGEBRA) $(KECCAK) -DGF_PRIME_MODULUS=7 -DEXTENSION_DEGREE=5 -DDEGREE_OF_CIRCULANCY=11

clean: 
	rm -f test_algebra
	rm -f test_bacuov


