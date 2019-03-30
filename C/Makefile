CC = gcc
CFLAGS = -I.
DEPS = bacuov.c gfp.c gfpcircm.h gfpm.c Keccak-readable-and-compact-c89.c bacuov.h gfpcirc.c gfpe.c gfpm.h test_algebra.c csprng.c gfpcirc.h gfpe.h gfpx.c test_bacuov.c csprng.h gfpcircm.c gfp.h gfpx.h gfpem.c gfpem.h api.h test_pack.c benchmark.c api.c
ALGEBRA = gfp.c gfpcircm.h gfpm.c gfpcirc.c gfpe.c gfpm.h gfpcirc.h gfpe.h gfpx.c gfpcircm.c gfp.h gfpx.h gfpem.c gfpem.h
KECCAK = Keccak-readable-and-compact-c89.c

algebra: $(DEPS)
	$(CC) -o test_algebra test_algebra.c $(ALGEBRA) $(KECCAK) -DGF_PRIME_MODULUS=7 -DEXTENSION_DEGREE=5 -DDEFINING_POLYNOMIAL=\"\x03\x01\x00\x00\x00\" -DDEGREE_OF_CIRCULANCY=11

bacuov: $(DEPS)
	$(CC) -o test_bacuov test_bacuov.c bacuov.c $(ALGEBRA) $(KECCAK) -DGF_PRIME_MODULUS=7 -DEXTENSION_DEGREE=5 -DDEFINING_POLYNOMIAL=\"\x03\x01\x00\x00\x00\" -DDEGREE_OF_CIRCULANCY=11 -DDEBUG

pack: $(DEPS)
	$(CC) -o test_pack test_pack.c bacuov.c $(ALGEBRA) $(KECCAK) -DGF_PRIME_MODULUS=7 -DGFP_NUMBITS=3 -DEXTENSION_DEGREE=5 -DDEFINING_POLYNOMIAL=\"\x03\x01\x00\x00\x00\" -DDEGREE_OF_CIRCULANCY=11 -DDEBUG

bench: $(DEPS)
	$(CC) -o benchmark benchmark.c bacuov.c api.c $(ALGEBRA) $(KECCAK) -DGF_PRIME_MODULUS=7 -DGFP_NUMBITS=3 -DEXTENSION_DEGREE=5 -DDEFINING_POLYNOMIAL=\"\x03\x01\x00\x00\x00\" -DDEGREE_OF_CIRCULANCY=11 -O3

benchII3: $(DEPS)
	$(CC) -o benchmarkII3 benchmark.c bacuov.c api.c $(ALGEBRA) $(KECCAK) -O3 -DGF_PRIME_MODULUS=3 -DGFP_NUMBITS=2 -DEXTENSION_DEGREE=12 -DDEFINING_POLYNOMIAL="\"\x02\x00\x01\x00\x00\x00\x00\x00\x00\x00\x00\x00\"" -DDEGREE_OF_CIRCULANCY=7 -DBACUOV_PARAM_O=8 -DBACUOV_PARAM_V=56

benchII7: $(DEPS)
	$(CC) -o benchmarkII7 benchmark.c bacuov.c api.c $(ALGEBRA) $(KECCAK) -O3 -DGF_PRIME_MODULUS=7 -DGFP_NUMBITS=3 -DEXTENSION_DEGREE=5 -DDEFINING_POLYNOMIAL="\"\x03\x01\x00\x00\x00\"" -DDEGREE_OF_CIRCULANCY=13 -DBACUOV_PARAM_O=5 -DBACUOV_PARAM_V=60

benchIV3: $(DEPS)
	$(CC) -o benchmarkIV3 benchmark.c bacuov.c api.c $(ALGEBRA) $(KECCAK) -O3 -DGF_PRIME_MODULUS=3 -DGFP_NUMBITS=2 -DEXTENSION_DEGREE=16 -DDEFINING_POLYNOMIAL="\"\x01\x00\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\"" -DDEGREE_OF_CIRCULANCY=7 -DBACUOV_PARAM_O=11 -DBACUOV_PARAM_V=84

benchIV7: $(DEPS)
	$(CC) -o benchmarkIV7 benchmark.c bacuov.c api.c $(ALGEBRA) $(KECCAK) -O3 -DGF_PRIME_MODULUS=7 -DGFP_NUMBITS=3 -DEXTENSION_DEGREE=16 -DDEFINING_POLYNOMIAL="\"\x03\x02\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\"" -DDEGREE_OF_CIRCULANCY=11 -DBACUOV_PARAM_O=7 -DBACUOV_PARAM_V=76

benchV3: $(DEPS)
	$(CC) -o benchmarkV3 benchmark.c bacuov.c api.c $(ALGEBRA) $(KECCAK) -O3 -DGF_PRIME_MODULUS=3 -DGFP_NUMBITS=2 -DEXTENSION_DEGREE=16 -DDEFINING_POLYNOMIAL="\"\x01\x00\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\"" -DDEGREE_OF_CIRCULANCY=7 -DBACUOV_PARAM_O=14 -DBACUOV_PARAM_V=104

benchV7: $(DEPS)
	$(CC) -o benchmarkV7 benchmark.c bacuov.c api.c $(ALGEBRA) $(KECCAK) -O3 -DGF_PRIME_MODULUS=7 -DGFP_NUMBITS=3 -DEXTENSION_DEGREE=8 -DDEFINING_POLYNOMIAL="\"\x03\x01\x00\x00\x00\x00\x00\x00\"" -DDEGREE_OF_CIRCULANCY=11 -DBACUOV_PARAM_O=9 -DBACUOV_PARAM_V=97

clean: 
	rm -f test_algebra
	rm -f test_bacuov
	rm -f test_pack


