#ifndef API_H
#define API_H

#include "bacuov.h"

#define CRYPTO_SECRETKEYBYTES (SECURITY_LEVEL/4 + ( BACUOV_PARAM_M*(BACUOV_PARAM_V*(BACUOV_PARAM_V+1)/2)*DEGREE_OF_CIRCULANCY  +  BACUOV_PARAM_M*(BACUOV_PARAM_V*BACUOV_PARAM_O)*DEGREE_OF_CIRCULANCY  +  BACUOV_PARAM_V*BACUOV_PARAM_O * DEGREE_OF_CIRCULANCY ))
#define CRYPTO_PUBLICKEYBYTES ( SECURITY_LEVEL/4 + bacuov_pack_length((BACUOV_PARAM_O*(BACUOV_PARAM_O+1)/2) * DEGREE_OF_CIRCULANCY * BACUOV_PARAM_M) )
#define CRYPTO_BYTES ( bacuov_pack_length((BACUOV_PARAM_O + BACUOV_PARAM_V) * DEGREE_OF_CIRCULANCY * EXTENSION_DEGREE) )

#define CRYPTO_ALGNAME "BACUOV"

int crypto_sign_keypair(unsigned char *pk, unsigned char *sk);
int crypto_sign(unsigned char *sm, unsigned long long *smlen, const unsigned char *m, unsigned long long mlen, const unsigned char *sk);
int crypto_sign_open(unsigned char *m, unsigned long long *mlen, const unsigned char *sm, unsigned long long smlen, const unsigned char *pk);

#endif

