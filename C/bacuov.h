#ifndef BACUOV_H
#define BACUOV_H

#ifndef GF_PRIME_MODULUS
#define GF_PRIME_MODULUS 7
#endif

// #ifndef GFPE_EXTENSION_DEGREE
// #define GFPE_DEFINING_POLYNOMIAL 0x23ae734f
// #define GFPE_EXTENSION_DEGREE 5
// #endif

#ifndef DEGREE_OF_CIRCULANCY
#define DEGREE_OF_CIRCULANCY 3
#endif

#ifndef SECURITY_LEVEL
#define SECURITY_LEVEL 256
#endif

#ifndef BACUOV_PARAM_O
#define BACUOV_PARAM_O 2
#endif

#ifndef BACUOV_PARAM_V
#define BACUOV_PARAM_V 7
#endif

#ifndef BACUOV_PARAM_M
#define BACUOV_PARAM_M (BACUOV_PARAM_O * DEGREE_OF_CIRCULANCY)
#endif

#include "gfp.h"
#include "gfpcirc.h"
#include "gfpcircm.h"

typedef struct
{
    gfpcircmatrix FFvv[BACUOV_PARAM_M];
    gfpcircmatrix FFvo[BACUOV_PARAM_M];
    gfpcircmatrix S;
} bacuov_secret_key;

typedef struct
{
    gfpcircmatrix PPoo[BACUOV_PARAM_M];
    unsigned char seed_PCT[SECURITY_LEVEL/4];
} bacuov_public_key;

typedef struct
{
    gfpcirc_element vector[(BACUOV_PARAM_O+BACUOV_PARAM_V)];
} bacuov_signature;


void bacuov_secret_key_init( bacuov_secret_key * sk );
void bacuov_secret_key_destroy( bacuov_secret_key sk );
void bacuov_public_key_init( bacuov_public_key * sk );
void bacuov_public_key_destroy( bacuov_public_key sk );
void bacuov_keygen( bacuov_secret_key * sk, bacuov_public_key * pk, unsigned char * randomness );
void bacuov_sign( bacuov_signature * sig, bacuov_secret_key sk, unsigned char * msg, int msg_len );
int bacuov_verify( bacuov_public_key pk, unsigned char * msg, int msg_len, bacuov_signature * sig );

#endif

