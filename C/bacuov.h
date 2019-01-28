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
#define DEGREE_OF_CIRCULANCY 11
#endif

#define SEC_LVL 128
#define BACUOV_PARAM_O 11
#define BACUOV_PARAM_V 101
#define BACUOV_PARAM_M (BACUOV_PARAM_O * DEGREE_OF_CIRCULANCY)

typedef struct
{
    unsigned char seed[SEC_LVL/4];
} bacuov_secret_key;

typedef struct
{
    gfpcirc_element forms_P[(BACUOV_PARAM_O+BAC_UOV_PARAM_V)*(BACUOV_PARAM_O+BAC_UOV_PARAM_V) * BACUOV_PARAM_M];
} bacuov_public_key;

typedef struct
{
    gfpcirc_element vector[(BACUOV_PARAM_O+BAC_UOV_PARAM_V)];
} bacuov_signature;

int bacuov_keygen( bacuov_secret_key * sk, bacuov_public_key * pk, unsigned char * randomness );
int bacuov_sign( bacuov_signature * sig, bacuov_secret_key sk, unsigned char * msg, int msg_len );
int bacuov_verify( bacuov_public_key pk, unsigned char * msg, int msg_len, bacuov_signature * sig );

#endif

