#include "bacuov.h"

#include <stdio.h>

void FIPS202_SHAKE256(const unsigned char *input, unsigned int inputByteLen, unsigned char *output, int outputByteLen);
void FIPS202_SHA3_256(const unsigned char *input, unsigned int inputByteLen, unsigned char *output);

#ifndef BACUOV_PARAM_V
#error "parameter V undefined"
#endif

/**
 * bacuov_generate__S
 * Generate the secret matrix S, given
 * the seed for S.
 */
void bacuov_generate_S( gfpcircmatrix * data, unsigned char * seed_S )
{
    unsigned char randomness[SECURITY_LEVEL/4 + 4];
    int i, j, k;
    int type;
    gfpcirc_element elm;
    unsigned char buffer[GFP_NUMBYTES * DEGREE_OF_CIRCULANCY * BACUOV_PARAM_V*BACUOV_PARAM_O];

    // copy over randomness
    for( i = 0 ; i < SECURITY_LEVEL/4 + sizeof(int) ; ++i )
    {
        randomness[i] = seed_S[i];
    };
    // append indicator
    type = 0; // S
    for( i = 0 ; i < 4 ; ++i )
    {
        randomness[SECURITY_LEVEL/4+i] = (unsigned char)(type >> (8*i));
    }

    // expand seed to buffer of mucho pseudorandomness
    FIPS202_SHAKE256(randomness, SECURITY_LEVEL/4 + 4, buffer, GFP_NUMBYTES*DEGREE_OF_CIRCULANCY * BACUOV_PARAM_V*BACUOV_PARAM_O);
    printf("buffer: ");
    for( i = 0 ; i < GFP_NUMBYTES*DEGREE_OF_CIRCULANCY * BACUOV_PARAM_V * BACUOV_PARAM_O ; ++i )
        printf("%02x", buffer[i]);
    printf("\n");

    // use pseudorandomness buffer to sample circulant ring elements
    k = 0;
    for( i = 0 ; i < BACUOV_PARAM_V ; ++i )
    {
        for( j = 0 ; j < BACUOV_PARAM_O ; ++j )
        {
            gfpcirc_random(&data->data[i*BACUOV_PARAM_O + j], buffer + (k++)*GFP_NUMBYTES*DEGREE_OF_CIRCULANCY);
        }
    }

    printf("S:\n");
    printf("%i x %i matrix (should be %i x %i)\n", data->height, data->width, BACUOV_PARAM_V, BACUOV_PARAM_O);
    gfpcircm_print(*data);
    printf("done with generating S.\n"); getchar();
}

/**
 * bacuov_generate_vinegar_coefficients
 * Generate the top-left VxV block of circulant ring elements of the
 * public quadratic form Pi.
 */
void bacuov_generate_vinegar_coefficients( gfpcircmatrix * Pi0V0V, unsigned char * seed, int index )
{
    unsigned char randomness[SECURITY_LEVEL/4 + 8];
    int i, j, k;
    int type;
    gfpcirc_element elm;
    unsigned char buffer[GFP_NUMBYTES*DEGREE_OF_CIRCULANCY * (BACUOV_PARAM_V * (BACUOV_PARAM_V + 1) / 2)];
    
    // copy over randomness
    for( i = 0 ; i < SECURITY_LEVEL/4 + sizeof(int) ; ++i )
    {
        randomness[i] = seed[i];
    };
    // append indicator
    // append indicators
    type = 1; // vinegars
    for( i = 0 ; i < 4 ; ++i )
    {
        randomness[SECURITY_LEVEL/4+i] = (unsigned char)(type >> (8*i));
        randomness[SECURITY_LEVEL/4 + 4 + i] = (unsigned char)(index >> (8*i));
    }

    // expand seed to buffer of mucho pseudorandomness
    FIPS202_SHAKE256(randomness, SECURITY_LEVEL/4 + 8, buffer, GFP_NUMBYTES*DEGREE_OF_CIRCULANCY * (BACUOV_PARAM_V * (BACUOV_PARAM_V + 1) / 2));

    //printf("buffer: ");
    //for( i = 0 ; i < GFP_NUMBYTES*DEGREE_OF_CIRCULANCY * BACUOV_PARAM_V * (BACUOV_PARAM_V+1) / 2 ; ++i )
    //    printf("%02x", buffer[i]);
    //printf(" <-- got it!\n");

    // use pseudorandomness buffer to sample circulant ring elements
    k = 0;
    for( i = 0 ; i < BACUOV_PARAM_V ; ++i )
    {
        gfpcirc_random(&Pi0V0V->data[i*BACUOV_PARAM_V+i], buffer + (k++)*GFP_NUMBYTES*DEGREE_OF_CIRCULANCY);
        for( j = i+1 ; j < BACUOV_PARAM_V ; ++j )
        {
            gfpcirc_random(&elm, buffer + (k++)*GFP_NUMBYTES*DEGREE_OF_CIRCULANCY);
            // divide elm by 2? not necessary
            gfpcirc_copy(&Pi0V0V->data[i*BACUOV_PARAM_V+j], &elm);
            gfpcirc_copy(&Pi0V0V->data[j*BACUOV_PARAM_V+i], &elm);
        }
    }
}

/**
 * bacuov_generate_linear_coefficients
 * Generate the top-right vinegar-oil part of the secret quadratic
 * form Pi.
 */
void bacuov_generate_linear_coefficients( gfpcircmatrix * Pi0VVN, unsigned char * seed, int index )
{
    unsigned char randomness[SECURITY_LEVEL/4 + 8];
    int type;
    int i, j, k;
    gfpcirc_element elm;
    unsigned char buffer[sizeof(gfpcirc_element) * (BACUOV_PARAM_V * BACUOV_PARAM_O)];
    
    // copy over randomness
    for( i = 0 ; i < SECURITY_LEVEL/4 ; ++i )
    {
        randomness[i] = seed[i];
    };
    // append indicators
    type = 2; // linears
    for( i = 0 ; i < 4 ; ++i )
    {
        randomness[SECURITY_LEVEL/4+i] = (unsigned char)(type >> (8*i));
        randomness[SECURITY_LEVEL/4 + 4 + i] = (unsigned char)(index >> (8*i));
    }

    // expand seed to buffer of mucho pseudorandomness
    FIPS202_SHAKE256(randomness, SECURITY_LEVEL/4 + 8, buffer, GFP_NUMBYTES * DEGREE_OF_CIRCULANCY * (BACUOV_PARAM_V * BACUOV_PARAM_O));

    printf("buffer: ");
    for( i = 0 ; i < GFP_NUMBYTES * DEGREE_OF_CIRCULANCY * BACUOV_PARAM_V * BACUOV_PARAM_O ; ++i )
        printf("%02x", buffer[i]);
    printf("\n");

    // use pseudorandomness buffer to sample circulant ring elements
    k = 0;
    for( i = 0 ; i < BACUOV_PARAM_V ; ++i )
    {
        for( j = 0 ; j < BACUOV_PARAM_O ; ++j )
        {
            gfpcirc_random(&Pi0VVN->data[i*BACUOV_PARAM_O+j], buffer + (k++)*GFP_NUMBYTES * DEGREE_OF_CIRCULANCY);
        }
    }
}

void bacuov_secret_key_init( bacuov_secret_key * sk )
{
    int i;
    sk->S = gfpcircm_init(BACUOV_PARAM_V, BACUOV_PARAM_O);
    for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
    {
        sk->FFvv[i] = gfpcircm_init(BACUOV_PARAM_V, BACUOV_PARAM_V);
        sk->FFvo[i] = gfpcircm_init(BACUOV_PARAM_V, BACUOV_PARAM_O);
    }
}

void bacuov_secret_key_destroy( bacuov_secret_key sk )
{
    int i;
    gfpcircm_destroy(sk.S);
    for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
    {
        gfpcircm_destroy(sk.FFvv[i]);
        gfpcircm_destroy(sk.FFvo[i]);
    }
}

void bacuov_public_key_init( bacuov_public_key * pk )
{
    int i;
    for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
    {
        pk->PPoo[i] = gfpcircm_init(BACUOV_PARAM_O, BACUOV_PARAM_O);
    }
}

void bacuov_public_key_destroy( bacuov_public_key pk )
{
    int i;
    for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
    {
        gfpcircm_destroy(pk.PPoo[i]);
    }
}

/**
 * bacuov_keygen
 * Generate secrey and public keypair for BAC UOV cryptosystem.
 */
void bacuov_keygen( bacuov_secret_key * sk, bacuov_public_key * pk, unsigned char * randomness )
{
    int i, j, k;
    unsigned char buffer[SECURITY_LEVEL/4 + SECURITY_LEVEL/4];
    unsigned char seed_S[SECURITY_LEVEL/4];
    gfpcircmatrix Pi0V0V, Pi0VVN;
    gfpcircmatrix tempOV, tempOO;

    printf("inside keygen with randomness ");
    for( i = 0 ; i < SECURITY_LEVEL/4 ; ++i )
    {
        printf("%02x", randomness[i]);
    }
    printf("\n");

    // initialize matrices
    Pi0V0V = gfpcircm_init(BACUOV_PARAM_V, BACUOV_PARAM_V);
    Pi0VVN = gfpcircm_init(BACUOV_PARAM_V, BACUOV_PARAM_O);
    tempOV = gfpcircm_init(BACUOV_PARAM_O, BACUOV_PARAM_V);
    tempOO = gfpcircm_init(BACUOV_PARAM_O, BACUOV_PARAM_O);

    // expand randomness into seeds
    FIPS202_SHAKE256(randomness, SECURITY_LEVEL/4, buffer, 2*SECURITY_LEVEL/4);


    // grab seed for S and PCT
    for( i = 0 ; i < SECURITY_LEVEL/4 ; ++i )
    {
        seed_S[i] = buffer[i];
        pk->seed_PCT[i] = buffer[SECURITY_LEVEL/4 + i];
    }
    printf("PCT seed: ");
    for( i = 0 ; i < SECURITY_LEVEL / 4 ; ++i )
        printf("%02x", pk->seed_PCT[i]);
    printf("\n");

    // expand seed for S into S proper
    bacuov_generate_S(&sk->S, seed_S);

    // for all m polynomials
    for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
    {
        // populate top-left and top-right blocks of Pi
        bacuov_generate_vinegar_coefficients(&Pi0V0V, pk->seed_PCT, i);
        bacuov_generate_linear_coefficients(&Pi0VVN, pk->seed_PCT, i);


        // copy Fi0V0V
        gfpcircm_copy(sk->FFvv[i], Pi0V0V);
        gfpcircm_flip(sk->FFvv[i]);


        // compute Fi0VVN
        gfpcircm_multiply(&sk->FFvo[i], Pi0V0V, sk->S);
        gfpcircm_flip(Pi0VVN);
        gfpcircm_add(sk->FFvo[i], sk->FFvo[i], Pi0VVN);
        gfpcircm_flip(Pi0VVN);

        // compute PiVNVN
        gfpcircm_transpose_multiply(&tempOV, sk->S, sk->FFvv[i]);
        gfpcircm_multiply(&pk->PPoo[i], tempOV, sk->S);
        gfpcircm_transpose_multiply(&tempOO, sk->S, sk->FFvo[i]);
        gfpcircm_flip(tempOO);
        gfpcircm_add(pk->PPoo[i], pk->PPoo[i], tempOO);
        gfpcircm_transpose(&tempOO);
        gfpcircm_add(pk->PPoo[i], pk->PPoo[i], tempOO);
    }


    // destroy matrices
    
    gfpcircm_destroy(Pi0V0V);
    gfpcircm_destroy(Pi0VVN);
    gfpcircm_destroy(tempOV);
    gfpcircm_destroy(tempOO);

    getchar();
}

