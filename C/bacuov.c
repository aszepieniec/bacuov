#include "bacuov.h"

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
void bacuov_generate_S( gfpcirc_element * data, unsigned char * seed_S )
{
    unsigned char randomness[SECURITY_LEVEL/4 + sizeof(int)];
    int i, j, k;
    int type;
    gfpcirc_element elm;
    unsigned char buffer[sizeof(gfpcirc_element) * BACUOV_PARAM_V*BACUOV_PARAM_O];
    
    // copy over randomness
    for( i = 0 ; i < SECURITY_LEVEL/4 + sizeof(int) ; ++i )
    {
        randomness[i] = seed_S[i];
    };
    // append indicator
    type = 0; // S
    for( i = 0 ; i < sizeof(unsigned int) ; ++i )
    {
        randomness[SECURITY_LEVEL/4+i] = (unsigned char)(type >> (8*i));
    }

    // expand seed to buffer of mucho pseudorandomness
    FIPS202_SHAKE256(randomness, SECURITY_LEVEL/4 + sizeof(int), buffer, sizeof(gfpcirc_element) * BACUOV_PARAM_V*BACUOV_PARAM_O);

    // use pseudorandomness buffer to sample circulant ring elements
    k = 0;
    for( i = 0 ; i < BACUOV_PARAM_V ; ++i )
    {
        for( j = 0 ; j < BACUOV_PARAM_O ; ++j )
        {
            gfpcirc_sample(&data[i*BACUOV_PARAM_O + j], buffer + (k++)*sizeof(gfpcirc_element));
        }
    }
}

/**
 * bacuov_generate_vinegar_coefficients
 * Generate the top-left VxV block of circulant ring elements of the
 * public quadratic form Pi.
 */
void bacuov_generate_vinegar_coefficients( gfpcirc_element * Pi0V0V, unsigned char * seed, int index )
{
    unsigned char randomness[SECURITY_LEVEL/4 + sizeof(int) + sizeof(int)];
    int i, j, k;
    int type;
    gfpcirc_element elm;
    unsigned char buffer[sizeof(gfpcirc_element) * (BACUOV_PARAM_V * (BACUOV_PARAM_V + 1) / 2)];
    
    // copy over randomness
    for( i = 0 ; i < SECURITY_LEVEL/4 + sizeof(int) ; ++i )
    {
        randomness[i] = seed[i];
    };
    // append indicator
    // append indicators
    type = 1; // vinegars
    for( i = 0 ; i < sizeof(unsigned int) ; ++i )
    {
        randomness[SECURITY_LEVEL/4+i] = (unsigned char)(type >> (8*i));
        randomness[SECURITY_LEVEL/4 + sizeof(int) + i] = (unsigned char)(index >> (8*i));
    }
    randomness[SECURITY_LEVEL/4] = (unsigned char)-1;

    // expand seed to buffer of mucho pseudorandomness
    FIPS202_SHAKE256(randomness, SECURITY_LEVEL/4 + sizeof(int)*2, buffer, sizeof(gfpcirc_element) * (BACUOV_PARAM_V * (BACUOV_PARAM_V + 1) / 2));

    // use pseudorandomness buffer to sample circulant ring elements
    k = 0;
    for( i = 0 ; i < BACUOV_PARAM_V ; ++i )
    {
        gfpcirc_sample(&Pi0V0V[i*BACUOV_PARAM_V+i], buffer + (k++)*sizeof(gfpcirc_element));
        for( j = i+1 ; j < BACUOV_PARAM_V ; ++j )
        {
            gfpcirc_sample(&elm, buffer + (k++)*sizeof(gfpcirc_element));
            // divide elm by 2? not necessary
            gfpcirc_copy(&Pi0V0V[i*BACUOV_PARAM_V+j], elm);
            gfpcirc_copy(&Pi0V0V[j*BACUOV_PARAM_V+i], elm);
        }
    }
}

/**
 * bacuov_generate_linear_coefficients
 * Generate the top-right vinegar-oil part of the secret quadratic
 * form Pi.
 */
void bacuov_generate_linear_coefficients( gfpcirc_element * Pi0VVN, unsigned char * seed, int index )
{
    unsigned char randomness[SECURITY_LEVEL/4 + sizeof(int) + sizeof(int)];
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
    for( i = 0 ; i < sizeof(unsigned int) ; ++i )
    {
        randomness[SECURITY_LEVEL/4+i] = (unsigned char)(type >> (8*i));
        randomness[SECURITY_LEVEL/4 + sizeof(int) + i] = (unsigned char)(index >> (8*i));
    }

    // expand seed to buffer of mucho pseudorandomness
    FIPS202_SHAKE256(randomness, SECURITY_LEVEL/4 + sizeof(int)*2, buffer, sizeof(gfpcirc_element) * (BACUOV_PARAM_V * BACUOV_PARAM_O));

    // use pseudorandomness buffer to sample circulant ring elements
    k = 0;
    for( i = 0 ; i < BACUOV_PARAM_V ; ++i )
    {
        for( j = 0 ; j < BACUOV_PARAM_O ; ++j )
        {
            gfpcirc_sample(&Pi0VVN[i*BACUOV_PARAM_O+j], buffer + (k++)*sizeof(gfpcirc_element));
        }
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
    gfpcirc_element Pi0V0V[BACUOV_PARAM_V * BACUOV_PARAM_V];
    gfpcirc_element Fi0V0V[BACUOV_PARAM_V * BACUOV_PARAM_V];
    gfpcirc_element Pi0VVN[BACUOV_PARAM_V * BACUOV_PARAM_O];
    gfpcirc_element Fi0VVN[BACUOV_PARAM_V * BACUOV_PARAM_O];
    gfpcirc_element S[BACUOV_PARAM_V * BACUOV_PARAM_O];
    gfpcirc_element temp_data[BACUOV_PARAM_O * BACUOV_PARAM_O];
    gfpcirc_element temp_data1[BACUOV_PARAM_O * BACUOV_PARAM_V];
    gfpcircmatrix Pi0V0V_matrix, Pi0VVN_matrix;
    gfpcircmatrix PiVNVN_matrix;
    gfpcircmatrix Fi0V0V_matrix, Fi0VVN_matrix;
    gfpcircmatrix S_matrix;
    gfpcircmatrix temp_matrix, temp_matrix1;

    // initialize matrices
    Pi0V0V_matrix.width = BACUOV_PARAM_V;
    Pi0V0V_matrix.height = BACUOV_PARAM_V;
    Pi0V0V_matrix.data = Pi0V0V;
    Pi0VVN_matrix.width = BACUOV_PARAM_V;
    Pi0VVN_matrix.height = BACUOV_PARAM_O;
    Pi0VVN_matrix.data = Pi0VVN;
    PiVNVN_matrix.width = BACUOV_PARAM_O;
    PiVNVN_matrix.height = BACUOV_PARAM_O;
    Fi0V0V_matrix.width = BACUOV_PARAM_V;
    Fi0V0V_matrix.height = BACUOV_PARAM_V;
    Fi0V0V_matrix.data = Fi0V0V;
    Fi0VVN_matrix.width = BACUOV_PARAM_V;
    Fi0VVN_matrix.height = BACUOV_PARAM_O;
    Fi0VVN_matrix.data = Fi0VVN;
    S_matrix.width = BACUOV_PARAM_O;
    S_matrix.height = BACUOV_PARAM_V;
    S_matrix.data = S;
    temp_matrix.width = BACUOV_PARAM_O;
    temp_matrix.height = BACUOV_PARAM_O;
    temp_matrix.data = temp_data;
    temp_matrix1.width = BACUOV_PARAM_V;
    temp_matrix1.height = BACUOV_PARAM_O;
    temp_matrix1.data = temp_data1;


    // expand randomness into seeds
    FIPS202_SHAKE256(randomness, SECURITY_LEVEL/4, buffer, 2*SECURITY_LEVEL/4);

    // grab seed for S and PCT
    for( i = 0 ; i < SECURITY_LEVEL/4 ; ++i )
    {
        seed_S[i] = buffer[i];
        pk->seed_PCT[i] = buffer[SECURITY_LEVEL/4 + i];
    }

    // expand seed for S into S proper
    bacuov_generate_S(S, seed_S);

    // for all m polynomials
    for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
    {
        // populate top-left and top-right blocks of Pi
        bacuov_generate_vinegar_coefficients(Pi0V0V, pk->seed_PCT, i);
        bacuov_generate_linear_coefficients(Pi0VVN, pk->seed_PCT, i);

        // copy Fi0V0V
        gfpcircm_copy(Fi0V0V_matrix, Pi0V0V_matrix);
        gfpcircm_flip(Fi0V0V_matrix);

        // compute Fi0VVN
        gfpcircm_multiply(&Fi0VVN_matrix, Pi0V0V_matrix, S_matrix);
        gfpcircm_flip(Pi0VVN_matrix);
        gfpcircm_add(Fi0VVN_matrix, Fi0VVN_matrix, Pi0VVN_matrix);
        gfpcircm_flip(Pi0VVN_matrix);

        // compute PiVNVN
        PiVNVN_matrix.data = pk->PP + sizeof(gfpcirc_element) * BACUOV_PARAM_O*BACUOV_PARAM_O * i;
        gfpcircm_transpose_multiply(&temp_matrix1, S_matrix, Fi0V0V_matrix);
        gfpcircm_multiply(&PiVNVN_matrix, temp_matrix1, S_matrix);
        gfpcircm_transpose_multiply(&temp_matrix, S_matrix, Fi0VVN_matrix);
        gfpcircm_flip(temp_matrix);
        gfpcircm_add(PiVNVN_matrix, PiVNVN_matrix, temp_matrix);
        gfpcircm_transpose(&temp_matrix);
        gfpcircm_add(PiVNVN_matrix, PiVNVN_matrix, temp_matrix);
    }
}

