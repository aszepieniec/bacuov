#include "bacuov.h"

void FIPS202_SHAKE256(const unsigned char *input, unsigned int inputByteLen, unsigned char *output, int outputByteLen);
void FIPS202_SHA3_256(const unsigned char *input, unsigned int inputByteLen, unsigned char *output);

/**
 * bacuov_generate_column_S
 * Generate a single indicated column of the secret matrix S, given
 * the seed for S.
 */
void bacuov_generate_column_S( gfpcirc_element * col, unsigned char * seed_S, int col_idx )
{
    unsigned char randomness[SEC_LVL/4 + sizeof(int) + sizeof(int)];
    int i, k;
    gfpcirc_element elm;
    unsigned char buffer[BYTES_PER_RINGELEMENT * BACUOV_PARAM_V];
    
    // copy over randomness
    for( i = 0 ; i < SEC_LVL/4 + sizeof(int) ; ++i )
    {
        randomness[i] = seed[i];
    };
    // append indicators
    type = 0 // S
    for( i = 0 ; i < sizeof(unsigned int) ; ++i )
    {
        randomness[SEC_LVL/4+i] = (unsigned char*)(&type)[i];
        randomness[SEC_LVL/4 + sizeof(int) + i] = (unsigned char*)(&col_idx)[i];
    }

    // expand seed to buffer of mucho pseudorandomness
    FIPS202_SHAKE256(randomness, SEC_LVL/4 + sizeof(int)*2, buffer, BYTES_PER_RINGELEMENT * BACUOV_PARAM_O);

    // use pseudorandomness buffer to sample circulant ring elements
    k = 0;
    for( i = 0 ; i < BACUOV_PARAM_V ; ++i )
    {
        gfpcirc_sample(&col[i], buffer + (k++)*BYTES_PER_RINGELEMENT);
    }
}

/**
 * bacuov_generate_vinegar_coefficients
 * Generate the top-left VxV block of circulant ring elements of the
 * public quadratic form Pi.
 */
void bacuov_generate_vinegar_coefficients( gfpcirc_element * Pi0V0V, unsigned char * seed, int index )
{
    unsigned char randomness[SEC_LVL/4 + sizeof(int) + sizeof(int)];
    int i, j, k;
    gfpcirc_element elm;
    unsigned char buffer[BYTES_PER_RINGELEMENT * (BACUOV_PARAM_V * (BACOUV_PARAM_V + 1) / 2)];
    
    // copy over randomness
    for( i = 0 ; i < SEC_LVL/4 + sizeof(int) ; ++i )
    {
        randomness[i] = seed[i];
    };
    // append indicator
    // append indicators
    type = 1 // vinegars
    for( i = 0 ; i < sizeof(unsigned int) ; ++i )
    {
        randomness[SEC_LVL/4+i] = (unsigned char*)(&type)[i];
        randomness[SEC_LVL/4 + sizeof(int) + i] = (unsigned char*)(&index)[i];
    }
    randomness[SEC_LVL/4] = (unsigned char)-1;

    // expand seed to buffer of mucho pseudorandomness
    FIPS202_SHAKE256(randomness, SEC_LVL/4 + sizeof(int)*2, buffer, BYTES_PER_RINGELEMENT * (BACUOV_PARAM_V * (BACUOV_PARAM_V + 1) / 2));

    // use pseudorandomness buffer to sample circulant ring elements
    k = 0;
    for( i = 0 ; i < BACUOV_PARAM_V ; ++i )
    {
        gfpcirc_sample(&Pi0V0V[i*BACUOV_PARAM_V+i], buffer + (k++)*BYTES_PER_RINGELEMENT);
        for( j = i+1 ; j < BACUOV_PARAM_V ; ++j )
        {
            gfpcirc_sample(&elm, buffer + (k++)*BYTES_PER_RINGELEMENT);
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
    unsigned char randomness[SEC_LVL/4 + sizeof(int) + sizeof(int)];
    int type;
    int i, j, k;
    gfpcirc_element elm;
    unsigned char buffer[BYTES_PER_RINGELEMENT * (BACUOV_PARAM_V * BACOUV_PARAM_O)];
    
    // copy over randomness
    for( i = 0 ; i < SEC_LVL/4 ; ++i )
    {
        randomness[i] = seed[i];
    };
    // append indicators
    type = 2 // linears
    for( i = 0 ; i < sizeof(unsigned int) ; ++i )
    {
        randomness[SEC_LVL/4+i] = (unsigned char*)(&type)[i];
        randomness[SEC_LVL/4 + sizeof(int) + i] = (unsigned char*)(&index)[i];
    }

    // expand seed to buffer of mucho pseudorandomness
    FIPS202_SHAKE256(randomness, SEC_LVL/4 + sizeof(int)*2, buffer, BYTES_PER_RINGELEMENT * (BACUOV_PARAM_V * BACUOV_PARAM_O));

    // use pseudorandomness buffer to sample circulant ring elements
    k = 0;
    for( i = 0 ; i < BACUOV_PARAM_V ; ++i )
    {
        for( j = 0 ; j < BACUOV_PARAM_O ; ++j )
        {
            gfpcirc_sample(&Pi0VVN[i*BACUOV_PARAM_O+j], buffer + (k++)*BYTES_PER_RINGELEMENT);
        }
    }
}

void bacuov_keygen( bacuov_secret_key * sk, bacuov_public_key * pk, unsigned char * randomness )
{
    int i, j, k;
    unsigned char buffer[SEC_LVL/4 + 1];
    unsigned char outputbuf[32];
    gfp_element Pi0V0V[BACUOV_PARAM_V * BACUOV_PARAM_V];
    gfp_element Pi0VVN[BACUOV_PARAM_V * BACUOV_PARAM_O];

    // initialize seed buffer
    for( i = 0 ; i < SEC_LVL/4 ; ++i )
    {
        buffer[i] = randomness[i];
    }

    // generate seed for sk
    buffer[SEC_LVL/4] = (unsigned char)-1;
    SHA3_256(buffer, SEC_LVL/4+1, outputbuf);
    for( i = 0 ; i < SEC_LVL/4 ; ++i )
    {
        sk->seed[i] = outputbuf[i];
    }

    // generate seed for pk
    buffer[SEC_LVL/4] = 1;
    SHA3_256(buffer, SEC_LVL/4+1, outputbuf);
    for( i = 0 ; i < SEC_LVL/4 ; ++i )
    {
        pk->seed[i] = outputbuf[i];
    };

    // for all m polynomials
    for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
    {
        // populate top-left and top-right blocks of Pi
        bacuov_generate_vinegar_coefficients(Pi0V0V, pk->seed, i);
        bacuov_generate_linear_coefficients(Pi0VVN, pk->seed, i);

        //
    }
}

