#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bacuov.h"

void FIPS202_SHAKE256(const unsigned char *input, unsigned int inputByteLen, unsigned char *output, int outputByteLen);

int main( int argc, char ** argv )
{
    unsigned long randomness;
    unsigned char * seed;
    unsigned char key_seed[SECURITY_LEVEL/4], key_seed2[SECURITY_LEVEL/4];
    int i;
    int equals;
    int num_trials, trial_index;
    int num_successes;
    int num_failures;

    bacuov_public_key pk;
    bacuov_secret_key sk;

    if( argc != 3 || strlen(argv[2]) % 2 != 0 )
    {
        printf("usage: ./test_bacuov [num trials, eg 13] [random seed, eg d13d13deadbeef]\n");
        return 0;
    }

    /* grab randomness */

    seed = malloc(strlen(argv[2])/2);
    for( i = 0 ; i < strlen(argv[2])/2 ; ++i )
    {
        sscanf(argv[2] + 2*i, "%2hhx", &seed[i]);
    }
    SHAKE256(seed, strlen(argv[2])/2, key_seed2, SECURITY_LEVEL/4);
    free(seed);

    /* grab trial number */
    num_trials = atoi(argv[1]);

    /* run trials */
    for( trial_index = 0 ; trial_index < num_trials ; ++trial_index )
    {
        SHAKE256(key_seed2, SECURITY_LEVEL/4, key_seed, SECURITY_LEVEL/4);
        for( i = 0 ; i < SECURITY_LEVEL/4 ; ++i )
        {
            key_seed2[i] = key_seed[i];
        }
        bacuov_keygen(&sk, &pk, seed);
    }

    /* report on results */
    num_successes = 0;
    printf("Ran %i trials with %i successes and %i failures.\n", num_trials, num_successes, num_failures);
    printf("Successes:\n");
    printf(" * %i total successes\n", num_successes);

    return 0;
}

