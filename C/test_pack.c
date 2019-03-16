#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bacuov.h"

int FIPS202_SHAKE256( unsigned char * source, int source_length, unsigned char * dest, int dest_length );

int main( int argc, char ** argv )
{
   
    unsigned char * seed;
    unsigned char randomness[32];
    int seed_length;
    unsigned int num_elements;
    unsigned char * elements;
    unsigned char * packed;
    unsigned char * unpacked;
    int equals;
    int num_trials, trial_index;
    int i;

    if( argc < 3 )
    {
        printf("Usage: ./test_pack [num_trials] [random seed, in hex]\n");
        return 0;
    }

    // get seed
    seed_length = strlen(argv[2])/2;
    seed = malloc(seed_length);
    for( i = 0 ; i < seed_length ; ++i )
    {
        sscanf(argv[2] + 2*i, "%2hhx", &seed[i]);
    }

    // get trial number
    num_trials = atoi(argv[1]);

    // run trials
    equals = 1;
    FIPS202_SHAKE256(seed, seed_length, randomness, 32);
    free(seed);

    for( trial_index = 0 ; trial_index < num_trials ; ++trial_index )
    {
        // generate elements
        num_elements = randomness[0] + randomness[1]*256;
        FIPS202_SHAKE256(randomness, 32, randomness, 32);
        elements = malloc(num_elements);
        FIPS202_SHAKE256(randomness, 32, elements, num_elements);
        for( i = 0 ; i < num_elements ; ++i )
        {
            elements[i] = elements[i] % GF_PRIME_MODULUS;
        }

        // pack
        packed = malloc(bacuov_pack_length(num_elements));
        bacuov_pack(packed, elements, num_elements);

        // unpack
        unpacked = malloc(num_elements);
        bacuov_unpack(unpacked, packed, num_elements);

        // test for equality
        for( i = 0 ; i < num_elements ; ++i )
        {
            if( elements[i] != unpacked[i] )
            {
                equals = 0;
            }
        }

        // free memory
        free(elements);
        free(packed);
        free(unpacked);

        if( equals == 0 )
        {
            break;
        }
    }

    if( equals == 0 )
    {
        printf("Ran %i trials ... ", trial_index+1);
        printf("failure :(\n");
    }
    else
    {
        printf("Ran %i trials ... ", trial_index);
        printf("success! \\o/\n");
    }

    return 1;
}

