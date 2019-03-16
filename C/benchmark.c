#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/random.h>
#include <time.h>
#include <x86intrin.h>

#include "api.h"

void randombytes( unsigned char * dest, unsigned int numbytes )
{
    getrandom(dest, numbytes, 0);
}

int main( int argc, char ** argv )
{
   
    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char sk[CRYPTO_SECRETKEYBYTES];
    unsigned char msg[32];
    unsigned char smsg[32 + CRYPTO_BYTES];
    unsigned long long smsglen, msglen;
    int valid;
    int num_trials, trial_index;
    int i;

    clock_t tick, tock;
    unsigned long long tack, tuck;
    clock_t time_keygen, time_signature, time_verification;
    unsigned long long cycles_keygen, cycles_signature, cycles_verification;

    // get trial number
    if( argc < 2 )
    {
        printf("Usage: ./benchmark [num_trials]\n");
        return 0;
    }
    num_trials = atoi(argv[1]);

    // run trials
    printf("Running %i trials ...\n", num_trials);

    valid = 1;
    time_keygen = 0;
    time_signature = 0;
    time_verification = 0;
    cycles_keygen = 0;
    cycles_signature = 0;
    cycles_verification = 0;
    for( trial_index = 0 ; trial_index < num_trials ; ++trial_index )
    {
        msglen = 32;
        randombytes(msg, msglen);

        tick = clock();
        tack = __rdtsc();
        crypto_sign_keypair(pk, sk);
        tock = clock();
        tuck = __rdtsc();
        time_keygen = time_keygen + (tock-tick);
        cycles_keygen = cycles_keygen + (tuck-tack);

        tick = clock();
        tack = __rdtsc();
        crypto_sign(smsg, &smsglen, msg, msglen, sk);
        tock = clock();
        tuck = __rdtsc();
        time_signature = time_signature + (tock-tick);
        cycles_signature = cycles_signature + (tuck-tack);

        tick = clock();
        tack = __rdtsc();
        valid = valid & (!crypto_sign_open(msg, &msglen, smsg, smsglen, pk));
        tock = clock();
        tuck = __rdtsc();
        time_verification = time_verification + (tock-tick);
        cycles_verification = cycles_verification + (tuck-tack);

        if( valid == 0 )
        {
            break;
        }
    }

    if( valid == 0 )
    {
        printf("Error! One or more signature were invalid.\n");
        return 0;
    }

    printf("Successfully ran %i trials.\nReport:\n", num_trials);
    printf("Keygen -- %f seconds, %llu cycles\n", (double)(time_keygen / CLOCKS_PER_SEC) / num_trials, cycles_keygen);
    printf("Sign -- %f seconds, %llu cycles\n", (double)(time_signature / CLOCKS_PER_SEC) / num_trials, cycles_signature);
    printf("Verify -- %f seconds, %llu cycles\n", (double)(time_verification / CLOCKS_PER_SEC) / num_trials, cycles_verification);

    return 1;
}

