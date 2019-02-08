#include "gfpe.h"

#include <stdio.h>
#include <stdlib.h>

void FIPS202_SHAKE256(const unsigned char *input, unsigned int inputByteLen, unsigned char *output, int outputByteLen);

int test_gfpe_inverse( unsigned long int * rand )
{
    unsigned char data[sizeof(unsigned long int)];
    unsigned char buff[1000];
    int i;
    gfpe_element a, b, c;
    gfpe_element ainv, binv, cinv;
    gfpe_element cc;
    int verbose;

    verbose = 0;
    //*rand = 12473469396624703114;
    //if( *rand == 12473469396624703114 )
    //if( *rand == 13849300305790542478 )
    //{
    //    printf("got right randomness?");
    //    getchar();
    //    verbose = 1;
    //}

    FIPS202_SHAKE256((unsigned char *)rand, sizeof(unsigned long int), buff, 100);

    printf("testing gfpe inverse with randomness %lu ... ", *rand);

    for( i = 0 ; i < sizeof(unsigned long int) ; ++i )
        ((unsigned char *)rand)[i] = buff[i];

    gfpe_random(&a, buff + sizeof(unsigned long int));
    gfpe_random(&b, buff + sizeof(unsigned long int) + EXTENSION_DEGREE);

    gfpe_multiply(&c, a, b);

    if( gfpe_is_zero(c) )
    {
        printf("(sampled zero element; can't try inverse)\n");
        return 1;
    }

    gfpe_inverse(&ainv, a);
    gfpe_inverse(&binv, b);

    gfpe_inverse(&cinv, c);

    gfpe_multiply(&cc, ainv, binv);

    if( gfpe_compare(cc, c) == 0 )
    {
        printf("success.\n");
        return 1;
    }
    else
    {
        printf("failure!\n");
    }

    return 1;
}

int main( int argc, char ** argv )
{
    unsigned long int random;
    int i;
    int cont;

    random = rand();

    cont = 1;
    for( i = 0 ; i < 1000 && cont == 1 ; ++i ) cont &= test_gfpe_inverse(&random);

    printf("All tests done.\n");
    if( cont == 1 )
    {
        printf("success \\o/\n");
    }
    else
    {
        printf("failure :(\n");
    }

    return 1;
}

