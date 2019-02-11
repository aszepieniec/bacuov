#include "gfpe.h"
#include "gfpem.h"

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

    gfpx modulus;

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

    modulus = gfpe_init_defining_polynomial();

    gfpe_inverse(&ainv, a);
    gfpe_inverse(&binv, b);


    gfpe_multiply(&cc, ainv, binv);

    gfpe_multiply(&cinv, cc, c);

    if( gfpe_is_one(cinv) == 1 )
    {
        printf("success.\n");
        gfpx_destroy(modulus);
        return 1;
    }
    else
    {
        printf("failure!\n");
    }

    gfpx_destroy(modulus);

    return 0;
}

int test_gfpem_inverse( unsigned long int * rand )
{
    unsigned char data[sizeof(unsigned long int)];
    unsigned char buff[EXTENSION_DEGREE * 7*7 * 3];
    int i;
    gfpematrix a, b, c;
    gfpematrix ainv, binv, cinv;
    gfpematrix cc;
    int verbose;
    int invertible;

    verbose = 0;

    FIPS202_SHAKE256((unsigned char *)rand, sizeof(unsigned long int), buff, sizeof(buff));

    printf("testing gfpematrix inverse with randomness %lu ... ", *rand);

    for( i = 0 ; i < sizeof(unsigned long int) ; ++i )
        ((unsigned char *)rand)[i] = buff[i];

    a = gfpem_init(7,7);
    b = gfpem_init(7,7);
    c = gfpem_init(7,7);
    ainv = gfpem_init(7,7);
    binv = gfpem_init(7,7);
    cinv = gfpem_init(7,7);
    cc = gfpem_init(7,7);

    gfpem_random(a, buff + sizeof(unsigned long int));
    gfpem_random(b, buff + sizeof(unsigned long int) + EXTENSION_DEGREE * 7*7);

    gfpem_multiply(&c, a, b);


    invertible = 1;
    invertible &= gfpem_inverse(ainv, a);
    invertible &= gfpem_inverse(binv, b);
    invertible &= gfpem_inverse(cinv, c);

    gfpem_multiply(&cc, binv, ainv);

    if( invertible == 0 )
    {
        printf("some matrix is not invertible :(\n");
    }
    else if( gfpem_equals(cc, cinv) == 1 )
    {
        printf("success.\n");
    }
    else
    {
        printf("failure!\n");

        gfpem_destroy(a);
        gfpem_destroy(b);
        gfpem_destroy(c);
        gfpem_destroy(ainv);
        gfpem_destroy(binv);
        gfpem_destroy(cinv);
        gfpem_destroy(cc);

        return 0;
    }

    gfpem_destroy(a);
    gfpem_destroy(b);
    gfpem_destroy(c);
    gfpem_destroy(ainv);
    gfpem_destroy(binv);
    gfpem_destroy(cinv);
    gfpem_destroy(cc);

    return 1;
}

int main( int argc, char ** argv )
{
    unsigned long int random;
    int i;
    int cont;

    random = rand();

    cont = 1;
    for( i = 0 ; i < 100 && cont == 1 ; ++i ) cont = cont & test_gfpe_inverse(&random);
    for( i = 0 ; i < 100 && cont == 1 ; ++i ) cont = cont & test_gfpem_inverse(&random);

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

