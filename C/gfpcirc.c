#include "gfpcirc.h"

#include <stdio.h>

#if GFP_NUMBYTES <= 4

gfpcirc_element gfpcirc( int castee )
{
    gfpcirc_element e;
    e = gfpcirc_init(sizeof(castee));
    e.data[0] = gfp(((castee % GF_PRIME_MODULUS) + GF_PRIME_MODULUS) % GF_PRIME_MODULUS);
    return e;
}

gfpcirc_element gfpcirc_init( unsigned int size )
{
    int i;
    gfpcirc_element e;

    for( i = 0 ; i < DEGREE_OF_CIRCULANCY ; ++i )
    {
        e.data[i] = gfp_init(0);
    }
    return e;
}

gfpcirc_element gfpcirc_clone( gfpcirc_element elm )
{
    return elm;
}

int gfpcirc_destroy( gfpcirc_element elm )
{
    return 1;
}

int gfpcirc_copy( gfpcirc_element * dest, gfpcirc_element * source )
{
    int i;
    for( i = 0 ; i < DEGREE_OF_CIRCULANCY ; ++i )
    {
        dest->data[i] = source->data[i];
    }
    return 1;
}

int gfpcirc_zero( gfpcirc_element* elm )
{
    int i;
    for( i = 0 ; i < DEGREE_OF_CIRCULANCY ; ++i )
    {
        gfp_zero( &(elm->data[i]) );
    }
    return 1;
}

int gfpcirc_one( gfpcirc_element* elm )
{
    int i;
    gfp_zero( &(elm->data[0]) );
    for( i = 1 ; i < DEGREE_OF_CIRCULANCY ; ++i )
    {
        gfp_zero( &(elm->data[i]) );
    }
    return 1;
}

int gfpcirc_random( gfpcirc_element* elm, unsigned char * randomness )
{
    int j;
    unsigned long int r = 0; /* TODO: invoke big integer help if necessary */
    for( j = 0 ; j < DEGREE_OF_CIRCULANCY * GFP_NUMBYTES ; ++j )
    {
        r = r * 256 + randomness[j];
    }
    for( j = 0 ; j < DEGREE_OF_CIRCULANCY ; ++j )
    {
        elm->data[j] = r % GF_PRIME_MODULUS;
        r = r / GF_PRIME_MODULUS;
    }
    return 1;
}

int gfpcirc_compare( gfpcirc_element lhs, gfpcirc_element rhs )
{
    int i;
    for( i = 0 ; i < DEGREE_OF_CIRCULANCY ; ++i )
    {
        if( (int)(lhs.data[i]) != (int)(rhs.data[i]) )
        {
            return 0;
        }
    }
    return 1;
}

int gfpcirc_is_one( gfpcirc_element elm )
{
    int i;
    if( elm.data[0] != 1 )
    {
        return 0;
    }
    for( i = 1 ; i < DEGREE_OF_CIRCULANCY ; ++i )
    {
        if( elm.data[i] != 0 )
        {
            return 0;
        }
    }
    return 1;
}

int gfpcirc_is_zero( gfpcirc_element elm )
{
    int i;
    for( i = 0 ; i < DEGREE_OF_CIRCULANCY ; ++i )
    {
        if( elm.data[i] != 0 )
        {
            return 0;
        }
    }
    return 1;
}

int gfpcirc_add( gfpcirc_element * res, gfpcirc_element lhs, gfpcirc_element rhs )
{
    int i;
    for( i = 0 ; i < DEGREE_OF_CIRCULANCY ; ++i )
    {
        res->data[i] = ((int)(lhs.data[i]) + (int)(rhs.data[i])) % GF_PRIME_MODULUS;
    }
    return 1;
}

int gfpcirc_subtract( gfpcirc_element * res, gfpcirc_element lhs, gfpcirc_element rhs )
{
    int i;
    for( i = 0 ; i < DEGREE_OF_CIRCULANCY ; ++i )
    {
        res->data[i] = (GF_PRIME_MODULUS + (int)(lhs.data[i]) - (int)(lhs.data[i])) % GF_PRIME_MODULUS;
    }
    return 1;
}

int gfpcirc_negate( gfpcirc_element * res, gfpcirc_element elm )
{
    int i;
    for( i = 0 ; i < DEGREE_OF_CIRCULANCY ; ++i )
    {
        res->data[i] = (GF_PRIME_MODULUS - (int)(elm.data[i])) % GF_PRIME_MODULUS;
    }
    return 1;
}

int gfpcirc_multiply( gfpcirc_element * res, gfpcirc_element lhs, gfpcirc_element rhs )
{
    gfp_element data[2*DEGREE_OF_CIRCULANCY];
    int i, j;
    for( i = 0 ; i < 2*DEGREE_OF_CIRCULANCY ; ++i )
    {
        data[i] = 0;
    }
    for( i = 0 ; i < DEGREE_OF_CIRCULANCY ; ++i )
    {
        for( j = 0 ; j < DEGREE_OF_CIRCULANCY ; ++j )
        {
            data[i+j] = (data[i+j] + (lhs.data[i] * rhs.data[j])) % GF_PRIME_MODULUS;
        }
    }
    /* reduce modulo x^l-1 */
    for( i = 0 ; i < DEGREE_OF_CIRCULANCY ; ++i )
    {
        res->data[i] = ((int)data[i] + (int)data[DEGREE_OF_CIRCULANCY+i]) % GF_PRIME_MODULUS;
    }

    return 1;
}

int gfpcirc_flip( gfpcirc_element * elm )
{
    int i;
    gfpcirc_element temp;
    for( i = 0 ; i < DEGREE_OF_CIRCULANCY ; ++i )
    {
        temp.data[i] = elm->data[i];
    }
    for( i = 0 ; i < DEGREE_OF_CIRCULANCY ; ++i )
    {
        elm->data[DEGREE_OF_CIRCULANCY-1-i] = temp.data[i];
    }
    return 1;
}

int gfpcirc_shift( gfpcirc_element * elm, int shamt )
{
    int i;
    gfpcirc_element temp;
    
    if( shamt == 0 )
    {
        return 1;
    }
    if( shamt < 0 )
    {
        i = (-shamt) / DEGREE_OF_CIRCULANCY;
        return gfpcirc_shift(elm, (i+1)*DEGREE_OF_CIRCULANCY + shamt);
    }
    if( shamt > DEGREE_OF_CIRCULANCY )
    {
        return gfpcirc_shift(elm, shamt % DEGREE_OF_CIRCULANCY);
    }

    for( i = 0 ; i < DEGREE_OF_CIRCULANCY ; ++i )
    {
        temp.data[i] = elm->data[i];
    }
    for( i = 0 ; i < DEGREE_OF_CIRCULANCY ; ++i )
    {
        elm->data[(i+shamt) % DEGREE_OF_CIRCULANCY] = temp.data[i];
    }
    return 1;
}

int gfpcirc_print( gfpcirc_element * elm )
{
    int i;
    printf("(");
    for( i = 0 ; i < DEGREE_OF_CIRCULANCY ; ++i )
    {
        printf("%i,", (int)elm->data[i]);
    }
    printf(")");
    return 1;
}

#endif


