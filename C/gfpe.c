#include "gfpe.h"

#include <stdio.h>
#include <stdlib.h>

#if GFP_NUMBYTES <= 4

gfpx gfpe_init_defining_polynomial()
{
    int i;
    gfpx elm;
    elm.data = malloc(EXTENSION_DEGREE+1);
    elm.degree = EXTENSION_DEGREE;
    for( i = 0 ; i < EXTENSION_DEGREE ; ++i )
    {
        elm.data[i] = (DEFINING_POLYNOMIAL >> (8*i)) & 0xff;
    }
    elm.data[EXTENSION_DEGREE] = 0x01;
    return elm;
}

gfpe_element gfpe( int castee )
{
    gfpe_element e;
    int i;
    e.data[0] = gfp(((castee % GF_PRIME_MODULUS) + GF_PRIME_MODULUS) % GF_PRIME_MODULUS);
    for( i = 1 ; i < EXTENSION_DEGREE ; ++i )
    {
        e.data[0] = 0;
    }
    return e;
}

gfpe_element gfpe_clone( gfpe_element elm )
{
    return elm;
}

int gfpe_copy( gfpe_element * dest, gfpe_element source )
{
    int i;
    for( i = 0 ; i < EXTENSION_DEGREE ; ++i )
    {
        dest->data[i] = source.data[i];
    }
    return 1;
}

int gfpe_zero( gfpe_element* elm )
{
    int i;
    for( i = 0 ; i < EXTENSION_DEGREE ; ++i )
    {
        elm->data[i] = 0;
    }
    return 1;
}

int gfpe_one( gfpe_element* elm )
{
    int i;
    elm->data[0] = 1;
    for( i = 1 ; i < EXTENSION_DEGREE ; ++i )
    {
        elm->data[i] = 0;
    }
    return 1;
}

int gfpe_random( gfpe_element* elm, unsigned char * randomness )
{
    int j;
    unsigned long int r = 0; 
    for( j = 0 ; j < EXTENSION_DEGREE * GFP_NUMBYTES ; ++j )
    {
        r = r * 256 + randomness[j];
    }
    for( j = 0 ; j < EXTENSION_DEGREE ; ++j )
    {
        elm->data[j] = r % GF_PRIME_MODULUS;
        r = r / GF_PRIME_MODULUS;
    }
    return 1;
}

int gfpe_compare( gfpe_element lhs, gfpe_element rhs )
{
    int i;
    for( i = 0 ; i < EXTENSION_DEGREE ; ++i )
    {
        if( (int)(lhs.data[i]) != (int)(rhs.data[i]) )
        {
            return 0;
        }
    }
    return 1;
}

int gfpe_is_one( gfpe_element elm )
{
    int i;
    if( elm.data[0] != 1 )
    {
        return 0;
    }
    for( i = 1 ; i < EXTENSION_DEGREE ; ++i )
    {
        if( elm.data[i] != 0 )
        {
            return 0;
        }
    }
    return 1;
}

int gfpe_is_zero( gfpe_element elm )
{
    int i;
    for( i = 0 ; i < EXTENSION_DEGREE ; ++i )
    {
        if( elm.data[i] != 0 )
        {
            return 0;
        }
    }
    return 1;
}

int gfpe_add( gfpe_element * res, gfpe_element lhs, gfpe_element rhs )
{
    int i;
    for( i = 0 ; i < EXTENSION_DEGREE ; ++i )
    {
        res->data[i] = ((int)(lhs.data[i]) + (int)(rhs.data[i])) % GF_PRIME_MODULUS;
    }
    return 1;
}

int gfpe_subtract( gfpe_element * res, gfpe_element lhs, gfpe_element rhs )
{
    int i;
    for( i = 0 ; i < EXTENSION_DEGREE ; ++i )
    {
        res->data[i] = (GF_PRIME_MODULUS + (int)(lhs.data[i]) - (int)(rhs.data[i])) % GF_PRIME_MODULUS;
    }
    return 1;
}

int gfpe_negate( gfpe_element * res, gfpe_element elm )
{
    int i;
    for( i = 0 ; i < EXTENSION_DEGREE ; ++i )
    {
        res->data[i] = (GF_PRIME_MODULUS - (int)(elm.data[i])) % GF_PRIME_MODULUS;
    }
    return 1;
}

int gfpe_multiply( gfpe_element * res, gfpe_element lhs, gfpe_element rhs )
{
    gfpx product;
    gfpx rem, quo, defpoly;
    int i, j;

    if( gfpe_is_zero(lhs) || gfpe_is_zero(rhs) )
    {
        gfpe_zero(res);
        return 1;
    }

    product = gfpx_init(2*EXTENSION_DEGREE);

    for( i = 0 ; i <= 2*EXTENSION_DEGREE ; ++i )
    {
        product.data[i] = 0;
    }
    for( i = 0 ; i < EXTENSION_DEGREE ; ++i )
    {
        for( j = 0 ; j < EXTENSION_DEGREE ; ++j )
        {
            product.data[i+j] = (product.data[i+j] + (lhs.data[i] * rhs.data[j])) % GF_PRIME_MODULUS;
        }
    }
    product.degree = 2*EXTENSION_DEGREE;

    // reduce modulo DEFINING_POLYNOMIAL
    rem = gfpx_init(0);
    quo = gfpx_init(0);
    defpoly = gfpe_init_defining_polynomial();
    gfpx_divide(&quo, &rem, product, defpoly);
    
    // copy rem over
    for( i = 0 ; i < EXTENSION_DEGREE ; ++i )
    {
        if( i <= rem.degree )
        {
            res->data[i] = rem.data[i];
        }
        else
        {
            res->data[i] = 0;
        }
    }

    // free allocated data
    gfpx_destroy(product);
    gfpx_destroy(rem);
    gfpx_destroy(quo);
    gfpx_destroy(defpoly);

    return 1;
}

int gfpe_inverse( gfpe_element * inv, gfpe_element elm )
{
    int i;
    gfpx a, b, g, element, defpoly;
    
    a = gfpx_init(0);
    b = gfpx_init(0);
    g = gfpx_init(0);
    element = gfpx_init(EXTENSION_DEGREE-1);
    defpoly = gfpe_init_defining_polynomial();

    for( i = 0 ; i < EXTENSION_DEGREE ; ++i )
    {
        element.data[i] = elm.data[i];
    }

    gfpx_xgcd(&a, &b, &g, element, defpoly);

    for( i = 0 ; i < EXTENSION_DEGREE ; ++i )
    {
        if( i <= a.degree )
        {
            inv->data[i] = a.data[i];
        }
        else
        {
            inv->data[i] = 0;
        }
    }

    gfpx_destroy(a);
    gfpx_destroy(b);
    gfpx_destroy(g);
    gfpx_destroy(element);
    gfpx_destroy(defpoly);
}

int gfpe_print( gfpe_element elm )
{
    int i;
    for( i = 0 ; i < EXTENSION_DEGREE ; ++i )
    {
        printf("%i", elm.data[i]);
        //if( elm.data[i] != 0 )
        //{
        //    printf("%i*x^%i", elm.data[i], i);
        //    if( i != EXTENSION_DEGREE - 1 && elm.data[i+1] != 0 )
        //    {
        //        printf(" + ");
        //    }
        //}
    }
    return 1;
}

#endif


