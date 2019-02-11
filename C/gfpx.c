#include "gfpx.h"
#include <stdlib.h>
#include <stdio.h>


/**
 * gfpx_init
 * Initialize a GF(p)[x] object of given degree. Allocate memory
 * and set to zero.
 */
gfpx gfpx_init( int deg )
{
    gfpx elm;
    elm.data = malloc(deg+1);
    elm.degree = deg;
    return elm;
}

/**
 * gfpx_zero
 * Set the given polynomial to zero.
 */
int gfpx_zero( gfpx* p )
{
    free(p->data);
    p->degree = 0;
    p->data = malloc(1);
    p->data[0] = 0;
    return 1;
}

/**
 * gfpx_one
 * Set the given polynomial to one.
 */
int gfpx_one( gfpx* p )
{
    free(p->data);
    p->degree = 0;
    p->data = malloc(1);
    p->data[0] = 1;
    return 1;
}

/**
 * gfpx_copy
 * Copy a GF(p)[x] element from one container to another, and
 * reinitialize as necessary.
 */
int gfpx_copy( gfpx* dest, gfpx source )
{
    int i;
    if( dest->degree != source.degree )
    {
        free(dest->data);
        dest->data = malloc(source.degree+1);
        dest->degree = source.degree;
    }
    for( i = 0 ; i <= source.degree ; ++i )
    {
        dest->data[i] = source.data[i];
    }
    return 1;
}

/**
 * gfpx_destroy
 * Destroy a GF(p)[x] object. Free memory.
 */
int gfpx_destroy( gfpx p )
{
    free(p.data);
}

/**
 * gfpx_add
 * Add two GF(p)[x] elements together.
 */
int gfpx_add( gfpx* dest, gfpx lhs, gfpx rhs )
{
    int i;
    unsigned char * data;
    if( rhs.degree > lhs.degree )
    {
        return gfpx_add(dest, rhs, lhs);
    }

    data = malloc(lhs.degree+1);

    for( i = 0 ; i <= rhs.degree ; ++i )
    {
        data[i] = (lhs.data[i] + rhs.data[i]) % GF_PRIME_MODULUS;
    }
    for( ; i <= lhs.degree ; ++i )
    {
        data[i] = lhs.data[i];
    }

    free(dest->data);
    dest->degree = lhs.degree;
    dest->data = data;

    while( data[dest->degree] == 0 && dest->degree > 0 )
    {
        dest->degree -= 1;
    }

    return 1;
}

/**
 * gfpx_subtract
 * Subtract the second GF(p)[x] element from the first.
 */
int gfpx_subtract( gfpx* dest, gfpx lhs, gfpx rhs )
{
    int i;
    unsigned char * data;
    int maxdeg;

    maxdeg = lhs.degree;
    if( rhs.degree > lhs.degree )
    {
        maxdeg = rhs.degree;
    }

    data = malloc(maxdeg+1);

    for( i = 0 ; i <= maxdeg ; ++i )
    {
        if( i <= lhs.degree && i <= rhs.degree )
        {
            data[i] = (GF_PRIME_MODULUS + lhs.data[i] - rhs.data[i]) % GF_PRIME_MODULUS;
        }
        else if( i > lhs.degree )
        {
            data[i] = (GF_PRIME_MODULUS - rhs.data[i]) % GF_PRIME_MODULUS;
        }
        else if( i > rhs.degree )
        {
            data[i] = lhs.data[i];
        }
    }

    free(dest->data);
    dest->degree = maxdeg;
    dest->data = data;

    while( data[dest->degree] == 0 && dest->degree > 0 )
    {
        dest->degree -= 1;
    }

    return 1;
}

/**
 * gfpx_multiply
 * Multiply two GF(p)[x] elements together.
 */
int gfpx_multiply( gfpx* dest, gfpx lhs, gfpx rhs )
{
    int i, j;
    int degree;
    unsigned char * data;
    gfp_element product;

    degree = lhs.degree + rhs.degree;
    data = malloc(degree + 1);

    for( i = 0 ; i <= degree ; ++i )
    {
        data[i] = 0;
    }

    for( i = 0 ; i <= lhs.degree ; ++i )
    {
        for( j = 0 ; j <= rhs.degree ; ++j )
        {
            gfp_multiply(&product, lhs.data[i], rhs.data[j]);
            data[i+j] = (data[i+j] + product) % GF_PRIME_MODULUS;
        }
    }

    free(dest->data);
    dest->data = data;
    dest->degree = degree;

    return 1;
}

/**
 * gfpx_equals
 * Decide if two elements of GF(p)[x] are equal, and return 1 if so.
 * (Return 0 otherwise.)
 */
int gfpx_equals( gfpx lhs, gfpx rhs )
{
    int i;
    int equal;
    if( lhs.degree != rhs.degree )
    {
        return 0;
    }
    equal = 1;
    for( i = 0 ; i <= lhs.degree ; ++i )
    {
        equal = equal & (lhs.data[i] == rhs.data[i]);
    }
    return equal;
}

/**
 * gfpx_is_zero
 * Determine if the given polynomial is equal to zero. Return one if
 * so, zero otherwise.
 */
int gfpx_is_zero( gfpx p )
{
    int zero;
    int i;
    zero = 1;
    for( i = 0 ; i <= p.degree ; ++i )
    {
        zero &= p.data[i] == 0;
    }
    return zero;
}

/**
 * gfpx_multiply_constant_shift
 * Multiply the polynomial with a constant and shift it (to the left,
 * i.e., towards higher degree). Satisfies:
 *  dest == constant * x^shift * poly
 */
int gfpx_multiply_constant_shift( gfpx* dest, gfpx poly, unsigned char constant, int shift )
{
    unsigned char * data;
    int i;
    int degree;

    degree = shift + poly.degree;

    data = malloc(degree+1);
    for( i = 0 ; i <= degree ; ++i )
    {
        data[i] = 0;
    }

    for( i = shift ; i <= degree ; ++i )
    {
        gfp_multiply(&data[i], poly.data[i-shift], constant);
    }

    free(dest->data);
    dest->data = data;
    dest->degree = degree;

    return 1;
}

/**
 * gfpx_divide
 * Divide one GF(p)[x] element by another and record the quotient
 * and remainder.
 */
int gfpx_divide( gfpx* quo, gfpx* rem, gfpx num, gfpx divisor )
{
    int i, j;
    unsigned char inv, compl;
    gfpx remainder, poly;
    unsigned char * quotient_data;

    /* make sure divisor leading coefficient is not zero */
    if( divisor.data[divisor.degree] == 0 )
    {
        poly.data = malloc(divisor.degree);
        for( i = 0 ; i < divisor.degree ; ++i )
        {
            poly.data[i] = divisor.data[i];
        }
        for( poly.degree = divisor.degree-1 ; poly.degree > 0 ; --poly.degree )
        {
            if( poly.data[poly.degree] != 0 )
            {
                break;
            }
        }
        gfpx_divide(quo, rem, num, poly);
        free(poly.data);
        return 1;
    }

    /* make sure numerator leading coefficient is not zero */
    if( num.data[num.degree] == 0 )
    {
        poly.data = malloc(num.degree);
        for( i = 0 ; i < num.degree ; ++i )
        {
            poly.data[i] = num.data[i];
        }
        for( poly.degree = num.degree-1 ; poly.degree > 0 ; --poly.degree )
        {
            if( poly.data[poly.degree] != 0 )
            {
                break;
            }
        }
        gfpx_divide(quo, rem, poly, divisor);
        free(poly.data);
        return 1;
    }

    /* make sure deg(divisor) <= deg(numerator) */
    if( divisor.degree > num.degree )
    {
        gfpx_zero(quo);
        gfpx_copy(rem, num);
        return 1;
    }

    /* filtered out edge cases, proceed with division already */
    remainder = gfpx_init(0);
    poly = gfpx_init(0);
    gfpx_copy(&remainder, num);
    if( num.degree - divisor.degree + 1 == -1 )
    {
        printf("should never get here.\n");
        getchar();
    }
    quotient_data = malloc(num.degree - divisor.degree + 1);
    for( i = 0 ; i <= num.degree - divisor.degree ; ++i )
    {
        quotient_data[i] = 0;
    }

    gfp_inverse(&inv, divisor.data[divisor.degree]);

    for( i = remainder.degree - divisor.degree ; i >= 0 ; --i )
    {
        if( remainder.degree < divisor.degree + i )
        {
            continue;
        }


        gfp_multiply(&compl, remainder.data[remainder.degree], inv);
        gfpx_multiply_constant_shift(&poly, divisor, compl, i);


        quotient_data[i] = compl;
        gfpx_subtract(&remainder, remainder, poly);

    }

    free(quo->data);
    quo->data = quotient_data;
    quo->degree = num.degree - divisor.degree;

    gfpx_copy(rem, remainder);
    gfpx_destroy(remainder);
    gfpx_destroy(poly);

    return 1;
}

/**
 * gfpx_xgcd
 * Compute the greatest common divisor g and Bezout coefficients a
 * and b for x and y using the extended Euclidean algorithm.
 */
int gfpx_xgcd( gfpx* a, gfpx* b, gfpx* g, gfpx x, gfpx y )
{
    gfpx s, old_s;
    gfpx t, old_t;
    gfpx r, old_r;
    gfpx quotient, remainder;
    gfpx temp;
    gfp_element lc;

    s = gfpx_init(0);
    old_s = gfpx_init(0);
    t = gfpx_init(0);
    old_t = gfpx_init(0);
    r = gfpx_init(0);
    old_r = gfpx_init(0);
    quotient = gfpx_init(0);
    remainder = gfpx_init(0);
    temp = gfpx_init(0);

    gfpx_zero(&s);
    gfpx_one(&old_s);
    gfpx_one(&t);
    gfpx_zero(&old_t);
    gfpx_copy(&r, y);
    gfpx_copy(&old_r, x);

    while( gfpx_is_zero(r) == 0 )
    {
        gfpx_divide(&quotient, &remainder, old_r, r);

        gfpx_copy(&old_r, r);
        gfpx_copy(&r, remainder);

        gfpx_multiply(&temp, quotient, s);
        gfpx_subtract(&temp, old_s, temp);
        gfpx_copy(&old_s, s);
        gfpx_copy(&s, temp);

        gfpx_multiply(&temp, quotient, t);
        gfpx_subtract(&temp, old_t, temp);
        gfpx_copy(&old_t, t);
        gfpx_copy(&t, temp);
    }

    lc = old_r.data[old_r.degree];
    gfp_inverse(&lc, lc);
    gfpx_multiply_constant_shift(a, old_s, lc, 0);
    gfpx_multiply_constant_shift(b, old_t, lc, 0);
    gfpx_multiply_constant_shift(g, old_r, lc, 0);

    gfpx_destroy(s);
    gfpx_destroy(old_s);
    gfpx_destroy(t);
    gfpx_destroy(old_t);
    gfpx_destroy(r);
    gfpx_destroy(old_r);
    gfpx_destroy(quotient);
    gfpx_destroy(remainder);
    gfpx_destroy(temp);

    return 1;
}

/**
 * gfpx_print
 * Cast the polynomial's coefficients to hex number and throw them to
 * stdout.
 */
int gfpx_print( gfpx p )
{
    int i;

    for( i = 0 ; i < p.degree ; ++i )
    {
        printf("%i*x^%i + ", p.data[i], i);
    }
    printf("%i*x^%i", p.data[p.degree], p.degree);

    return 1;
}

