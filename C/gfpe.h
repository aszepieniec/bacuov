#ifndef GFPE
#define GFPE

#include "gfpx.h"

#ifndef EXTENSION_DEGREE
#warning "Extension degree not defined, setting to 11."
#define EXTENSION_DEGREE 11
#endif

#ifndef DEFINING_POLYNOMIAL
#warning "No defining polynomial provided, setting to 2."
#define DEFINING_POLYNOMIAL 0x02 // without the monic leading coefficient
#endif

typedef struct 
{
    gfp_element data[EXTENSION_DEGREE];
} gfpe_element;

gfpx gfpe_init_defining_polynomial();

gfpe_element gfpe( int castee );
gfpe_element gfpe_clone( gfpe_element elm );
int gfpe_zero( gfpe_element* elm );
int gfpe_one( gfpe_element* elm );
int gfpe_random( gfpe_element* elm, unsigned char * randomness );
int gfpe_random_invertible( gfpe_element* elm, unsigned char * randomness );
int gfpe_compare( gfpe_element lhs, gfpe_element rhs );
int gfpe_copy( gfpe_element * dest, gfpe_element source );
int gfpe_add( gfpe_element * res, gfpe_element lhs, gfpe_element rhs );
int gfpe_subtract( gfpe_element * res, gfpe_element lhs, gfpe_element rhs );
int gfpe_negate( gfpe_element * res, gfpe_element elm );
int gfpe_multiply( gfpe_element * res, gfpe_element lhs, gfpe_element rhs );
int gfpe_divide( gfpe_element * quo, gfpe_element numerator, gfpe_element divisor );
int gfpe_inverse( gfpe_element * inv, gfpe_element elm );
int gfpe_print( gfpe_element elm );
int gfpe_is_one( gfpe_element elm );
int gfpe_is_zero( gfpe_element elm );

#endif

