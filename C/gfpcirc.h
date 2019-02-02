#ifndef GFPCIRC
#define GFPCIRC

#include "gfp.h"

#ifndef DEGREE_OF_CIRCULANCY
#warning "Degree of circulancy not defined, setting to 11."
#define DEGREE_OF_CIRCULANCY 11
#endif

typedef struct 
{
    gfp_element data[DEGREE_OF_CIRCULANCY];
} gfpcirc_element;


gfpcirc_element gfpcirc( int castee );
gfpcirc_element gfpcirc_init( unsigned int size );
gfpcirc_element gfpcirc_clone( gfpcirc_element elm );
int gfpcirc_destroy( gfpcirc_element elm );
int gfpcirc_zero( gfpcirc_element* elm );
int gfpcirc_one( gfpcirc_element* elm );
int gfpcirc_random( gfpcirc_element* elm, unsigned char * randomness );
int gfpcirc_random_invertible( gfpcirc_element* elm, unsigned char * randomness );
int gfpcirc_compare( gfpcirc_element lhs, gfpcirc_element rhs );
int gfpcirc_copy( gfpcirc_element * dest, gfpcirc_element * source );
int gfpcirc_add( gfpcirc_element * res, gfpcirc_element lhs, gfpcirc_element rhs );
int gfpcirc_subtract( gfpcirc_element * res, gfpcirc_element lhs, gfpcirc_element rhs );
int gfpcirc_negate( gfpcirc_element * res, gfpcirc_element elm );
int gfpcirc_multiply( gfpcirc_element * res, gfpcirc_element lhs, gfpcirc_element rhs );
int gfpcirc_divide( gfpcirc_element * quo, gfpcirc_element numerator, gfpcirc_element divisor );
int gfpcirc_flip( gfpcirc_element * elm );
int gfpcirc_shift( gfpcirc_element * elm, int shamt );
int gfpcirc_print( gfpcirc_element * elm );
int gfpcirc_is_one( gfpcirc_element elm );
int gfpcirc_is_zero( gfpcirc_element elm );

#endif

