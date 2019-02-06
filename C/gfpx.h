#ifndef GFPX_H
#define GFPX_H

#include "gfp.h"

typedef struct
{
    unsigned char * data;
    int degree;
} gfpx;

gfpx gfpx_init( int deg );
int gfpx_copy( gfpx* dest, gfpx source );
int gfpx_destroy( gfpx p );

int gfpx_add( gfpx* dest, gfpx lhs, gfpx rhs );
int gfpx_multiply( gfpx* dest, gfpx lhs, gfpx rhs );
int gfpx_multiply_constant_shift( gfpx* dest, gfpx poly, unsigned char constant, int shift );
int gfpx_divide( gfpx* quo, gfpx* rem, gfpx num, gfpx den );
int gfpx_xgcd( gfpx* a, gfpx* b, gfpx* g, gfpx x, gfpx y );

int gfpx_equals( gfpx lhs, gfpx rhs );
int gfpx_is_zero( gfpx op );

int gfpx_print( gfpx p );

#endif

