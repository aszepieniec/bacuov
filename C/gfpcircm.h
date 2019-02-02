#ifndef GFPCIRCM_H
#define GFPCIRCM_H

/**
 * GFPMCIRC.H
 * Routines and structures for matrices over circulant rings, i.e.
 * Fq[x] / <x^l - 1>
 */

#include "gfpcirc.h"

typedef struct
{
    unsigned  int width;
    unsigned  int height;
    gfpcirc_element * data;
} gfpcircmatrix;

gfpcircmatrix gfpcircm( unsigned  int height, unsigned  int width, gfpcirc_element * pdata );

gfpcircmatrix gfpcircm_init( unsigned  int height, unsigned  int width );
int gfpcircm_destroy( gfpcircmatrix fm );
int gfpcircm_copy( gfpcircmatrix dest, gfpcircmatrix source );
gfpcircmatrix gfpcircm_clone( gfpcircmatrix source );

int gfpcircm_eye( gfpcircmatrix mat );
int gfpcircm_is_eye( gfpcircmatrix mat );
int gfpcircm_equals( gfpcircmatrix lhs, gfpcircmatrix rhs );
int gfpcircm_zeros( gfpcircmatrix mat );
int gfpcircm_random( gfpcircmatrix mat, unsigned char * randomness );
int gfpcircm_random_upper_triangular( gfpcircmatrix mat, unsigned char * randomness  );
int gfpcircm_transpose( gfpcircmatrix * mat );
int gfpcircm_add( gfpcircmatrix dest, gfpcircmatrix left, gfpcircmatrix right );
int gfpcircm_multiply( gfpcircmatrix * dest, gfpcircmatrix left, gfpcircmatrix right );
int gfpcircm_multiply_transpose( gfpcircmatrix * dest, gfpcircmatrix left, gfpcircmatrix rightT );
int gfpcircm_transpose_multiply( gfpcircmatrix * dest, gfpcircmatrix leftT, gfpcircmatrix right );
int gfpcircm_multiply_constant( gfpcircmatrix dest, gfpcircmatrix source, gfpcirc_element constant );
int gfpcircm_sum( gfpcircmatrix dest, gfpcircmatrix left_matrix, gfpcircmatrix right_matrix );
int gfpcircm_weighted_sum( gfpcircmatrix dest, gfpcirc_element left_constant, gfpcircmatrix left_matrix, gfpcirc_element right_constant, gfpcircmatrix right_matrix );
int gfpcircm_flip( gfpcircmatrix mat );
int gfpcircm_flipshift( gfpcircmatrix mat, int shamt );
int gfpcircm_shiftflip( gfpcircmatrix mat, int shamt );
int gfpcircm_shift( gfpcircmatrix mat, int shamt );
int gfpcircm_negate( gfpcircmatrix mat );

int gfpcircm_stack( gfpcircmatrix dest, gfpcircmatrix top, gfpcircmatrix bottom );
int gfpcircm_cat( gfpcircmatrix dest, gfpcircmatrix left, gfpcircmatrix right );
int gfpcircm_slice( gfpcircmatrix dest, gfpcircmatrix mat, unsigned  int row_start, unsigned  int col_start );

int gfpcircm_print( gfpcircmatrix mat );
int gfpcircm_print_transpose( gfpcircmatrix mat );

#endif

