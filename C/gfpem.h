#ifndef GFPEM_H
#define GFPEM_H

/**
 * GFPEM.H
 * Routines and structures for matrices over extensions finite fields.
 */

#include "gfpe.h"

typedef struct
{
    unsigned  int width;
    unsigned  int height;
    gfpe_element * data;
} gfpematrix;

gfpematrix gfpem( unsigned  int height, unsigned  int width, gfpe_element * pdata );

gfpematrix gfpem_init( unsigned  int height, unsigned  int width );
int gfpem_destroy( gfpematrix fm );
int gfpem_copy( gfpematrix dest, gfpematrix source );
gfpematrix gfpem_clone( gfpematrix source );

int gfpem_eye( gfpematrix mat );
int gfpem_is_eye( gfpematrix mat );
int gfpem_equals( gfpematrix lhs, gfpematrix rhs );
int gfpem_zeros( gfpematrix mat );
int gfpem_random( gfpematrix mat, unsigned char * randomness );
int gfpem_random_upper_triangular( gfpematrix mat, unsigned char * randomness  );
int gfpem_transpose( gfpematrix * mat );
int gfpem_add( gfpematrix dest, gfpematrix left, gfpematrix right );
int gfpem_multiply( gfpematrix * dest, gfpematrix left, gfpematrix right );
int gfpem_multiply_transpose( gfpematrix * dest, gfpematrix left, gfpematrix rightT );
int gfpem_transpose_multiply( gfpematrix * dest, gfpematrix leftT, gfpematrix right );
int gfpem_multiply_constant( gfpematrix dest, gfpematrix source, gfpe_element constant );
int gfpem_sum( gfpematrix dest, gfpematrix left_matrix, gfpematrix right_matrix );
int gfpem_weighted_sum( gfpematrix dest, gfpe_element left_constant, gfpematrix left_matrix, gfpe_element right_constant, gfpematrix right_matrix );
int gfpem_rowop( gfpematrix mat, unsigned  int destrow, unsigned  int sourcerow, gfpe_element constant, unsigned  int offset );
int gfpem_fliprows( gfpematrix mat, unsigned  int destrow, unsigned  int sourcerow );
int gfpem_scalerow( gfpematrix mat, unsigned  int rowidx, gfpe_element constant );
int gfpem_redech( gfpematrix mat );
int gfpem_solve( gfpematrix coeffs, gfpematrix target, gfpematrix solution, gfpematrix * kernel );
int gfpem_inspan( gfpematrix vec, gfpematrix mat );

int gfpem_stack( gfpematrix dest, gfpematrix top, gfpematrix bottom );
int gfpem_cat( gfpematrix dest, gfpematrix left, gfpematrix right );
int gfpem_slice( gfpematrix dest, gfpematrix mat, unsigned  int row_start, unsigned  int col_start );

int gfpem_inverse( gfpematrix dest, gfpematrix mat );

int gfpem_print( gfpematrix mat );
int gfpem_print_transpose( gfpematrix mat );

#endif

