#ifndef GFP
#define GFP

#ifndef GF_PRIME_MODULUS
#warning "Prime modulus undefined; setting to 31."
#define GF_PRIME_MODULUS 31
#endif

#ifndef GFP_NUMBITES
#define GFP_NUMBITS 8
#endif

#ifndef GFP_NUMBYTES
#define GFP_NUMBYTES 1
#endif

typedef unsigned char gfp_element;

gfp_element gfp( int castee );
gfp_element gfp_init( unsigned int size );
gfp_element gfp_clone( gfp_element elm );
int gfp_destroy( gfp_element elm );
int gfp_copy( gfp_element* dest, gfp_element source );
int gfp_zero( gfp_element* elm );
int gfp_one( gfp_element* elm );
int gfp_random( gfp_element* elm, unsigned char * randomness );
int gfp_random_invertible( gfp_element* elm, unsigned char * randomness );
int gfp_compare( gfp_element lhs, gfp_element rhs );
int gfp_copy( gfp_element * dest, gfp_element source );
int gfp_add( gfp_element * res, gfp_element lhs, gfp_element rhs );
int gfp_subtract( gfp_element * res, gfp_element lhs, gfp_element rhs );
int gfp_negate( gfp_element * res, gfp_element elm );
int gfp_multiply( gfp_element * res, gfp_element lhs, gfp_element rhs );
int gfp_divide( gfp_element * quo, gfp_element numerator, gfp_element divisor );
int gfp_inverse( gfp_element * res, gfp_element elm );
int gfp_print( gfp_element elm );
int gfp_is_one( gfp_element elm );
int gfp_is_zero( gfp_element elm );

#endif

