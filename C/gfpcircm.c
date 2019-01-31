#include "gfpcircm.h"

#include <stdlib.h>
#include <stdio.h>

/**
 * gfpcircm
 * Create gfpcircmatrix object with given buffer. Usefule for reusing
 * the same data line.
 */
gfpcircmatrix gfpcircm( unsigned  int height, unsigned  int width, gfpcirc_element* pdata )
{
    int i, j;
    gfpcircmatrix mat;
    mat.height = height;
    mat.width = width;
    mat.data = pdata;
    return mat;
}

/**
 * gfpcircm_init
 * Create a field matrix object and allocate space for it.
 */
gfpcircmatrix gfpcircm_init( unsigned  int height, unsigned  int width )
{
    int i, j;
    gfpcircmatrix mat;
    mat.width = width;
    mat.height = height;
    mat.data = malloc(width*height*sizeof(gfpcirc_element));
    for( i = 0 ; i < height ; ++i )
    {
        for( j = 0 ; j < width ; ++j )
        {
            mat.data[i*width + j] = gfpcirc_init(0);
        }
    }
    return mat;
}

/**
 * gfpcircm_destroy
 * Deallocates space allocated to a field matrix object. Call this
 * function before closing the scope where the field matrix object
 * was initialized.
 * @return
 *  * 1 if success
 */
int gfpcircm_destroy( gfpcircmatrix fm )
{
    free(fm.data);
    fm.width = 0;
    fm.height = 0;
    return 1;
}

/**
 * gfpcircm_copy
 * Copy the contents of one matrix to another. Does not allocate
 * memory for the new object; you must do that yourself! (Or use
 * gfpcircm_clone instead.)
 * @promise
 *  * dest.width >= source.width
 *  * dest.height >= source.height
 * @return
 *  * 1 if success
 */
int gfpcircm_copy( gfpcircmatrix dest, gfpcircmatrix source )
{
    unsigned int i, j;
    for( i = 0 ; i < source.height ; ++i )
    {
        for( j = 0 ; j < source.width ; ++j )
        {
            gfpcirc_copy(&dest.data[i*dest.width + j], &source.data[i*source.width + j]);
        }
    }

    return 1;
}

/**
 * gfpcircm_clone
 * Copy one matrix into a new object. Don't forget to call
 * gfpcircm_destroy at the end of scope!
 */
gfpcircmatrix gfpcircm_clone( gfpcircmatrix source )
{
    gfpcircmatrix mat;
    mat = gfpcircm_init(source.height, source.width);
    gfpcircm_copy(mat, source);
    return mat;
}

/**
 * gfpcircm_eye
 * Set a matrix to the identity matrix.
 * @param
 *  * eye : matrix object to set to identity
 * @returns
 *  * 1 if success, 0 otherwise
 */
int gfpcircm_eye( gfpcircmatrix eye )
{
    unsigned  int i, j;
    for( i = 0 ; i < eye.height ; ++i )
    {
        for( j = 0 ; j < eye.width ; ++j )
        {
            gfpcirc_zero(&eye.data[i*eye.width + j]);
        }
        gfpcirc_one(&eye.data[i*eye.width + i]);
    }
    return 1;
}

/**
 * Decide whether the given matrix is an identity matrix or, for
 * rectangular matrices, whether the main diagonal has ones and all
 * the rest is zero.
 * @return
 *  * 1 if identity, 0 otherwise
 */
int gfpcircm_is_eye( gfpcircmatrix eye )
{
    unsigned int i, j;
    int b = 1;
    gfpcirc_element one, zero;
    
    one = gfpcirc_init(1);
    zero = gfpcirc_init(1);
    gfpcirc_one(&one);
    gfpcirc_zero(&zero);

    for( i = 0 ; i < eye.height ; ++i )
    {
        for( j = 0 ; j < eye.width ; ++j )
        {
            if( i == j )
            {
                b = b & gfpcirc_compare(eye.data[i*eye.width + j], one);
            }
            else
            {
                b = b & gfpcirc_compare(eye.data[i*eye.width + j], zero);
            }
        }
    }

    gfpcirc_destroy(one);
    gfpcirc_destroy(zero);

    return b;
}

/**
 * gfpcircm_equals
 * Test two matrices for equality.
 * @return
 *  * 1 if equal, 0 otherwise
 */
int gfpcircm_equals( gfpcircmatrix lhs, gfpcircmatrix rhs )
{
    unsigned int i, j;
    int b;

    if( lhs.width != rhs.width || lhs.height != rhs.height )
    {
        return 0;
    }

    b = 1;
    for( i = 0 ; i < lhs.height ; ++i )
    {
        for( j = 0 ; j < lhs.width ; ++j )
        {
            b = b & gfpcirc_compare(lhs.data[lhs.width*i + j], rhs.data[rhs.width*i + j]);
        }
    }

    return b;
}

/**
 * gfpcircm_zero
 * Sets a matrix to all zeros.
 * @param
 *  * zero : matrix object to set to zero
 * @returns
 *  * 1 if success
 */
int gfpcircm_zeros( gfpcircmatrix zero )
{
    unsigned  int i, j;
    for( i = 0 ; i < zero.height ; ++i )
    {
        for( j = 0 ; j < zero.width ; ++j )
        {
            gfpcirc_zero(&zero.data[i*zero.width + j]);
        }
    }
    return 1;
}

/**
 * gfpcircm_random
 * Put random values into the matrix.
 * @params
 *  * random : matrix objects with data already allocated and whose
 *    elements are to be assigned random values
 *  * randomness : pointer to large-enough string of random bytes
 *    "large enough" means n*n*sizeof(unsigned int)
 * @result
 *  * random <-$- matrix_space(random.height, random.width)
 */
int gfpcircm_random( gfpcircmatrix random, unsigned char * randomness )
{
    unsigned  int i, j;
    unsigned int l;
    unsigned int * rand = (unsigned int *)randomness;
    unsigned int num_limbs;
    num_limbs = (GFP_NUMBITS + sizeof(unsigned long int)*8 - 1) / (sizeof(unsigned long int)*8);
   
    l = 0;
    for( i = 0 ; i < random.height ; ++i )
    {
        for( j = 0 ; j < random.width ; ++j )
        {
            gfpcirc_random(&random.data[i*random.width + j], randomness + l*num_limbs*sizeof(unsigned long int)); // TODO: BYTES PER GFPCIRC ELEMENT
            l = l + 1;
        }
    }

    return 1;
}

/**
 * gfpcircm_random_upper_triangular
 * Set the matrix to a random upper triangular with ones on the
 * diagonal.
 * @params
 *  * random : matrix objects with data already allocated and whose
 *    elements above the diagonal are to be assigned random values;
 *    whereas the elements above the diagonal are to be 0 and the
 *    elements on the diagonal are to be 1.
 *  * randomness : pointer to large-enough string of random bytes
 *    "large enough" means n*(n-1)/2*sizeof(unsigned int)
 * @result
 *  * random <-$- matrix_space(random.height, random.width)
 *    subject to
 *    forall i . random[i,i] = 1
 *    forall i, j>i . random[i,j] = 0
 */
int gfpcircm_random_upper_triangular( gfpcircmatrix random, unsigned char * randomness )
{
    unsigned  int i, j;
    unsigned int l;

    l = 0;
    for( i = 0 ; i < random.height ; ++i )
    {
        for( j = 0 ; j < i ; ++j )
        {
            gfpcirc_random(&random.data[i*random.width + j], &randomness[(l++)*(GFP_NUMBYTES+1)]);
        }
        gfpcirc_one(&random.data[i*random.width + i]);
        for( j = i+1 ; j < random.width ; ++j )
        {
            gfpcirc_zero(&random.data[i*random.width + j]);
        }
    }

    return 1;
}

/**
 * gfpcircm_transpose
 */
int gfpcircm_transpose( gfpcircmatrix * trans )
{
    unsigned int a;
    unsigned int i, j;

    gfpcircmatrix T;

    T = gfpcircm_init(trans->height, trans->width);
    gfpcircm_copy(T, *trans);

    a = trans->width;
    trans->width = trans->height;
    trans->height = a;

    for( i = 0 ; i < trans->height ; ++i )
    {
        for( j = 0 ; j < trans->width ; ++j )
        {
            gfpcirc_copy(&trans->data[i*trans->width + j], &T.data[j*T.width + i]);
        }
    }

    gfpcircm_destroy(T);

    return 1;
}

/**
 * gfpcircm_multiply
 * Multiplies two matrices and stores the result in the third, which
 * should be pre-allocated.
 * @params
 *  * dest : field matrix object to store the matrix product
 *  * left, right : field matrix object representing left- and right-
 *    hand-sides respectively.
 * @return
 *  * 1 if success, 0 otherwise
 */
int gfpcircm_multiply( gfpcircmatrix * dest, gfpcircmatrix left, gfpcircmatrix right )
{
    unsigned  int i, j, k;
    gfpcirc_element prod, sum, lsum, lhs, rhs;
    gfpcirc_element * data;

    prod = gfpcirc_init(1);
    sum = gfpcirc_init(1);
    lsum = gfpcirc_init(1);

    data = malloc(sizeof(gfpcirc_element) * dest->width * dest->height);
    for( i = 0 ; i < dest->height ; ++i )
    {
        for( j = 0 ; j < dest->width ; ++j )
        {
            data[i*dest->width + j] = gfpcirc_init(0);
        }
    }

    #ifdef DEBUG
        if( dest->height != left.height || dest->width != right.width || left.width != right.height )
        {
            printf("in gfpcircm_multiply: trying to multiply matrices with unmatched dimensions: %ix%i * %ix%i = %ix%i\n", left.height, left.width, right.height, right.width, dest->height, dest->width);
            return 0;
        }
    #endif

    for( i = 0 ; i < left.height ; ++i )
    {
        for( j = 0 ; j < right.width ; ++j )
        {
            gfpcirc_zero(&sum);
            for( k = 0 ; k < left.width ; ++k )
            {
                gfpcirc_copy(&lsum, &sum);
                lhs = left.data[i*left.width + k];
                rhs = right.data[k*right.width + j];
                gfpcirc_multiply(&prod, lhs, rhs);
                gfpcirc_add(&sum, lsum, prod);
            }
            gfpcirc_copy(&data[i*dest->width + j], &sum);
        }
    }

    for( i = 0 ; i < dest->height ; ++i )
    {
        for( j = 0 ; j < dest->width ; ++j )
        {
            lhs = dest->data[i*dest->width + j];
            gfpcirc_destroy(lhs);
        }
    }
    free(dest->data);
    dest->data = data;

    gfpcirc_destroy(prod);
    gfpcirc_destroy(sum);
    gfpcirc_destroy(lsum);

    return 1;
}

/**
 * gfpcircm_multiply_transpose
 * Multiplies the left hand side matrix with the transpose of the
 * hand side matrix and stores the result in the third, which
 * should be pre-allocated.
 * @params
 *  * dest : field matrix object to store the matrix product
 *  * left, rightT : field matrix object representing left- and right-
 *    hand-sides respectively.
 * @return
 *  * 1 if success, 0 otherwise
 */
int gfpcircm_multiply_transpose( gfpcircmatrix * dest, gfpcircmatrix left, gfpcircmatrix rightT )
{
    unsigned  int i, j, k;
    gfpcirc_element prod, sum, lsum;
    gfpcirc_element * data;

    prod = gfpcirc_init(1);
    sum = gfpcirc_init(1);
    lsum = gfpcirc_init(1);

    #ifdef DEBUG
        if( dest->height != left.height || dest->width != rightT.height || left.width != rightT.width )
        {
            printf("in gfpcircm_multiply_transpose: trying to multiply matrices with unmatched dimensions: %ix%i * (%ix%i)^T = %ix%i\n", left.height, left.width, rightT.height, rightT.width, dest->height, dest->width);
            return 0;
        }
    #endif

    data = malloc(sizeof(gfpcirc_element) * dest->width * dest->height);
    for( i = 0 ; i < dest->height ; ++i )
    {
        for( j = 0 ; j < dest->width ; ++j )
        {
            data[i*dest->width + j] = gfpcirc_init(0);
        }
    }

    for( i = 0 ; i < left.height ; ++i )
    {
        for( j = 0 ; j < rightT.height ; ++j )
        {
            gfpcirc_zero(&sum);
            for( k = 0 ; k < left.width ; ++k )
            {
                gfpcirc_copy(&lsum, &sum);
                gfpcirc_multiply(&prod, left.data[i*left.width + k], rightT.data[j*rightT.width + k]);
                gfpcirc_add(&sum, lsum, prod);
            }
            gfpcirc_copy(&data[i*dest->width + j], &sum);
        }
    }

    for( i = 0 ; i < dest->height ; ++i )
    {
        for( j = 0 ; j < dest->width ; ++j )
        {
            gfpcirc_destroy(dest->data[i*dest->width + j]);
        }
    }
    free(dest->data);
    dest->data = data;

    gfpcirc_destroy(prod);
    gfpcirc_destroy(sum);
    gfpcirc_destroy(lsum);

    return 1;
}

/**
 * gfpcircm_transpose_multiply
 * Multiplies the transpose of the left hand side matrix with the
 * hand side matrix and stores the result in the third, which
 * should be pre-allocated.
 * @params
 *  * dest : field matrix object to store the matrix product
 *  * leftT, right : field matrix object representing left- and right-
 *    hand-sides respectively.
 * @return
 *  * 1 if success, 0 otherwise
 */
int gfpcircm_transpose_multiply( gfpcircmatrix * dest, gfpcircmatrix leftT, gfpcircmatrix right )
{
    unsigned  int i, j, k;
    gfpcirc_element prod, sum, lsum, rhs, lhs;
    gfpcirc_element * data;

    prod = gfpcirc_init(1);
    sum = gfpcirc_init(1);
    lsum = gfpcirc_init(1);

    #ifdef DEBUG
        if( dest->height != leftT.width || dest->width != right.width || leftT.height != right.height )
        {
            printf("in gfpcircm_transpose_multiply: trying to multiply matrices with unmatched dimensions: (%ix%i)^T * %ix%i = %ix%i\n", leftT.height, leftT.width, right.height, right.width, dest->height, dest->width);
            return 0;
        }
    #endif

    data = malloc(sizeof(gfpcirc_element)*dest->width*dest->height);
    for( i = 0 ; i < dest->height ; ++i )
    {
        for( j = 0 ; j < dest->width ; ++j )
        {
            data[i*dest->width + j] = gfpcirc_init(0);
        }
    }

    for( i = 0 ; i < leftT.width ; ++i )
    {
        for( j = 0 ; j < right.width ; ++j )
        {
            gfpcirc_zero(&sum);
            for( k = 0 ; k < leftT.height ; ++k )
            {
                gfpcirc_copy(&lsum, &sum);
                lhs = leftT.data[k*leftT.width + i];
                rhs = right.data[k*right.width + j];
                gfpcirc_multiply(&prod, lhs, rhs);
                gfpcirc_add(&sum, lsum, prod);
            }
            gfpcirc_copy(&data[i*dest->width + j], &sum);
        }
    }

    for( i = 0 ; i < dest->height ; ++i )
    {
        for( j = 0 ; j < dest->width ; ++j )
        {
            lhs = dest->data[i*dest->width + j];
            gfpcirc_destroy(lhs);
        }
    }
    free(dest->data);
    dest->data = data;

    gfpcirc_destroy(prod);
    gfpcirc_destroy(sum);
    gfpcirc_destroy(lsum);

    return 1;
}

/**
 * gfpcircm_multiply_constant
 * Multiply the matrix with a constant.
 * @return
 *  * 1 if success
 */
int gfpcircm_multiply_constant( gfpcircmatrix dest, gfpcircmatrix source, gfpcirc_element constant )
{
    unsigned  int i, j;
    gfpcirc_element lhs, rhs;

#ifdef DEBUG
    if( dest.width != source.width || dest.height != source.height )
    {
        printf("gfpcircm_multiply_constant: cannot multiply matrix with constant because dimensions of destination do not match those of source! %ix%i <- %ix%i\n", dest.height, dest.width, source.height, source.width);
        return 0;
    }
#endif

    for( i = 0 ; i < dest.height ; ++i )
    {
        for( j = 0 ; j < dest.width ; ++j )
        {
            gfpcirc_multiply(&dest.data[i*dest.width + j], source.data[i*dest.width + j], constant);
        }
    }
    return 1;
}

/**
 * gfpcircm_add
 * Add one matrix to another and store the result in a third. The
 * third matrix must be preallocated.
 * @params
 *  * dest : the matrix object to store the result in
 *  * left, right : the matrix objects to add together; they must 
 *    have the same dimensions!
 * @return
 *  * 1 if success, 0 otherwise
 */
int gfpcircm_add( gfpcircmatrix dest, gfpcircmatrix left, gfpcircmatrix right )
{
    unsigned  int i, j;
    gfpcirc_element lhs, rhs;

    #ifdef DEBUG
        if( dest.width != left.width || left.width != right.width || dest.height != left.height || left.height != right.height )
        {
            printf("in gfpcircm_add: trying to add matrices of incompatible dimensions! %ix%i + %ix%i = %ix%i\n", left.height, left.width, right.height, right.width, dest.height, dest.width);
            return 0;
        }
    #endif
    
    for( i = 0 ; i < dest.height ; ++i )
    {
        for( j = 0 ; j < dest.width ; ++j )
        {
            lhs = left.data[i*left.width + j];
            rhs = right.data[i*right.width + j];
            gfpcirc_add(&dest.data[i*dest.width + j], lhs, rhs);
        }
    }

    return 1;
}

/**
 * gfpcircm_weighted_sum
 * Compute the weighted sum of two matrices, and store the result in
 * a third one. This third matrix must be pre-allocated.
 * @params:
 *  * dest : the matrix object to store the result into
 *  * left_constant, right_constant : gfpcirc_elements that represent
 *    the field elements to weight the left and right matrices with
 *  * left_matrix, right_matrix : the two matrix objects to add
 *    together.
 * @return
 * 1 if success, 0 otherwise
 */
int gfpcircm_weighted_sum( gfpcircmatrix dest, gfpcirc_element left_constant, gfpcircmatrix left_matrix, gfpcirc_element right_constant, gfpcircmatrix right_matrix )
{
    unsigned  int i, j;

    gfpcirc_element lhs, rhs;

    lhs = gfpcirc_init(1);
    rhs = gfpcirc_init(1);

    #ifdef DEBUG
        if( dest.width != left_matrix.width || left_matrix.width != right_matrix.width || dest.height != left_matrix.height || left_matrix.height != right_matrix.height )
        {
            printf("in gfpcircm_weighted_sum: trying to add matrices of incompatible dimensions! %ix%i + %ix%i = %ix%i\n", left_matrix.height, left_matrix.width, right_matrix.height, right_matrix.width, dest.height, dest.width);
            return 0;
        }
    #endif
    
    for( i = 0 ; i < dest.height ; ++i )
    {
        for( j = 0 ; j < dest.width ; ++j )
        {
            gfpcirc_multiply(&lhs, left_matrix.data[i*left_matrix.width + j], left_constant);
            gfpcirc_multiply(&rhs, right_matrix.data[i*right_matrix.width + j], right_constant);
            gfpcirc_add(&dest.data[i*dest.width + j], lhs, rhs);
        }
    }

    gfpcirc_destroy(lhs);
    gfpcirc_destroy(rhs);

    return 1;
}

/**
 * gfpcircm_flip
 * Flip the elements of the matrix according to the following rule:
 * all coefficients are reversed, except the last.
 */
int gfpcircm_flip( gfpcircmatrix mat )
{
    unsigned int i, j;
    int k;
    gfpcirc_element elm;
    for( i = 0 ; i < mat.height ; ++i )
    {
        for( j = 0 ; j < mat.width ; ++j )
        {
            for( k = 0 ; k < DEGREE_OF_CIRCULANCY-1 ; ++k )
            {
                elm.data[k] = mat.data[i*mat.width + j].data[k];
            }
            for( k = 0 ; k < DEGREE_OF_CIRCULANCY-1 ; ++k )
            {
                mat.data[i*mat.width + j].data[k] = elm.data[DEGREE_OF_CIRCULANCY-1-k];
            }
        }
    }
}

/**
 * gfpcircm_stack
 * Stacks one matrix on top of another, and stores the result in the third
 * @params
 *  * mat : matrix object to store the result into
 *  * top, bottom : matrix objects to stack
 * @return
 *  * 1 if success, 0 otherwise
 */
int gfpcircm_stack( gfpcircmatrix mat, gfpcircmatrix top, gfpcircmatrix bottom )
{
    unsigned  int i, j;

    #ifdef DEBUG
        if( mat.width != top.width || top.width != bottom.width || mat.height != top.height + bottom.height )
        {
            printf("in gfpcircm_stack: cannot stack matrices of conflicting dimensions! %ix%i stack %ix%i = %ix%i\n", top.height, top.width, bottom.height, bottom.width, mat.height, mat.width);
            return 0;
        }
    #endif
    for( i = 0 ; i < top.height ; ++i )
    {
        for( j = 0 ; j < top.width ; ++j )
        {
            gfpcirc_copy(&mat.data[i*top.width + j], &top.data[i*top.width + j]);
        }
    }
    for( i = 0 ; i < bottom.height ; ++i )
    {
        for( j = 0 ; j < bottom.width ; ++j )
        {
            gfpcirc_copy(&mat.data[(i+top.height)*mat.width + j], &bottom.data[i*bottom.width + j]);
        }
    }

    return 1;
}

/**
 * gfpcircm_cat
 * Concatenates one matrix to another, and stores the result in a
 * third matrix.
 * @params
 *  * res : matrix object that contains the result
 *  * left, right : matrix objects to be stacked on the left, and
 *    right, respectively
 * @return
 *  * 1 if success, 0 otherwise
 */
int gfpcircm_cat( gfpcircmatrix res, gfpcircmatrix left, gfpcircmatrix right )
{
    unsigned  int i, j;

    #ifdef DEBUG
        if( res.height != left.height || left.height != right.height || res.width != left.width + right.width )
        {
            printf("in gfpcircm_cat: cannot concatenate two matrices of conflicting dimensions! %ix%i cat %ix%i = %ix%i\n", left.height, left.width, right.height, right.width, res.height, res.width);
            return 0;
        }
    #endif

    for( i = 0 ; i < res.height ; ++i )
    {
        for( j = 0 ; j < left.width ; ++j )
        {
            gfpcirc_copy(&res.data[i*res.width + j], &left.data[i*left.width + j]);
        }
        for( j = 0 ; j < right.width ; ++j )
        {
            gfpcirc_copy(&res.data[i*res.width + left.width + j], &right.data[i*right.width + j]);
        }
    }
    return 1;
}

/**
 * gfpcircm_slice
 * Slice a submatrix out of another matrix.
 * @params:
 *  * dest : the matrix to store the result in; the height and width
 *    of this matrix determine the area of elements to be copied
 *  * source : the matrix where the slice comes from
 *  * row_start, col_start : the indices of the row and column at
 *    which to start the slice
 * @return
 *  * 1 if success, 0 otherwise
 */
int gfpcircm_slice( gfpcircmatrix dest, gfpcircmatrix source, unsigned  int row_start, unsigned  int col_start )
{
    unsigned  int i, j;
    #ifdef DEBUG
        if( source.width < col_start + dest.width || source.height < row_start + dest.height )
        {
            printf("in gfpcircm_slice: cannot grab slice because slice size exceeds bounds! slicing %ix%i submatrix starting at (%i,%i) from %ix%i matrix\n", dest.height, dest.width, row_start, col_start, source.height, source.width);
            return 0;
        }
    #endif
    for( i = 0 ; i < dest.height ; ++i )
    {
        for( j = 0 ; j < dest.width ; ++j )
        {
            gfpcirc_copy(&dest.data[i*dest.width + j], &source.data[(i+row_start)*source.width + col_start + j]);
        }
    }
    return 0;
}

/**
 * gfpcircm_print
 * Use printf to print the matrix to stdout.
 */
int gfpcircm_print( gfpcircmatrix mat )
{
    unsigned int i, j;
    printf("[");
    for( i = 0 ; i < mat.height ; ++i )
    {
        printf("[");
        for( j = 0 ; j < mat.width ; ++j )
        {
            gfpcirc_print(&mat.data[i*mat.width+j]);
            if( j < mat.width - 1 )
            {
                printf(",");
                printf(" ");
            }
        }
        printf("]");
        if( i < mat.height - 1 )
        {
            printf(",");
            printf("\n");
        }
    }
    printf("]\n");
}

/**
 * gfpcircm_print_transpose
 * Use printf to print the transpose of the matrix to stdout.
 */
int gfpcircm_print_transpose( gfpcircmatrix mat )
{
    gfpcircmatrix temp;
    temp = gfpcircm_clone(mat);
    gfpcircm_transpose(&temp);
    gfpcircm_print(temp);
    gfpcircm_destroy(temp);
    return 1;
}

