#include "gfpem.h"

#include <stdlib.h>
#include <stdio.h>

/**
 * gfpem
 * Create gfpematrix object with given buffer. Usefule for reusing
 * the same data line.
 */
gfpematrix gfpem( unsigned  int height, unsigned  int width, gfpe_element* pdata )
{
    int i, j;
    gfpematrix mat;
    mat.height = height;
    mat.width = width;
    mat.data = pdata;
    return mat;
}

/**
 * gfpem_init
 * Create a field matrix object and allocate space for it.
 */
gfpematrix gfpem_init( unsigned  int height, unsigned  int width )
{
    int i, j;
    gfpematrix mat;
    mat.width = width;
    mat.height = height;
    mat.data = malloc(width*height*sizeof(gfpe_element));
    /*
    printf("created gfpem object with data member set to memory address %#010x\n", mat.data);
    */
    return mat;
}

/**
 * gfpem_destroy
 * Deallocates space allocated to a field matrix object. Call this
 * function before closing the scope where the field matrix object
 * was initialized.
 * @return
 *  * 1 if success
 */
int gfpem_destroy( gfpematrix fm )
{
    int i, j;
    free(fm.data);
    fm.width = 0;
    fm.height = 0;
    return 1;
}

/**
 * gfpem_copy
 * Copy the contents of one matrix to another. Does not allocate
 * memory for the new object; you must do that yourself! (Or use
 * gfpem_clone instead.)
 * @promise
 *  * dest.width >= source.width
 *  * dest.height >= source.height
 * @return
 *  * 1 if success
 */
int gfpem_copy( gfpematrix dest, gfpematrix source )
{
    unsigned int i, j;
    for( i = 0 ; i < source.height ; ++i )
    {
        for( j = 0 ; j < source.width ; ++j )
        {
            gfpe_copy(&dest.data[i*dest.width + j], source.data[i*source.width + j]);
        }
    }

    return 1;
}

/**
 * gfpem_clone
 * Copy one matrix into a new object. Don't forget to call
 * gfpem_destroy at the end of scope!
 */
gfpematrix gfpem_clone( gfpematrix source )
{
    gfpematrix mat;
    mat = gfpem_init(source.height, source.width);
    gfpem_copy(mat, source);
    return mat;
}

/**
 * gfpem_eye
 * Set a matrix to the identity matrix.
 * @param
 *  * eye : matrix object to set to identity
 * @returns
 *  * 1 if success, 0 otherwise
 */
int gfpem_eye( gfpematrix eye )
{
    unsigned  int i, j;
    for( i = 0 ; i < eye.height ; ++i )
    {
        for( j = 0 ; j < eye.width ; ++j )
        {
            gfpe_zero(&eye.data[i*eye.width + j]);
        }
        gfpe_one(&eye.data[i*eye.width + i]);
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
int gfpem_is_eye( gfpematrix eye )
{
    unsigned int i, j;
    int b = 1;
    gfpe_element one, zero;
    
    gfpe_one(&one);
    gfpe_zero(&zero);

    for( i = 0 ; i < eye.height ; ++i )
    {
        for( j = 0 ; j < eye.width ; ++j )
        {
            if( i == j )
            {
                b = b & gfpe_compare(eye.data[i*eye.width + j], one);
            }
            else
            {
                b = b & gfpe_compare(eye.data[i*eye.width + j], zero);
            }
        }
    }

    return b;
}

/**
 * gfpem_equals
 * Test two matrices for equality.
 * @return
 *  * 1 if equal, 0 otherwise
 */
int gfpem_equals( gfpematrix lhs, gfpematrix rhs )
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
            b = b & gfpe_compare(lhs.data[lhs.width*i + j], rhs.data[rhs.width*i + j]);
        }
    }

    return b;
}

/**
 * gfpem_zero
 * Sets a matrix to all zeros.
 * @param
 *  * zero : matrix object to set to zero
 * @returns
 *  * 1 if success
 */
int gfpem_zeros( gfpematrix zero )
{
    unsigned  int i, j;
    for( i = 0 ; i < zero.height ; ++i )
    {
        for( j = 0 ; j < zero.width ; ++j )
        {
            gfpe_zero(&zero.data[i*zero.width + j]);
        }
    }
    return 1;
}

/**
 * gfpem_random
 * Put random values into the matrix.
 * @params
 *  * random : matrix objects with data already allocated and whose
 *    elements are to be assigned random values
 *  * randomness : pointer to large-enough string of random bytes
 *    "large enough" means n*n*sizeof(unsigned int)
 * @result
 *  * random <-$- matrix_space(random.height, random.width)
 */
int gfpem_random( gfpematrix random, unsigned char * randomness )
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
            gfpe_random(&random.data[i*random.width + j], randomness + l*num_limbs*sizeof(unsigned long int));
            l = l + 1;
        }
    }

    return 1;
}

/**
 * gfpem_random_upper_triangular
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
int gfpem_random_upper_triangular( gfpematrix random, unsigned char * randomness )
{
    unsigned  int i, j;
    unsigned int l;

    l = 0;
    for( i = 0 ; i < random.height ; ++i )
    {
        for( j = 0 ; j < i ; ++j )
        {
            gfpe_random(&random.data[i*random.width + j], &randomness[(l++)*(GFP_NUMBYTES+1)]);
        }
        gfpe_one(&random.data[i*random.width + i]);
        for( j = i+1 ; j < random.width ; ++j )
        {
            gfpe_zero(&random.data[i*random.width + j]);
        }
    }

    return 1;
}

/**
 * gfpem_transpose_square
 * Perform a matrix transposition in situ.
 */
int gfpem_transpose( gfpematrix * trans )
{
    unsigned int a;
    unsigned int i, j;

    gfpematrix T;

    T = gfpem_init(trans->height, trans->width);
    gfpem_copy(T, *trans);

    a = trans->width;
    trans->width = trans->height;
    trans->height = a;

    for( i = 0 ; i < trans->height ; ++i )
    {
        for( j = 0 ; j < trans->width ; ++j )
        {
            gfpe_copy(&trans->data[i*trans->width + j], T.data[j*T.width + i]);
        }
    }

    gfpem_destroy(T);

    return 1;
}

/**
 * gfpem_multiply
 * Multiplies two matrices and stores the result in the third, which
 * should be pre-allocated.
 * @params
 *  * dest : field matrix object to store the matrix product
 *  * left, right : field matrix object representing left- and right-
 *    hand-sides respectively.
 * @return
 *  * 1 if success, 0 otherwise
 * NOTE. Modular reduction is applied only once, after computing the
 * inner product between left row and right column and storing the
 * intermediate result in an int. If the int is four bytes and the
 * field element takes up a full byte, then in the worst case every
 * product requires two bytes and every 2^k additions requires k
 * extra bits. Since we have 16 bits to spare we have at most 2^16
 * additions which is anyway the max. height and width of matrices
 * that can be stored in a  unsigned int.
 */
int gfpem_multiply( gfpematrix * dest, gfpematrix left, gfpematrix right )
{
    unsigned  int i, j, k;
    gfpe_element prod, sum, lsum;
    gfpe_element * data;

    data = malloc(sizeof(gfpe_element) * dest->width * dest->height);

    #ifdef DEBUG
        if( dest->height != left.height || dest->width != right.width || left.width != right.height )
        {
            printf("in gfpem_multiply: trying to multiply matrices with unmatched dimensions: %ix%i * %ix%i = %ix%i\n", left.height, left.width, right.height, right.width, dest->height, dest->width);
            return 0;
        }
    #endif

    for( i = 0 ; i < left.height ; ++i )
    {
        for( j = 0 ; j < right.width ; ++j )
        {
            gfpe_zero(&sum);
            for( k = 0 ; k < left.width ; ++k )
            {
                gfpe_copy(&lsum, sum);
                gfpe_multiply(&prod, left.data[i*left.width + k], right.data[k*right.width + j]);
                gfpe_add(&sum, lsum, prod);
            }
            gfpe_copy(&data[i*dest->width + j], sum);
        }
    }

    free(dest->data);
    dest->data = data;
    return 1;
}

/**
 * gfpem_multiply_transpose
 * Multiplies the left hand side matrix with the transpose of the
 * hand side matrix and stores the result in the third, which
 * should be pre-allocated.
 * @params
 *  * dest : field matrix object to store the matrix product
 *  * left, rightT : field matrix object representing left- and right-
 *    hand-sides respectively.
 * @return
 *  * 1 if success, 0 otherwise
 * NOTE. Modular reduction is applied only once, after computing the
 * inner product between left row and right column and storing the
 * intermediate result in an int. If the int is four bytes and the
 * field element takes up a full byte, then in the worst case every
 * product requires two bytes and every 2^k additions requires k
 * extra bits. Since we have 16 bits to spare we have at most 2^16
 * additions which is anyway the max. height and width of matrices
 * that can be stored in a  unsigned int.
 */
int gfpem_multiply_transpose( gfpematrix * dest, gfpematrix left, gfpematrix rightT )
{
    unsigned  int i, j, k;
    gfpe_element prod, sum, lsum;
    gfpe_element * data;

    #ifdef DEBUG
        if( dest->height != left.height || dest->width != rightT.height || left.width != rightT.width )
        {
            printf("in gfpem_multiply_transpose: trying to multiply matrices with unmatched dimensions: %ix%i * (%ix%i)^T = %ix%i\n", left.height, left.width, rightT.height, rightT.width, dest->height, dest->width);
            return 0;
        }
    #endif

    data = malloc(sizeof(gfpe_element) * dest->width * dest->height);

    for( i = 0 ; i < left.height ; ++i )
    {
        for( j = 0 ; j < rightT.height ; ++j )
        {
            gfpe_zero(&sum);
            for( k = 0 ; k < left.width ; ++k )
            {
                gfpe_copy(&lsum, sum);
                gfpe_multiply(&prod, left.data[i*left.width + k], rightT.data[j*rightT.width + k]);
                gfpe_add(&sum, lsum, prod);
            }
            gfpe_copy(&data[i*dest->width + j], sum);
        }
    }

    free(dest->data);
    dest->data = data;
    return 1;
}

/**
 * gfpem_transpose_multiply
 * Multiplies the transpose of the left hand side matrix with the
 * hand side matrix and stores the result in the third, which
 * should be pre-allocated.
 * @params
 *  * dest : field matrix object to store the matrix product
 *  * leftT, right : field matrix object representing left- and right-
 *    hand-sides respectively.
 * @return
 *  * 1 if success, 0 otherwise
 * NOTE. Modular reduction is applied only once, after computing the
 * inner product between left row and right column and storing the
 * intermediate result in an int. If the int is four bytes and the
 * field element takes up a full byte, then in the worst case every
 * product requires two bytes and every 2^k additions requires k
 * extra bits. Since we have 16 bits to spare we have at most 2^16
 * additions which is anyway the max. height and width of matrices
 * that can be stored in a  unsigned int.
 */
int gfpem_transpose_multiply( gfpematrix * dest, gfpematrix leftT, gfpematrix right )
{
    unsigned  int i, j, k;
    gfpe_element prod, sum, lsum;
    gfpe_element * data;

    #ifdef DEBUG
        if( dest->height != leftT.width || dest->width != right.width || leftT.height != right.height )
        {
            printf("in gfpem_transpose_multiply: trying to multiply matrices with unmatched dimensions: (%ix%i)^T * %ix%i = %ix%i\n", leftT.height, leftT.width, right.height, right.width, dest->height, dest->width);
            return 0;
        }
    #endif

    data = malloc(sizeof(gfpe_element)*dest->width*dest->height);

    for( i = 0 ; i < leftT.width ; ++i )
    {
        for( j = 0 ; j < right.width ; ++j )
        {
            gfpe_zero(&sum);
            for( k = 0 ; k < leftT.height ; ++k )
            {
                gfpe_copy(&lsum, sum);
                gfpe_multiply(&prod, leftT.data[k*leftT.width + i], right.data[k*right.width + j]);
                gfpe_add(&sum, lsum, prod);
            }
            gfpe_copy(&data[i*dest->width + j], sum);
        }
    }
    free(dest->data);
    dest->data = data;

    return 1;
}

/**
 * gfpem_multiply_constant
 * Multiply the matrix with a constant.
 * @return
 *  * 1 if success
 */
int gfpem_multiply_constant( gfpematrix dest, gfpematrix source, gfpe_element constant )
{
    unsigned  int i, j;
    gfpe_element lhs, rhs;

#ifdef DEBUG
    if( dest.width != source.width || dest.height != source.height )
    {
        printf("gfpem_multiply_constant: cannot multiply matrix with constant because dimensions of destination do not match those of source! %ix%i <- %ix%i\n", dest.height, dest.width, source.height, source.width);
        return 0;
    }
#endif

    for( i = 0 ; i < dest.height ; ++i )
    {
        for( j = 0 ; j < dest.width ; ++j )
        {
            gfpe_multiply(&dest.data[i*dest.width + j], source.data[i*dest.width + j], constant);
        }
    }
    return 1;
}

/**
 * gfpem_add
 * Add one matrix to another and store the result in a third. The
 * third matrix must be preallocated.
 * @params
 *  * dest : the matrix object to store the result in
 *  * left, right : the matrix objects to add together; they must 
 *    have the same dimensions!
 * @return
 *  * 1 if success, 0 otherwise
 */
int gfpem_add( gfpematrix dest, gfpematrix left, gfpematrix right )
{
    unsigned  int i, j;

    #ifdef DEBUG
        if( dest.width != left.width || left.width != right.width || dest.height != left.height || left.height != right.height )
        {
            printf("in gfpem_add: trying to add matrices of incompatible dimensions! %ix%i + %ix%i = %ix%i\n", left.height, left.width, right.height, right.width, dest.height, dest.width);
            return 0;
        }
    #endif
    
    for( i = 0 ; i < dest.height ; ++i )
    {
        for( j = 0 ; j < dest.width ; ++j )
        {
            gfpe_add(&dest.data[i*dest.width + j], left.data[i*left.width + j], right.data[i*right.width + j]);
        }
    }

    return 1;
}

/**
 * gfpem_weighted_sum
 * Compute the weighted sum of two matrices, and store the result in
 * a third one. This third matrix must be pre-allocated.
 * @params:
 *  * dest : the matrix object to store the result into
 *  * left_constant, right_constant : gfpe_elements that represent
 *    the field elements to weight the left and right matrices with
 *  * left_matrix, right_matrix : the two matrix objects to add
 *    together.
 * @return
 * 1 if success, 0 otherwise
 */
int gfpem_weighted_sum( gfpematrix dest, gfpe_element left_constant, gfpematrix left_matrix, gfpe_element right_constant, gfpematrix right_matrix )
{
    unsigned  int i, j;

    gfpe_element lhs, rhs;

    #ifdef DEBUG
        if( dest.width != left_matrix.width || left_matrix.width != right_matrix.width || dest.height != left_matrix.height || left_matrix.height != right_matrix.height )
        {
            printf("in gfpem_weighted_sum: trying to add matrices of incompatible dimensions! %ix%i + %ix%i = %ix%i\n", left_matrix.height, left_matrix.width, right_matrix.height, right_matrix.width, dest.height, dest.width);
            return 0;
        }
    #endif
    
    for( i = 0 ; i < dest.height ; ++i )
    {
        for( j = 0 ; j < dest.width ; ++j )
        {
            gfpe_multiply(&lhs, left_matrix.data[i*left_matrix.width + j], left_constant);
            gfpe_multiply(&rhs, right_matrix.data[i*right_matrix.width + j], right_constant);
            gfpe_add(&dest.data[i*dest.width + j], lhs, rhs);
        }
    }

    return 1;
}

/**
 * gfpem_rowop
 * Perform a row operation on the given matrix, i.e., add one row,
 * weighted by a constant, to another.
 * @params
 *  * mat : the matrix object to operate on
 *  * destrow, sourcerow : unsigned  ints representing indices of the rows to operate on
 *  * constant : gfpe_element representing the right constant
 *  * offset : unsigned  int, represents the number of zeros to skip before applying the row operation
 * @returns
 *  * 1 if success
 */
int gfpem_rowop( gfpematrix mat, unsigned  int destrow, unsigned  int sourcerow, gfpe_element constant, unsigned  int offset )
{
    unsigned  int j;
    gfpe_element prod, sum;

    for( j = offset ; j < mat.width ; ++j )
    {
        gfpe_multiply(&prod, mat.data[sourcerow*mat.width + j], constant);
        gfpe_add(&sum, mat.data[destrow*mat.width + j], prod);
        gfpe_copy(&mat.data[destrow*mat.width + j], sum);
    }

    return 1;
}

/**
 * gfpem_scalerow
 * Scales a single row in the matrix with a given constant.
 * @params
 *  * mat : the matrix object to operate on
 *  * rowidx : index of the row to scale
 *  * constant : gfpe_element -- the field element to multiply the
 *    row with
 * @returns
 *  * 1 if success, 0 otherwise
 */
int gfpem_scalerow( gfpematrix mat, unsigned  int rowidx, gfpe_element constant )
{
    unsigned  int j;
    gfpe_element temp;

    for( j = 0 ; j < mat.width ; ++j )
    {
        gfpe_multiply(&temp, mat.data[rowidx*mat.width + j], constant);
        gfpe_copy(&mat.data[rowidx*mat.width + j], temp);
    }

    return 1;
}

/**
 * gfpem_fliprows
 * Flip two rows in the given matrix.
 * @params
 *  * mat : the matrix object to operate on
 *  * destrow, sourcerow : the indices of the rows to flip
 * @return
 *  * 1 if success
 */
int gfpem_fliprows( gfpematrix mat, unsigned  int destrow, unsigned  int sourcerow )
{
    unsigned  int j;
    gfpe_element a;

    for( j = 0 ; j < mat.width ; ++j )
    {
        gfpe_copy(&a, mat.data[destrow*mat.width + j]);
        gfpe_copy(&mat.data[destrow*mat.width + j], mat.data[sourcerow*mat.width + j]);
        gfpe_copy(&mat.data[sourcerow*mat.width + j], a);
    }

    return 1;
}

/**
 * gfpem_redech
 * Reduce the given matrix to reduced row echelon form using row
 * operations.
 * @params
 *  * mat : the matrix object to work on
 * @return
 *  1 if success
 */
int gfpem_redech( gfpematrix mat )
{
    unsigned  int col, row, i;
    gfpe_element inv;
    gfpe_element prod;
    gfpe_element diff;

    row = 0;
    for( col = 0 ; col < mat.width ; ++col )
    {
        for( i = row ; i < mat.height ; ++i )
        {
            if( gfpe_is_zero(mat.data[i*mat.width + col]) != 1 )
            {
                if( i != row )
                {
                    gfpem_fliprows(mat, i, row);
                }
                break;
            }
        }

        if( i == mat.height )
        {
            continue;
        }
        gfpe_inverse(&inv, mat.data[row*mat.width + col]);
       
        if( gfpe_is_one(inv) != 1 )
        {
            gfpem_scalerow(mat, row, inv);
        }

        for( i = 0 ; i < mat.height ; ++i )
        {
            if( i == row )
            {
                continue;
            }
            gfpe_negate(&diff, mat.data[i*mat.width + col]);
            gfpem_rowop(mat, i, row, diff, col);
        }

        row = row + 1;

        if( row == mat.height )
        {
            break;
        }
    }

    return 1;
}

/**
 * gfpem_solve
 * Solve a matrix equation of the form Ax = b for x up to a term in
 * the kernel of A. This routine also initializes a kernel matrix,
 * whose rows form a basis for the kernel of A.
 * @params
 *  * coeffs : a mxn gfpematrix object representing the coefficient
 *    matrix A
 *  * target : a mx1 gfpematrix object representing the b vector
 *  * solution : a nx1 gfpematrix object to store one solution into
 *  * kernel : an uninitialized gfpematrix object whose columns will
 *    span the kernel of A
 * @post
 *  * for all kernel.width x 1 vectors "random" holds:
 *          coeffs * (solution + kernel * random) = target
 * @return
 *  * 1 if a solution exists, 0 otherwise
 */
int gfpem_solve( gfpematrix coeffs, gfpematrix target, gfpematrix solution, gfpematrix * kernel )
{
    /* declare variables for echelon reduction */
    unsigned  int col, row, i, j;
    gfpe_element inv, zero, one, minusone, neg;
    gfpematrix mat;

    /* declare variables for pivot tracking */
    unsigned  int *pivots;
    unsigned  int *npivots;
    unsigned  int num_pivots;
    unsigned  int num_npivots;
    int have_solution;

    gfpe_zero(&zero);
    gfpe_one(&one);
    gfpe_negate(&minusone, one);

    /* initialize variables for pivot tracking */
    num_pivots = 0;
    num_npivots = 0;
    pivots = malloc(sizeof(unsigned  int) * (coeffs.width+1));
    npivots = malloc(sizeof(unsigned  int) * (coeffs.width+1));

    /* initialize mat and copy coeffs and target to it */
    mat = gfpem_init(coeffs.height, coeffs.width+1);
    /*gfpem_copy(mat, coeffs);*/
    for( i = 0 ; i < mat.height ; ++i )
    {
        for( j = 0 ; j < coeffs.width ; ++j )
        {
            gfpe_copy(&mat.data[i*mat.width + j], coeffs.data[i*coeffs.width + j]);
        }
    }
    for( i = 0 ; i < mat.height ; ++i )
    {
        gfpe_copy(&mat.data[i*mat.width + mat.width - 1], target.data[i*target.width + 0]);
    }

    /* perform row echelon reduction */
    row = 0;
    for( col = 0 ; col < mat.width ; ++col )
    {
        for( i = row ; i < mat.height ; ++i )
        {
            /* if the leading element is different from zero, use it as pivot element */
            if( gfpe_compare(mat.data[i*mat.width + col], zero) != 1 )
            {
                if( i != row )
                {
                    gfpem_fliprows(mat, i, row);
                }
                break;
            }
        }

        if( i == mat.height )
        {
            /* non pivot */
            npivots[num_npivots++] = col;
            continue;
        }

        /* pivot */
        pivots[num_pivots++] = col;

        gfpe_inverse(&inv, mat.data[row*mat.width + col]);
        
        /* rescale row if necessary */
        if( gfpe_compare(inv, one) != 1 )
        {
            gfpem_scalerow(mat, row, inv);
        }

        for( i = 0 ; i < mat.height ; ++i )
        {
            if( i == row )
            {
                continue;
            }
            gfpe_negate(&neg, mat.data[i*mat.width + col]);
            gfpem_rowop(mat, i, row, neg, col);
        }

        row = row + 1;

        if( row == mat.height )
        {
            for( i = col+1 ; i < mat.width ; ++i )
            {
                npivots[num_npivots++] = i;
            }
            break;
        }

    }

    /* read out solution if the system is consistent */
    have_solution = (pivots[num_pivots-1] != mat.width-1);
    for( i = 0 ; i < mat.width-1 ; ++i )
    {
        gfpe_zero(&solution.data[i*solution.width]);
    }
    if( have_solution == 1 )
    {
        for( i = 0 ; i < num_pivots ; ++i )
        {
            gfpe_copy(&solution.data[pivots[i]*solution.width], mat.data[i*mat.width + mat.width - 1]);
        }
    }

    /* read out kernel, if it exists */
    if( num_npivots > 1 )
    {
        *kernel = gfpem_init(mat.width-1, num_npivots-1);
        gfpem_zeros(*kernel);
        for( j = 0 ; j < num_npivots-1 ; ++j )
        {
            gfpe_copy(&(kernel->data[npivots[j]*kernel->width + j]), minusone);
            for( i = 0 ; i < num_pivots && pivots[i] < npivots[j] ; ++i )
            {
                gfpe_copy(&(kernel->data[pivots[i]*kernel->width + j]), mat.data[i*mat.width + npivots[j]]);
            }
        }
    }
    else
    {
        kernel->width = 0;
    }

    /* free allocated memory */
    gfpem_destroy(mat);
    free(pivots);
    free(npivots);

    return have_solution;
}

/**
 * gfpem_inspan
 * Decide if a vector is in the column span of a matrix.
 * @params:
 *  * vec : nx1 matrix
 *  * mat : nxm matrix
 * @returns:
 *  * 1 if vec is in colspan(mat); 0 otherwise
 */
int gfpem_inspan( gfpematrix vec, gfpematrix mat )
{
    gfpematrix solution, kernel;
    int i, j;
    int success;

    #ifdef DEBUG
        if( mat.height != vec.height || vec.width != 1 )
        {
            printf("cannot decide if %ix%i vector is in colspan of %ix%i matrix because of dimension mismatch.\n", vec.height, vec.width, mat.height, mat.width);
            return 0;
        }
    #endif

    /* try and solve the system; if the system is consistent, the vector
     * lies in the span of the matrix */
    solution = gfpem_init(mat.width, 1);
    success = gfpem_solve(mat, vec, solution, &kernel);
    gfpem_destroy(solution);
    if( kernel.width > 0 )
        gfpem_destroy(kernel);

    return success;
}

/**
 * gfpem_stack
 * Stacks one matrix on top of another, and stores the result in the third
 * @params
 *  * mat : matrix object to store the result into
 *  * top, bottom : matrix objects to stack
 * @return
 *  * 1 if success, 0 otherwise
 */
int gfpem_stack( gfpematrix mat, gfpematrix top, gfpematrix bottom )
{
    unsigned  int i, j;

    #ifdef DEBUG
        if( mat.width != top.width || top.width != bottom.width || mat.height != top.height + bottom.height )
        {
            printf("in gfpem_stack: cannot stack matrices of conflicting dimensions! %ix%i stack %ix%i = %ix%i\n", top.height, top.width, bottom.height, bottom.width, mat.height, mat.width);
            return 0;
        }
    #endif
    for( i = 0 ; i < top.height ; ++i )
    {
        for( j = 0 ; j < top.width ; ++j )
        {
            gfpe_copy(&mat.data[i*mat.width + j], top.data[i*top.width + j]);
        }
    }
    for( i = 0 ; i < bottom.height ; ++i )
    {
        for( j = 0 ; j < bottom.width ; ++j )
        {
            gfpe_copy(&mat.data[(i+top.height)*mat.width + j], bottom.data[i*bottom.width + j]);
        }
    }

    return 1;
}

/**
 * gfpem_cat
 * Concatenates one matrix to another, and stores the result in a
 * third matrix.
 * @params
 *  * res : matrix object that contains the result
 *  * left, right : matrix objects to be stacked on the left, and
 *    right, respectively
 * @return
 *  * 1 if success, 0 otherwise
 */
int gfpem_cat( gfpematrix res, gfpematrix left, gfpematrix right )
{
    unsigned  int i, j;

    #ifdef DEBUG
        if( res.height != left.height || left.height != right.height || res.width != left.width + right.width )
        {
            printf("in gfpem_cat: cannot concatenate two matrices of conflicting dimensions! %ix%i cat %ix%i = %ix%i\n", left.height, left.width, right.height, right.width, res.height, res.width);
            return 0;
        }
    #endif

    for( i = 0 ; i < res.height ; ++i )
    {
        for( j = 0 ; j < left.width ; ++j )
        {
            gfpe_copy(&res.data[i*res.width + j], left.data[i*left.width + j]);
        }
        for( j = 0 ; j < right.width ; ++j )
        {
            gfpe_copy(&res.data[i*res.width + left.width + j], right.data[i*right.width + j]);
        }
    }
    return 1;
}

/**
 * gfpem_slice
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
int gfpem_slice( gfpematrix dest, gfpematrix source, unsigned  int row_start, unsigned  int col_start )
{
    unsigned  int i, j;
    #ifdef DEBUG
        if( source.width < col_start + dest.width || source.height < row_start + dest.height )
        {
            printf("in gfpem_slice: cannot grab slice because slice size exceeds bounds! slicing %ix%i submatrix starting at (%i,%i) from %ix%i matrix\n", dest.height, dest.width, row_start, col_start, source.height, source.width);
            return 0;
        }
    #endif
    for( i = 0 ; i < dest.height ; ++i )
    {
        for( j = 0 ; j < dest.width ; ++j )
        {
            gfpe_copy(&dest.data[i*dest.width + j], source.data[(i+row_start)*source.width + col_start + j]);
        }
    }
    return 0;
}

/**
 * gfpem_inverse
 * Compute the matrix inverse of mat, store the result in inv.
 * @return
 *  * 1 if success
 */
int gfpem_inverse( gfpematrix inv, gfpematrix mat )
{
    unsigned int i, j;
    unsigned  int catwidth;
    int invertible;
    gfpematrix concat;
   
    catwidth = inv.width + mat.width;

    /* Set inv to the identity matrix. */
    for( i = 0 ; i < inv.height ; ++i )
    {
        for( j = 0 ; j < inv.width ; ++j )
        {
            gfpe_zero(&inv.data[i*inv.width + j]);
        }
        gfpe_one(&inv.data[i*inv.width + i]);
    }

    /* Concatenate mat with identity */
    concat = gfpem_init(mat.height, catwidth);
    gfpem_cat(concat, mat, inv);

    /* row-reduce concat to echelon form */
    gfpem_redech(concat);

    /* test if main diagonal has only ones, because otherwise the
     * matrix is not invertible */
    invertible = 1;
    for( i = 0 ; i < inv.height ; ++i )
    {
        invertible = invertible & gfpe_is_one(concat.data[i*concat.width + i]);
    }

    if( 0 == invertible )
    {
        gfpem_destroy(concat);
        return 0;
    }

    /* select rightmost square from concat */
    for( i = 0 ; i < inv.height ; ++i )
    {
        for( j = 0 ; j < inv.width ; ++j )
        {
            gfpe_copy(&inv.data[i*inv.width + j], concat.data[i*concat.width + mat.width + j]);
        }
    }

    /* free concat */
    gfpem_destroy(concat);

    return 1;
}

/**
 * gfpem_print
 * Use printf to print the matrix to stdout.
 */
int gfpem_print( gfpematrix mat )
{
    unsigned int i, j;
    printf("[");
    for( i = 0 ; i < mat.height ; ++i )
    {
        printf("[");
        for( j = 0 ; j < mat.width ; ++j )
        {
            gfpe_print(mat.data[i*mat.width+j]);
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
 * gfpem_print_transpose
 * Use printf to print the transpose of the matrix to stdout.
 */
int gfpem_print_transpose( gfpematrix mat )
{
    gfpematrix temp;
    temp = gfpem_clone(mat);
    gfpem_transpose(&temp);
    gfpem_print(temp);
    gfpem_destroy(temp);
    return 1;
}

