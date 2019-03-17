#include "bacuov.h"
#include "gfpem.h"
#include "gfpm.h"

#include <stdio.h>
#include <stdlib.h>

void FIPS202_SHAKE256(const unsigned char *input, unsigned int inputByteLen, unsigned char *output, int outputByteLen);
void FIPS202_SHA3_256(const unsigned char *input, unsigned int inputByteLen, unsigned char *output);

#ifndef BACUOV_PARAM_V
#error "parameter V undefined"
#endif

/**
 * bacuov_generate__S
 * Generate the secret matrix S, given
 * the seed for S.
 */
void bacuov_generate_S( gfpcircmatrix * data, unsigned char * seed_S )
{
    unsigned char randomness[SECURITY_LEVEL/4 + 4];
    int i, j, k;
    int type;
    gfpcirc_element elm;
    unsigned char buffer[GFP_NUMBYTES * DEGREE_OF_CIRCULANCY * BACUOV_PARAM_V*BACUOV_PARAM_O];

    // copy over randomness
    for( i = 0 ; i < SECURITY_LEVEL/4 + sizeof(int) ; ++i )
    {
        randomness[i] = seed_S[i];
    };
    // append indicator
    type = 0; // S
    for( i = 0 ; i < 4 ; ++i )
    {
        randomness[SECURITY_LEVEL/4+i] = (unsigned char)(type >> (8*i));
    }

    // expand seed to buffer of mucho pseudorandomness
    FIPS202_SHAKE256(randomness, SECURITY_LEVEL/4 + 4, buffer, GFP_NUMBYTES*DEGREE_OF_CIRCULANCY * BACUOV_PARAM_V*BACUOV_PARAM_O);

    // use pseudorandomness buffer to sample circulant ring elements
    k = 0;
    for( i = 0 ; i < BACUOV_PARAM_V ; ++i )
    {
        for( j = 0 ; j < BACUOV_PARAM_O ; ++j )
        {
            for( k = 0 ; k < DEGREE_OF_CIRCULANCY ; ++k )
			{
				data->data[i*BACUOV_PARAM_O + j].data[k] = buffer[(i*BACUOV_PARAM_O + j)*DEGREE_OF_CIRCULANCY + k] % GF_PRIME_MODULUS;
			}
        }
    }
}

/**
 * bacuov_generate_vinegar_coefficients
 * Generate the top-left VxV block of circulant ring elements of the
 * public quadratic form Pi.
 */
void bacuov_generate_vinegar_coefficients( gfpcircmatrix * Pi0V0V, unsigned char * seed, int index )
{
    unsigned char randomness[SECURITY_LEVEL/4 + 8];
    int i, j, k, K;
    int type;
    gfpcirc_element elm;
    unsigned char buffer[GFP_NUMBYTES*DEGREE_OF_CIRCULANCY * (BACUOV_PARAM_V * (BACUOV_PARAM_V + 1) / 2)];
    
    // copy over randomness
    for( i = 0 ; i < SECURITY_LEVEL/4 + sizeof(int) ; ++i )
    {
        randomness[i] = seed[i];
    };
    // append indicator
    // append indicators
    type = 1; // vinegars
    for( i = 0 ; i < 4 ; ++i )
    {
        randomness[SECURITY_LEVEL/4+i] = (unsigned char)(type >> (8*i));
        randomness[SECURITY_LEVEL/4 + 4 + i] = (unsigned char)(index >> (8*i));
    }

    // expand seed to buffer of mucho pseudorandomness
    FIPS202_SHAKE256(randomness, SECURITY_LEVEL/4 + 8, buffer, GFP_NUMBYTES*DEGREE_OF_CIRCULANCY * (BACUOV_PARAM_V * (BACUOV_PARAM_V + 1) / 2));

    // use pseudorandomness buffer to sample circulant ring elements
	K = 0;
    for( i = 0 ; i < BACUOV_PARAM_V ; ++i )
    {
		for( k = 0 ; k < DEGREE_OF_CIRCULANCY ; ++k )
		{
        	Pi0V0V->data[i*BACUOV_PARAM_V+i].data[k] = buffer[K + k] % GF_PRIME_MODULUS;
		}
		K = K + DEGREE_OF_CIRCULANCY;
        for( j = i+1 ; j < BACUOV_PARAM_V ; ++j )
        {
			for( k = 0 ; k < DEGREE_OF_CIRCULANCY ; ++k )
			{
            	elm.data[k] = buffer[K + k] % GF_PRIME_MODULUS;
			}
			K = K + DEGREE_OF_CIRCULANCY;
            // divide elm by 2? not necessary
            gfpcirc_copy(&Pi0V0V->data[i*BACUOV_PARAM_V+j], &elm);
            gfpcirc_copy(&Pi0V0V->data[j*BACUOV_PARAM_V+i], &elm);
        }
    }
}

/**
 * bacuov_generate_linear_coefficients
 * Generate the top-right vinegar-oil part of the secret quadratic
 * form Pi.
 */
void bacuov_generate_linear_coefficients( gfpcircmatrix * Pi0VVN, unsigned char * seed, int index )
{
    unsigned char randomness[SECURITY_LEVEL/4 + 8];
    int type;
    int i, j, k;
    gfpcirc_element elm;
    unsigned char buffer[sizeof(gfpcirc_element) * (BACUOV_PARAM_V * BACUOV_PARAM_O)];
    
    // copy over randomness
    for( i = 0 ; i < SECURITY_LEVEL/4 ; ++i )
    {
        randomness[i] = seed[i];
    };
    // append indicators
    type = 2; // linears
    for( i = 0 ; i < 4 ; ++i )
    {
        randomness[SECURITY_LEVEL/4+i] = (unsigned char)(type >> (8*i));
        randomness[SECURITY_LEVEL/4 + 4 + i] = (unsigned char)(index >> (8*i));
    }

    // expand seed to buffer of mucho pseudorandomness
    FIPS202_SHAKE256(randomness, SECURITY_LEVEL/4 + 8, buffer, GFP_NUMBYTES * DEGREE_OF_CIRCULANCY * (BACUOV_PARAM_V * BACUOV_PARAM_O));

    // use pseudorandomness buffer to sample circulant ring elements
    k = 0;
    for( i = 0 ; i < BACUOV_PARAM_V ; ++i )
    {
        for( j = 0 ; j < BACUOV_PARAM_O ; ++j )
        {
            for( k = 0 ; k < DEGREE_OF_CIRCULANCY ; ++k )
            {
            	Pi0VVN->data[i*BACUOV_PARAM_O+j].data[k] = buffer[(i*BACUOV_PARAM_O+j)*DEGREE_OF_CIRCULANCY + k] % GF_PRIME_MODULUS;
            }
        }
    }
}

void bacuov_secret_key_init( bacuov_secret_key * sk )
{
    int i;
    sk->S = gfpcircm_init(BACUOV_PARAM_V, BACUOV_PARAM_O);
    for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
    {
        sk->FFvv[i] = gfpcircm_init(BACUOV_PARAM_V, BACUOV_PARAM_V);
        sk->FFvo[i] = gfpcircm_init(BACUOV_PARAM_V, BACUOV_PARAM_O);
    }
}

void bacuov_secret_key_destroy( bacuov_secret_key sk )
{
    int i;
    gfpcircm_destroy(sk.S);
    for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
    {
        gfpcircm_destroy(sk.FFvv[i]);
        gfpcircm_destroy(sk.FFvo[i]);
    }
}

void bacuov_public_key_init( bacuov_public_key * pk )
{
    int i;
    for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
    {
        pk->PPoo[i] = gfpcircm_init(BACUOV_PARAM_O, BACUOV_PARAM_O);
    }
}

void bacuov_public_key_destroy( bacuov_public_key pk )
{
    int i;
    for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
    {
        gfpcircm_destroy(pk.PPoo[i]);
    }
}

/**
 * bacuov_keygen
 * Generate secrey and public keypair for BAC UOV cryptosystem.
 */
void bacuov_keygen( bacuov_secret_key * sk, bacuov_public_key * pk, unsigned char * randomness )
{
    int i, j, k;
    unsigned char buffer[SECURITY_LEVEL/4 + SECURITY_LEVEL/4];
    unsigned char seed_S[SECURITY_LEVEL/4];
    gfpcircmatrix Pi0V0V, Pi0VVN;
    gfpcircmatrix tempVV, tempVO, tempOV, matVO, matOO;

    // initialize matrices
    Pi0V0V = gfpcircm_init(BACUOV_PARAM_V, BACUOV_PARAM_V);
    Pi0VVN = gfpcircm_init(BACUOV_PARAM_V, BACUOV_PARAM_O);
    tempVV = gfpcircm_init(BACUOV_PARAM_V, BACUOV_PARAM_V);
    tempVO = gfpcircm_init(BACUOV_PARAM_V, BACUOV_PARAM_O);
    tempOV = gfpcircm_init(BACUOV_PARAM_O, BACUOV_PARAM_V);
    //matVO = gfpcircm_init(BACUOV_PARAM_V, BACUOV_PARAM_O);
    matOO = gfpcircm_init(BACUOV_PARAM_O, BACUOV_PARAM_O);

    // copy over randomness to secret key
    for( i = 0 ; i < SECURITY_LEVEL/4 ; ++i )
    {
        sk->seed[i] = randomness[i];
    }

    // expand randomness into seeds
    FIPS202_SHAKE256(randomness, SECURITY_LEVEL/4, buffer, 2*SECURITY_LEVEL/4);

    // grab seed for S and PCT
    for( i = 0 ; i < SECURITY_LEVEL/4 ; ++i )
    {
        seed_S[i] = buffer[i];
        pk->seed_PCT[i] = buffer[SECURITY_LEVEL/4 + i];
    }

    // expand seed for S into S proper
    bacuov_generate_S(&sk->S, seed_S);

    // for all m polynomials
    for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
    {
        // populate top-left and top-right blocks of Pi
        bacuov_generate_vinegar_coefficients(&Pi0V0V, pk->seed_PCT, i);
        bacuov_generate_linear_coefficients(&Pi0VVN, pk->seed_PCT, i);

        // copy FF[i][0:V,0:V]
        gfpcircm_copy(sk->FFvv[i], Pi0V0V);

        // compute FF[i][0:V,V:N]
        gfpcircm_multiply(&sk->FFvo[i], Pi0V0V, sk->S);
        gfpcircm_negate(sk->FFvo[i]);

        //gfpcircm_copy(matVO, Pi0VVN);
        //gfpcircm_add(sk->FFvo[i], sk->FFvo[i], matVO);
	gfpcircm_add(sk->FFvo[i], sk->FFvo[i], Pi0VVN);

	// TODO: get rid of more temp variables
	
        // compute PP[i][V:N,V:N]
        gfpcircm_transpose(&sk->S);
        gfpcircm_copy(tempVO, sk->FFvo[i]);
        gfpcircm_multiply(&matOO, sk->S, tempVO);
        gfpcircm_transpose(&matOO);

        gfpcircm_copy(tempVV, sk->FFvv[i]);
        gfpcircm_multiply(&tempOV, sk->S, tempVV);
        gfpcircm_transpose(&sk->S);
        gfpcircm_multiply(&pk->PPoo[i], tempOV, sk->S);

        gfpcircm_add(pk->PPoo[i], pk->PPoo[i], matOO);
        gfpcircm_transpose(&matOO);
        gfpcircm_add(pk->PPoo[i], pk->PPoo[i], matOO);
        gfpcircm_transpose(&matOO);
    }

    // destroy matrices
    gfpcircm_destroy(Pi0V0V);
    gfpcircm_destroy(Pi0VVN);
    gfpcircm_destroy(tempVV);
    gfpcircm_destroy(tempVO);
    gfpcircm_destroy(tempOV);
    gfpcircm_destroy(matOO);
    //gfpcircm_destroy(matVO);
}

/**
 * gfpcircm_press_anticirculant
 * Map a matrix over the circulant ring down to the base ring using
 * anticirculant expansion.
 */
void gfpcircm_press_anticirculant( gfpmatrix * dest, gfpcircmatrix source )
{
	int i, j, k, l, d;
	d = DEGREE_OF_CIRCULANCY;
	for( i = 0 ; i < source.height ; ++i )
	{
		for( j = 0 ; j < source.width ; ++j )
		{
			// copy first row
			for( k = 0 ; k < d ; ++k )
			{
				dest->data[i*d * dest->width + j*d + d-1-k] = source.data[i*source.width + j].data[k];
			}
			// infer remaining rows
			for( l = 1 ; l < d ; ++l )
			{
				for( k = 0 ; k < d-1 ; ++k )
				{
					dest->data[(i*d+l) * dest->width + j*d + k] = dest->data[(i*d+l-1) * dest->width + j*d + k + 1];
				}
				dest->data[(i*d+l) * dest->width + j*d + d-1] = dest->data[(i*d+l-1) * dest->width + j*d];
			}
		}
	}
}


/**
 * gfpcircm_press_circulant
 * Map a matrix over the circulant ring down to the base ring using
 * anticirculant expansion.
 */
void gfpcircm_press_circulant( gfpmatrix * dest, gfpcircmatrix source )
{
	int i, j, k, l, d;
	d = DEGREE_OF_CIRCULANCY;
	for( i = 0 ; i < source.height ; ++i )
	{
		for( j = 0 ; j < source.width ; ++j )
		{
			// copy first row
			for( k = 0 ; k < d ; ++k )
			{
				dest->data[i*d * dest->width + j*d + k] = source.data[i*source.width + j].data[k];
			}
			// infer remaining rows
			for( l = 1 ; l < d ; ++l )
			{
				dest->data[(i*d+l) * dest->width + j*d] = dest->data[(i*d+l-1) * dest->width + j*d + d-1];
				for( k = 0 ; k < d-1 ; ++k )
				{
					dest->data[(i*d+l) * dest->width + j*d + k+1] = dest->data[(i*d+l-1) * dest->width + j*d + k];
				}
			}
		}
	}
}

/**
 * gfpem_multiply_base
 * Multiply a matrix over extension of GF(p), with a matrix over
 * GF(p) and store the result in a designated variable.
 */
void gfpem_multiply_base( gfpematrix * dest, gfpematrix lhs, gfpmatrix rhs )
{
	int i, j, k, l;
    unsigned long int temp[EXTENSION_DEGREE];
	
	for( i = 0 ; i < lhs.height ; ++i )
	{
		for( j = 0 ; j < rhs.width ; ++j )
		{
            for( l = 0 ; l < EXTENSION_DEGREE ; ++l )
            {
                temp[l] = 0;
            }
			for( k = 0 ; k < lhs.width ; ++k )
			{
                for( l = 0 ; l < EXTENSION_DEGREE ; ++l )
                {
                    temp[l] = temp[l] + lhs.data[i*lhs.width + k].data[l] * rhs.data[k*rhs.width + j];
                }
			}
            for( l = 0 ; l < EXTENSION_DEGREE ; ++l )
            {
                dest->data[i*dest->width + j].data[l] = temp[l] % GF_PRIME_MODULUS;
            }
		}
	}
}

/**
 * gfpem_base_multiply
 * Multiply a matrix over GF(p), with a matrix over an extension of
 * GF(p) and store the result in a designated variable.
 */
void gfpem_base_multiply( gfpematrix * dest, gfpmatrix lhs, gfpematrix rhs )
{
	int i, j, k, l;
    unsigned long temp[EXTENSION_DEGREE];
	
	for( i = 0 ; i < lhs.height ; ++i )
	{
		for( j = 0 ; j < rhs.width ; ++j )
		{
            for( l = 0 ; l < EXTENSION_DEGREE ; ++l )
            {
                temp[l] = 0;
            }
			for( k = 0 ; k < lhs.width ; ++k )
			{
                for( l = 0 ; l < EXTENSION_DEGREE ; ++l )
                {
                    temp[l] = temp[l] + lhs.data[i*lhs.width + k] * rhs.data[k*rhs.width + j].data[l];
                }
			}
            for( l = 0 ; l < EXTENSION_DEGREE ; ++l )
            {
                dest->data[i*dest->width + j].data[l] = temp[l] % GF_PRIME_MODULUS;
            }
		}
	}
}

void bacuov_sign( bacuov_signature * sig, bacuov_secret_key sk, const unsigned char * msg, int msg_len )
{

    unsigned char randomness[sizeof(sk.seed) + 32 + 1];
    unsigned char array[BACUOV_PARAM_M * EXTENSION_DEGREE];
	unsigned char vinegar_array[BACUOV_PARAM_V*DEGREE_OF_CIRCULANCY * EXTENSION_DEGREE];
	unsigned int trial_index;
    int i, j, k;
	gfpcirc_element elm;
	gfpe_element element;
    gfpe_element inner_product, temp;
    gfpematrix target;
	gfpematrix coeff, inverse;
	gfpematrix xv;
	gfpematrix xo;
	gfpematrix coeff_line;
	gfpematrix FFvoi;
	gfpmatrix FFvoi_temp, FFvvi;
    gfpmatrix S_, S__;
    gfpematrix xv_FFvvi;
    gfpcircmatrix S;
    gfpematrix x, s;
    gfpe_element signature[(BACUOV_PARAM_O+BACUOV_PARAM_V)*DEGREE_OF_CIRCULANCY];
    gfpematrix update;

    // init matrices
    target = gfpem_init(BACUOV_PARAM_M, 1);
	coeff = gfpem_init(BACUOV_PARAM_M, BACUOV_PARAM_M);
	inverse = gfpem_init(BACUOV_PARAM_M, BACUOV_PARAM_M);
	xv = gfpem_init(1, BACUOV_PARAM_V*DEGREE_OF_CIRCULANCY);
	xo = gfpem_init(BACUOV_PARAM_M, 1);
	coeff_line = gfpem_init(1, BACUOV_PARAM_M);
    FFvoi = gfpem_init(BACUOV_PARAM_V*DEGREE_OF_CIRCULANCY, BACUOV_PARAM_O*DEGREE_OF_CIRCULANCY);
	FFvoi_temp = gfpm_init(BACUOV_PARAM_V*DEGREE_OF_CIRCULANCY, BACUOV_PARAM_O*DEGREE_OF_CIRCULANCY);
    FFvvi = gfpm_init(BACUOV_PARAM_V*DEGREE_OF_CIRCULANCY, BACUOV_PARAM_V*DEGREE_OF_CIRCULANCY);
    xv_FFvvi = gfpem_init(1, BACUOV_PARAM_V*DEGREE_OF_CIRCULANCY);
    S_ = gfpm_init((BACUOV_PARAM_O+BACUOV_PARAM_V)*DEGREE_OF_CIRCULANCY, (BACUOV_PARAM_O+BACUOV_PARAM_V)*DEGREE_OF_CIRCULANCY);
    S__ = gfpm_init((BACUOV_PARAM_V)*DEGREE_OF_CIRCULANCY, (BACUOV_PARAM_O)*DEGREE_OF_CIRCULANCY);
    S = gfpcircm_init(BACUOV_PARAM_O+BACUOV_PARAM_V, BACUOV_PARAM_O+BACUOV_PARAM_V);
    x = gfpem_init((BACUOV_PARAM_V+BACUOV_PARAM_O)*DEGREE_OF_CIRCULANCY, 1);
    s = gfpem_init((BACUOV_PARAM_V+BACUOV_PARAM_O)*DEGREE_OF_CIRCULANCY, 1);
    update = gfpem_init(BACUOV_PARAM_V*DEGREE_OF_CIRCULANCY, 1);

    // get randomness from
    //  - randomness used for key generation
    //  - hash of document
    // (if a true RNG is available, that is preferrable.)
    for( i = 0 ; i < SECURITY_LEVEL/4 ; ++i )
    {
        randomness[i] = sk.seed[i];
    }

    FIPS202_SHA3_256(msg, msg_len, randomness + SECURITY_LEVEL/4);

    // hash to vector
    FIPS202_SHAKE256(msg, msg_len, array, BACUOV_PARAM_M * EXTENSION_DEGREE);
    for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
    {
        for( j = 0 ; j < EXTENSION_DEGREE ; ++j )
        {
            target.data[i].data[j] = array[i*EXTENSION_DEGREE + j] % GF_PRIME_MODULUS;
        }
    }

	// try random assignment for vinegar variables until we have an invertible coefficient matrix
	for( trial_index = 0 ; trial_index < 256 ; ++trial_index )
	{
		// use SHAKE to extract randomness from which to sample vinegar variables
		randomness[sizeof(randomness)-1] = trial_index & 0xff;
		FIPS202_SHAKE256(randomness, sizeof(randomness), vinegar_array, sizeof(vinegar_array));
		for( i = 0 ; i < BACUOV_PARAM_V*DEGREE_OF_CIRCULANCY ; ++i )
		{
			for( j = 0 ; j < EXTENSION_DEGREE ; ++j )
			{
				xv.data[i].data[j] = vinegar_array[i*EXTENSION_DEGREE + j] % GF_PRIME_MODULUS;
			}
		}

		// compile coefficient matrix line by line
		for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
		{
			gfpcircm_press_anticirculant(&FFvoi_temp, sk.FFvo[i]);
			gfpem_multiply_base(&coeff_line, xv, FFvoi_temp);

			// copy line (and multiply by two, while you're at it)
			for( j = 0 ; j < BACUOV_PARAM_O*DEGREE_OF_CIRCULANCY ; ++j )
			{
				for( k = 0 ; k < EXTENSION_DEGREE ; ++k )
				{
					coeff.data[i*BACUOV_PARAM_M + j].data[k] = (2 * coeff_line.data[j].data[k]) % GF_PRIME_MODULUS;
				}
			}
		}
	
		if( gfpem_inverse(inverse, coeff) == 1 )
			break;
	}

	// update target vector
	for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
    {
        // press vinegar-vinegar part of FF[i] down to base field
        gfpcircm_press_anticirculant(&FFvvi, sk.FFvv[i]);

        // multiply vinegar vector with vinegar part of FF[i]
        // store the result in temporary variable
        gfpem_multiply_base(&xv_FFvvi, xv, FFvvi);

        // compute inner product between that temporary variable
        // and the vinegar vector
        gfpe_zero(&inner_product);
        for( j = 0 ; j < BACUOV_PARAM_V*DEGREE_OF_CIRCULANCY ; ++j )
        {
            gfpe_zero(&temp);
            gfpe_multiply(&temp, xv_FFvvi.data[j], xv.data[j]);
            gfpe_add(&inner_product, inner_product, temp);
        }

        // subtract that from the target
        gfpe_subtract(&target.data[i], target.data[i], inner_product);
    }

    // find oil vector
    gfpem_multiply(&xo, inverse, target);

    // get signature from S^-1 * x
    // target for optimization: compute only the change to xv, as xo remains the same
    gfpcircm_eye(S);
    for( i = 0 ; i < BACUOV_PARAM_V ; ++i )
    {
        for( j = BACUOV_PARAM_V ; j < BACUOV_PARAM_V + BACUOV_PARAM_O ; ++j )
        {
			gfpcirc_copy(&elm, &sk.S.data[i*sk.S.width + j-BACUOV_PARAM_V]);
			gfpcirc_flr(&elm);
            gfpcirc_subtract(&S.data[i * S.width + j], S.data[i * S.width + j], elm);
        }
    }

    gfpcircm_press_circulant(&S_, S);

    for( i = 0 ; i < BACUOV_PARAM_V * DEGREE_OF_CIRCULANCY ; ++i )
    {
        gfpe_copy(&x.data[i], xv.data[i]);
    }
    for( i = 0  ; i < BACUOV_PARAM_O * DEGREE_OF_CIRCULANCY ; ++i )
    {
        gfpe_copy(&x.data[BACUOV_PARAM_V*DEGREE_OF_CIRCULANCY+i], xo.data[i]);
    }

    gfpem_base_multiply(&s, S_, x);

    for( i = 0 ; i < (BACUOV_PARAM_O + BACUOV_PARAM_V)*DEGREE_OF_CIRCULANCY ; ++i )
    {
        gfpe_copy(&sig->vector[i], s.data[i]);
    }

    // destroy matrices
    gfpem_destroy(target);
	gfpem_destroy(coeff);
	gfpem_destroy(inverse);
	gfpem_destroy(xv);
	gfpem_destroy(xo);
	gfpem_destroy(coeff_line);
	gfpem_destroy(FFvoi);
	gfpm_destroy(FFvoi_temp);
	gfpm_destroy(FFvvi);
    gfpcircm_destroy(S);
	gfpm_destroy(S_);
	gfpm_destroy(S__);
	gfpem_destroy(xv_FFvvi);
    gfpem_destroy(x);
    gfpem_destroy(s);
    gfpem_destroy(update);
}

int bacuov_verify( bacuov_public_key pk, const unsigned char * msg, int msg_len, bacuov_signature * sig )
{
    gfpcircmatrix Pi0V0V, Pi0VVN;
    gfpmatrix Pi0V0V_pressed, Pi0VVN_pressed, PiVNVN_pressed, Pi0N0N_pressed;
    gfpematrix evaluation, target;
    gfpematrix signature;
    gfpematrix temp_row, temp_inner_product;
    int valid;
    int i, j, k;
    unsigned char array[BACUOV_PARAM_O*DEGREE_OF_CIRCULANCY * EXTENSION_DEGREE];

    // init matrices
    Pi0V0V = gfpcircm_init(BACUOV_PARAM_V, BACUOV_PARAM_V);
    Pi0VVN = gfpcircm_init(BACUOV_PARAM_V, BACUOV_PARAM_O);
    Pi0V0V_pressed = gfpm_init(BACUOV_PARAM_V * DEGREE_OF_CIRCULANCY, BACUOV_PARAM_V * DEGREE_OF_CIRCULANCY);
    Pi0VVN_pressed = gfpm_init(BACUOV_PARAM_V * DEGREE_OF_CIRCULANCY, BACUOV_PARAM_O * DEGREE_OF_CIRCULANCY);
    PiVNVN_pressed = gfpm_init(BACUOV_PARAM_O * DEGREE_OF_CIRCULANCY, BACUOV_PARAM_O * DEGREE_OF_CIRCULANCY);
    Pi0N0N_pressed = gfpm_init((BACUOV_PARAM_V+BACUOV_PARAM_O) * DEGREE_OF_CIRCULANCY, (BACUOV_PARAM_V+BACUOV_PARAM_O) * DEGREE_OF_CIRCULANCY);
    evaluation = gfpem_init(BACUOV_PARAM_O * DEGREE_OF_CIRCULANCY, 1);
    target = gfpem_init(BACUOV_PARAM_O * DEGREE_OF_CIRCULANCY, 1);
    signature = gfpem_init((BACUOV_PARAM_V+BACUOV_PARAM_O) * DEGREE_OF_CIRCULANCY, 1);
    temp_row = gfpem_init(1, (BACUOV_PARAM_V+BACUOV_PARAM_O)*DEGREE_OF_CIRCULANCY);
    temp_inner_product = gfpem_init(1,1);

    // copy signature
    for( i = 0 ; i < (BACUOV_PARAM_V + BACUOV_PARAM_O) * DEGREE_OF_CIRCULANCY ; ++i )
    {
        gfpe_copy(&signature.data[i], sig->vector[i]);
    }

    // get P matrices from PRG and public key, and evaluate while at it
    for( i = 0 ; i < BACUOV_PARAM_O * DEGREE_OF_CIRCULANCY ; ++i )
    {
        // populate top-left and top-right blocks of Pi
        bacuov_generate_vinegar_coefficients(&Pi0V0V, pk.seed_PCT, i);
        bacuov_generate_linear_coefficients(&Pi0VVN, pk.seed_PCT, i);   

        // press down matrix parts
        gfpcircm_press_anticirculant(&Pi0V0V_pressed, Pi0V0V);
        gfpcircm_press_anticirculant(&Pi0VVN_pressed, Pi0VVN);
        gfpcircm_press_anticirculant(&PiVNVN_pressed, pk.PPoo[i]);

        // copy matrix parts into the big thing
        // starting with top-left vinegar-vinegar block ...
        for( j = 0 ; j < BACUOV_PARAM_V*DEGREE_OF_CIRCULANCY ; ++j )
        {
            for( k = 0 ; k < BACUOV_PARAM_V*DEGREE_OF_CIRCULANCY ; ++k )
            {
                Pi0N0N_pressed.data[j*Pi0N0N_pressed.width + k] = Pi0V0V_pressed.data[j*Pi0V0V_pressed.width + k];
            }
        }
        // ... then the vinegar-oil and oil-vinegar blocks ...
        for( j = 0 ; j < BACUOV_PARAM_V*DEGREE_OF_CIRCULANCY ; ++j )
        {
            for( k = 0 ; k < BACUOV_PARAM_O*DEGREE_OF_CIRCULANCY ; ++k )
            {
                Pi0N0N_pressed.data[j*Pi0N0N_pressed.width + BACUOV_PARAM_V*DEGREE_OF_CIRCULANCY + k] = Pi0VVN_pressed.data[j*Pi0VVN_pressed.width + k];
                Pi0N0N_pressed.data[(BACUOV_PARAM_V * DEGREE_OF_CIRCULANCY + k)*Pi0N0N_pressed.width + j] = Pi0VVN_pressed.data[j*Pi0VVN_pressed.width + k];
            }
        }
        // ... and lastly, the oil-oil block
        for( j = 0 ; j < BACUOV_PARAM_O*DEGREE_OF_CIRCULANCY ; ++j )
        {
            for( k = 0 ; k < BACUOV_PARAM_O*DEGREE_OF_CIRCULANCY ; ++k )
            {
                Pi0N0N_pressed.data[(BACUOV_PARAM_V*DEGREE_OF_CIRCULANCY+j)*Pi0N0N_pressed.width + BACUOV_PARAM_V*DEGREE_OF_CIRCULANCY + k] = PiVNVN_pressed.data[j*PiVNVN_pressed.width + k];
            }
        }

        // compute vector-matrix-vector product // <-- target for optimization
        signature.height = 1;
        signature.width = (BACUOV_PARAM_V+BACUOV_PARAM_O)*DEGREE_OF_CIRCULANCY;
        gfpem_multiply_base(&temp_row, signature, Pi0N0N_pressed);
        signature.width = 1;
        signature.height = (BACUOV_PARAM_V+BACUOV_PARAM_O)*DEGREE_OF_CIRCULANCY;
        gfpem_multiply(&temp_inner_product, temp_row, signature); 
        gfpe_copy(&evaluation.data[i], temp_inner_product.data[0]);
    }

    // get target hash-of-document
    FIPS202_SHAKE256(msg, msg_len, array, sizeof(array));
    for( i = 0 ; i < BACUOV_PARAM_O*DEGREE_OF_CIRCULANCY ; ++i )
    {
        for( j = 0 ; j < EXTENSION_DEGREE ; ++j )
        {
            target.data[i].data[j] = array[i*EXTENSION_DEGREE + j] % GF_PRIME_MODULUS;
        }
    }

    // compare evaluation to target
    valid = gfpem_equals(evaluation, target);

    // destroy matrices
    gfpcircm_destroy(Pi0V0V);
    gfpcircm_destroy(Pi0VVN);
    gfpm_destroy(Pi0V0V_pressed);
    gfpm_destroy(Pi0VVN_pressed);
    gfpm_destroy(PiVNVN_pressed);
    gfpm_destroy(Pi0N0N_pressed);
    gfpem_destroy(evaluation);
    gfpem_destroy(target);
    gfpem_destroy(signature);
    gfpem_destroy(temp_row);
    gfpem_destroy(temp_inner_product);

    return valid;
}

/**
 * bacuov_pack
 * Turn a list of field elements, represented as bytes, into a much
 * denser string of bytes.
 */
void bacuov_pack( unsigned char * dest_buffer, unsigned char * source_elements, int num_elements )
{
    int i, j;
    int position;
    unsigned char bit;

    position = 0;

    // set to zero
    for( i = 0 ; i < bacuov_pack_length(num_elements) ; ++i )
    {
        dest_buffer[i] = 0;
    }

    // for each element ...
    for( i = 0 ; i < num_elements ; ++i )
    {
        // for each bit in the element
        for( j = 0 ; j < GFP_NUMBITS ; ++j )
        {
            bit = (source_elements[i] >> j) & 1;
            // set the matching bit in the destination buffer
            bit = bit << (position % 8);
            dest_buffer[position/8] = dest_buffer[position/8] | bit;
            position = position + 1;
        }
    }
}

/**
 * bacuov_unpack
 * Turn a dense stream of bytes representing field elements, into a
 * list of such field elements where each byte represents one indi-
 * vidually.
 */
void bacuov_unpack( unsigned char * dest_elements, const unsigned char * source_buffer, int num_elements )
{
    int i, j;
    int position;
    unsigned char bit;

    position = 0;

    // set destination to zero
    for( i = 0 ; i < num_elements ; ++i )
    {
        dest_elements[i] = 0;
    }

    // for each element ...
    for( i = 0 ; i < num_elements ; ++i )
    {
        // for each bit in the element
        for( j = 0 ; j < GFP_NUMBITS ; ++j )
        {
            bit = (source_buffer[position/8] >> (position % 8)) & 1;
            bit = bit << (j % GFP_NUMBITS);
            dest_elements[i] = dest_elements[i] ^ bit;
            position = position + 1;
        }
    }
}

