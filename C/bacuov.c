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
    //printf("buffer for S: ");
    //for( i = 0 ; i < GFP_NUMBYTES*DEGREE_OF_CIRCULANCY * BACUOV_PARAM_V * BACUOV_PARAM_O ; ++i )
    //    printf("%02x", buffer[i]);
    //printf("\n");

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

    //printf("S:\n");
    //printf("%i x %i matrix (should be %i x %i)\n", data->height, data->width, BACUOV_PARAM_V, BACUOV_PARAM_O);
    //gfpcircm_print(*data);
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

    //printf("buffer: ");
    //for( i = 0 ; i < GFP_NUMBYTES*DEGREE_OF_CIRCULANCY * BACUOV_PARAM_V * (BACUOV_PARAM_V+1) / 2 ; ++i )
    //    printf("%02x", buffer[i]);
    //printf(" <-- got it!\n");

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

    //printf("buffer: ");
    //for( i = 0 ; i < GFP_NUMBYTES * DEGREE_OF_CIRCULANCY * BACUOV_PARAM_V * BACUOV_PARAM_O ; ++i )
    //    printf("%02x", buffer[i]);
    //printf("\n");

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

    //printf("inside keygen with randomness ");
    //for( i = 0 ; i < SECURITY_LEVEL/4 ; ++i )
    //{
    //    printf("%02x", randomness[i]);
    //}
    //printf("\n");

    // initialize matrices
    Pi0V0V = gfpcircm_init(BACUOV_PARAM_V, BACUOV_PARAM_V);
    Pi0VVN = gfpcircm_init(BACUOV_PARAM_V, BACUOV_PARAM_O);
    tempVV = gfpcircm_init(BACUOV_PARAM_V, BACUOV_PARAM_V);
    tempVO = gfpcircm_init(BACUOV_PARAM_V, BACUOV_PARAM_O);
    tempOV = gfpcircm_init(BACUOV_PARAM_O, BACUOV_PARAM_V);
    matVO = gfpcircm_init(BACUOV_PARAM_V, BACUOV_PARAM_O);
    matOO = gfpcircm_init(BACUOV_PARAM_O, BACUOV_PARAM_O);

    // copy over randomness to secret key
    for( i = 0 ; i < SECURITY_LEVEL/4 ; ++i )
    {
        sk->seed[i] = randomness[i];
    }

    // expand randomness into seeds
	printf("SHAKE input: ");
	for( i = 0 ; i < SECURITY_LEVEL/4 ; ++i )
		printf("%02x", randomness[i]);
	printf("\n");
    FIPS202_SHAKE256(randomness, SECURITY_LEVEL/4, buffer, 2*SECURITY_LEVEL/4);
	printf("SHAKE output: ");
	for( i = 0 ; i < 2*SECURITY_LEVEL/4 ; ++i )
		printf("%02x", buffer[i]);
	printf("\n");


    // grab seed for S and PCT
    for( i = 0 ; i < SECURITY_LEVEL/4 ; ++i )
    {
        seed_S[i] = buffer[i];
        pk->seed_PCT[i] = buffer[SECURITY_LEVEL/4 + i];
    }
    //printf("PCT seed: ");
    //for( i = 0 ; i < SECURITY_LEVEL / 4 ; ++i )
    //    printf("%02x", pk->seed_PCT[i]);
    //printf("\n"); getchar();

    // expand seed for S into S proper
    bacuov_generate_S(&sk->S, seed_S);
	//printf("S:\n");
	//gfpcircm_print(sk->S);
	//getchar();


    // for all m polynomials
    for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
    {
        // populate top-left and top-right blocks of Pi
        bacuov_generate_vinegar_coefficients(&Pi0V0V, pk->seed_PCT, i);
        bacuov_generate_linear_coefficients(&Pi0VVN, pk->seed_PCT, i);

		//printf("Pi0VVN:\n");
		//gfpcircm_print(Pi0VVN);
		//getchar();

        // copy FF[i][0:V,0:V]
        gfpcircm_copy(sk->FFvv[i], Pi0V0V);
        gfpcircm_shiftflip(sk->FFvv[i],1);

        // compute FF[i][0:V,V:N]
        gfpcircm_multiply(&sk->FFvo[i], Pi0V0V, sk->S);
        gfpcircm_shiftflip(sk->FFvo[i],2);
        gfpcircm_negate(sk->FFvo[i]);

        gfpcircm_copy(matVO, Pi0VVN);
        gfpcircm_shiftflip(matVO, 1);
        gfpcircm_add(sk->FFvo[i], sk->FFvo[i], matVO);

		//printf("FFov:\n");
		//gfpcircm_print(sk->FFvo[i]);
		//getchar();

        // compute PP[i][V:N,V:N]
        gfpcircm_transpose(&sk->S);
        gfpcircm_copy(tempVO, sk->FFvo[i]);
        gfpcircm_flip(tempVO);
        gfpcircm_multiply(&matOO, sk->S, tempVO);
        gfpcircm_transpose(&matOO);

        gfpcircm_copy(tempVV, sk->FFvv[i]);
        gfpcircm_flipshift(tempVV, 1);
        gfpcircm_multiply(&tempOV, sk->S, tempVV);
        gfpcircm_transpose(&sk->S);
        gfpcircm_multiply(&pk->PPoo[i], tempOV, sk->S);

        gfpcircm_add(pk->PPoo[i], pk->PPoo[i], matOO);
        gfpcircm_transpose(&matOO);
        gfpcircm_add(pk->PPoo[i], pk->PPoo[i], matOO);
        gfpcircm_transpose(&matOO);
    }

    /*
    printf("FF[0][0:V, 0:V]:\n");
    gfpcircm_print(sk->FFvv[0]);
    printf("FF[0][0:V, V:N]:\n");
    gfpcircm_print(sk->FFvo[0]);
    printf("PP[0][V:N, V:N]:\n");
    gfpcircm_print(pk->PPoo[0]);
    */

    // destroy matrices
    
    gfpcircm_destroy(Pi0V0V);
    gfpcircm_destroy(Pi0VVN);
    gfpcircm_destroy(tempVV);
    gfpcircm_destroy(tempVO);
    gfpcircm_destroy(tempOV);
    gfpcircm_destroy(matOO);
    gfpcircm_destroy(matVO);

}

void gfpcircm_press( gfpmatrix * dest, gfpcircmatrix source )
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
 * gfpem_multiply_base
 * Multiply a matrix over extension of GF(p), with a matrix over
 * GF(p) and store the result in a designated variable.
 */
void gfpem_multiply_base( gfpematrix * dest, gfpematrix lhs, gfpmatrix rhs )
{
	gfpe_element r, t;
	int i, j, k;
	
	for( k = 0 ; k < EXTENSION_DEGREE ; ++k )
	{
		r.data[k] = 0;
	}

	for( i = 0 ; i < lhs.height ; ++i )
	{
		for( j = 0 ; j < rhs.width ; ++j )
		{
			gfpe_zero(&dest->data[i*dest->width + j]);
			for( k = 0 ; k < lhs.width ; ++k )
			{
				r.data[0] = rhs.data[k*rhs.width + j];
				gfpe_multiply(&t, lhs.data[i*lhs.width + k], r); // <-- target for optimization
				gfpe_add(&dest->data[i*dest->width + j], dest->data[i*dest->width + j], t);
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
	gfpe_element l, t;
	int i, j, k;
	
	for( k = 0 ; k < EXTENSION_DEGREE ; ++k )
	{
		l.data[k] = 0;
	}

	for( i = 0 ; i < lhs.height ; ++i )
	{
		for( j = 0 ; j < rhs.width ; ++j )
		{
			gfpe_zero(&dest->data[i*dest->width + j]);
			for( k = 0 ; k < lhs.width ; ++k )
			{
				l.data[0] = lhs.data[i*lhs.width + k];
				gfpe_multiply(&t, l, rhs.data[k*rhs.width + j]); // <-- target for optimization
				gfpe_add(&dest->data[i*dest->width + j], dest->data[i*dest->width + j], t);
			}
		}
	}
}

void bacuov_sign( bacuov_signature * sig, bacuov_secret_key sk, unsigned char * msg, int msg_len )
{

    unsigned char array[BACUOV_PARAM_M * EXTENSION_DEGREE + SECURITY_LEVEL/4 + 1];
	unsigned char vinegar_array[BACUOV_PARAM_V*DEGREE_OF_CIRCULANCY * EXTENSION_DEGREE];
	unsigned int trial_index;
    int i, j, k;
	gfpe_element element;
    gfpe_element inner_product, temp;
    gfpematrix target;
	gfpematrix coeff, inverse;
	gfpematrix xv;
	gfpematrix xo;
	gfpematrix coeff_line;
	gfpematrix FFvoi;
	gfpmatrix FFvoi_temp, FFvvi;
    gfpmatrix S_;
    gfpematrix xv_FFvvi;
    gfpcircmatrix S;
    gfpematrix x, s;

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
    S = gfpcircm_init(BACUOV_PARAM_O+BACUOV_PARAM_V, BACUOV_PARAM_O+BACUOV_PARAM_V);
    x = gfpem_init((BACUOV_PARAM_V+BACUOV_PARAM_O)*DEGREE_OF_CIRCULANCY, 1);
    s = gfpem_init((BACUOV_PARAM_V+BACUOV_PARAM_O)*DEGREE_OF_CIRCULANCY, 1);

    // hash to vector
    FIPS202_SHAKE256(msg, msg_len, array, BACUOV_PARAM_M * EXTENSION_DEGREE);
    for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
    {
        for( j = 0 ; j < EXTENSION_DEGREE ; ++j )
        {
            target.data[i].data[j] = array[i*EXTENSION_DEGREE + j] % GF_PRIME_MODULUS;
        }
    }

	// prepare array: add in randomness from secret key
	for( i = 0 ; i < SECURITY_LEVEL / 4 ; ++i )
	{
		array[BACUOV_PARAM_M * EXTENSION_DEGREE + i] = sk.seed[i];
	}

	// try random assignment for vinegar variables until we have an invertible coefficient matrix
	for( trial_index = 0 ; trial_index < 256 ; ++trial_index )
	{
		// use SHAKE to extract randomness from which to sample vinegar variables
		array[BACUOV_PARAM_M * EXTENSION_DEGREE + SECURITY_LEVEL/4] = trial_index & 0xff;
		FIPS202_SHAKE256(array, sizeof(array), vinegar_array, sizeof(vinegar_array));
		for( i = 0 ; i < BACUOV_PARAM_V*DEGREE_OF_CIRCULANCY ; ++i )
		{
			for( j = 0 ; j < EXTENSION_DEGREE ; ++j )
			{
				xv.data[i].data[j] = vinegar_array[i*EXTENSION_DEGREE + j] % GF_PRIME_MODULUS;
			}
		}
		//printf("vinegar variables:\n");
		//gfpem_print(xv);

		// compile coefficient matrix line by line
		for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
		{
			gfpcircm_press(&FFvoi_temp, sk.FFvo[i]);
			gfpem_multiply_base(&coeff_line, xv, FFvoi_temp);

			//if( i >= 1 )
			//{
			//printf("FFvoi:\n");
			//gfpm_print(FFvoi_temp);
			//printf("xv:\n");
			//gfpem_print(xv);
			////gfpcircm_print(sk.FFvo[i]);
			//getchar();
			//}
	
			//if( i >= 1 )
			//{
			//printf("line of coefficient matrix:\n");
			//gfpem_print(coeff_line);
			//getchar();
			//}

			// copy line (and multiply by two, while you're at it)
			for( j = 0 ; j < BACUOV_PARAM_O*DEGREE_OF_CIRCULANCY ; ++j )
			{
				for( k = 0 ; k < EXTENSION_DEGREE ; ++k )
				{
					coeff.data[i*BACUOV_PARAM_M + j].data[k] = (2 * coeff_line.data[j].data[k]) % GF_PRIME_MODULUS;
				}
			}
		}

		//printf("coefficient matrix:\n");
		//gfpem_print(coeff);
		
		if( gfpem_inverse(inverse, coeff) == 1 )
			break;
	}

	// update target vector
	for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
    {
        // press vinegar-vinegar part of FF[i] down to base field
        gfpcircm_press(&FFvvi, sk.FFvv[i]);

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
    gfpcircm_eye(S);
    for( i = 0 ; i < BACUOV_PARAM_V ; ++i )
    {
        for( j = BACUOV_PARAM_V ; j < BACUOV_PARAM_V + BACUOV_PARAM_O ; ++j )
        {
            gfpcirc_subtract(&S.data[i,j], S.data[i,j], sk.S.data[i, j-BACUOV_PARAM_V]);
        }
    }
    gfpcircm_press(&S_, S);
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
	gfpem_destroy(xv_FFvvi);
    gfpem_destroy(x);
    gfpem_destroy(s);
}

