#include "api.h"
#include <stdio.h>

void randombytes( unsigned char * dest, unsigned int length );

/**
 * crypto_sign_keypair
 * Generate keys; output as serialized byte strings.
 */
int crypto_sign_keypair(unsigned char *pk_buffer, unsigned char *sk_buffer)
{
    unsigned char seed[SECURITY_LEVEL/4];
    unsigned char unpacked_pk_buffer[BACUOV_PARAM_M*(BACUOV_PARAM_O*(BACUOV_PARAM_O+1)/2)*DEGREE_OF_CIRCULANCY];
    bacuov_secret_key sk;
    bacuov_public_key pk;
    unsigned int i, j, k, l, m, n;
    
    // get randomness 
    randombytes(seed, SECURITY_LEVEL/4);

    // generate keys as mathematical objects
    bacuov_secret_key_init(&sk);
    bacuov_public_key_init(&pk);
    bacuov_keygen(&sk, &pk, seed);

    // serialize mathematical objects: first secret key ...
    for( i = 0 ; i < SECURITY_LEVEL/4 ; ++i )
    {
        sk_buffer[i] = sk.seed[i];
    }
    //temp_buffer_length =   BACUOV_PARAM_M*(BACUOV_PARAM_V*(BACUOV_PARAM_V+1)/2)*DEGREE_OF_CIRCULANCY // FFvv
    //                     + BACUOV_PARAM_M*(BACUOV_PARAM_V*BACKUOV_PARAM_O)*DEGREE_OF_CIRCULANCY      // FFvo
    //                     + BACUOV_PARAM_V*BACUOV_PARAM_O * DEGREE_OF_CIRCULANCY;                     // S
    //temp_buffer = malloc(temp_buffer_length);
    m = 0;
    for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
    {
        for( j = 0 ; j < BACUOV_PARAM_V ; ++j )
        {
            for( k = j ; k < BACUOV_PARAM_V ; ++k )
            {
                for( l = 0 ; l < DEGREE_OF_CIRCULANCY ; ++l )
                {
                    sk_buffer[SECURITY_LEVEL/4 + m++] = sk.FFvv[i].data[j*BACUOV_PARAM_V + k].data[l];
                }
            }
        }
    }
    n = 0;
    for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
    {
        for( j = 0 ; j < BACUOV_PARAM_V ; ++j )
        {
            for( k = 0 ; k < BACUOV_PARAM_O ; ++k )
            {
                for( l = 0 ; l < DEGREE_OF_CIRCULANCY ; ++l )
                {
                    sk_buffer[SECURITY_LEVEL/4 + m + n++] = sk.FFvo[i].data[j*BACUOV_PARAM_O+k].data[l];
                }
            }
        }
    }
    for( i = 0 ; i < BACUOV_PARAM_V ; ++i )
    {
        for( j = 0 ; j < BACUOV_PARAM_O ; ++j )
        {
            for( k = 0 ; k < DEGREE_OF_CIRCULANCY ; ++k )
            {
                sk_buffer[SECURITY_LEVEL/4 + m + n + i*BACUOV_PARAM_O*DEGREE_OF_CIRCULANCY + j*DEGREE_OF_CIRCULANCY + k] = sk.S.data[i*BACUOV_PARAM_O+j].data[k];
            }
        }
    }

    // and then the public key, with packing
    m = 0;
    for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
    {
        for( j = 0 ; j < BACUOV_PARAM_O ; ++j )
        {
            for( k = j ; k < BACUOV_PARAM_O ; ++k )
            {
                for( l = 0 ; l < DEGREE_OF_CIRCULANCY ; ++l )
                {
                    unpacked_pk_buffer[m++] = pk.PPoo[i].data[j*BACUOV_PARAM_O + k].data[l];
                }
            }
        }
    }
    bacuov_pack(pk_buffer, unpacked_pk_buffer, sizeof(unpacked_pk_buffer));
    for( i = 0 ; i < SECURITY_LEVEL/4 ; ++i )
    {
        pk_buffer[bacuov_pack_length(sizeof(unpacked_pk_buffer))+i] = pk.seed_PCT[i];
    }

    // free memory and return
    bacuov_public_key_destroy(pk);
    bacuov_secret_key_destroy(sk);

    return 0;
}

int crypto_sign(unsigned char *sm, unsigned long long *smlen, const unsigned char *msg, unsigned long long msglen, const unsigned char *sk_buffer)
{
    unsigned int i, j, k, l, m, n;
    bacuov_secret_key sk;
    bacuov_signature signature;

    bacuov_secret_key_init(&sk);

    // deserialize secret key
    for( i = 0 ; i < SECURITY_LEVEL/4 ; ++i )
    {
        sk.seed[i] = sk_buffer[i];
    }
    //temp_buffer_length =   BACUOV_PARAM_M*(BACUOV_PARAM_V*(BACUOV_PARAM_V+1)/2)*DEGREE_OF_CIRCULANCY // FFvv
    //                     + BACUOV_PARAM_M*(BACUOV_PARAM_V*BACKUOV_PARAM_O)*DEGREE_OF_CIRCULANCY      // FFvo
    //                     + BACUOV_PARAM_V*BACUOV_PARAM_O * DEGREE_OF_CIRCULANCY;                     // S
    //temp_buffer = malloc(temp_buffer_length);
    m = 0;
    for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
    {
        for( j = 0 ; j < BACUOV_PARAM_V ; ++j )
        {
            for( k = j ; k < BACUOV_PARAM_V ; ++k )
            {
                for( l = 0 ; l < DEGREE_OF_CIRCULANCY ; ++l )
                {
                    sk.FFvv[i].data[j*BACUOV_PARAM_V + k].data[l] = sk_buffer[SECURITY_LEVEL/4 + m];
                    sk.FFvv[i].data[k*BACUOV_PARAM_V + j].data[l] = sk_buffer[SECURITY_LEVEL/4 + m];
                    m = m + 1;
                }
            }
        }
    }
    n = 0;
    for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
    {
        for( j = 0 ; j < BACUOV_PARAM_V ; ++j )
        {
            for( k = 0 ; k < BACUOV_PARAM_O ; ++k )
            {
                for( l = 0 ; l < DEGREE_OF_CIRCULANCY ; ++l )
                {
                    sk.FFvo[i].data[j*BACUOV_PARAM_O+k].data[l] = sk_buffer[SECURITY_LEVEL/4 + m + n++];
                }
            }
        }
    }
    for( i = 0 ; i < BACUOV_PARAM_V ; ++i )
    {
        for( j = 0 ; j < BACUOV_PARAM_O ; ++j )
        {
            for( k = 0 ; k < DEGREE_OF_CIRCULANCY ; ++k )
            {
                sk.S.data[i*BACUOV_PARAM_O+j].data[k] = sk_buffer[SECURITY_LEVEL/4 + m + n + i*BACUOV_PARAM_O*DEGREE_OF_CIRCULANCY + j*DEGREE_OF_CIRCULANCY + k];
            }
        }
    }

    // sign
    bacuov_sign(&signature, sk, msg, msglen);

    // serialize signature
    bacuov_pack(sm + msglen, (unsigned char *)signature.vector, (BACUOV_PARAM_O+BACUOV_PARAM_V) * DEGREE_OF_CIRCULANCY * EXTENSION_DEGREE);
    for( i = 0 ; i < msglen ; ++i )
    {
        sm[i] = msg[i];
    }
    *smlen = msglen + CRYPTO_BYTES;

    // free memory and return
    bacuov_secret_key_destroy(sk);

    return 0;
}

int crypto_sign_open(unsigned char *msg, unsigned long long *msglen, const unsigned char *sm, unsigned long long smlen, const unsigned char *pk_buffer)
{
    bacuov_public_key pk;
    bacuov_signature signature;
    int signature_is_valid;
    unsigned char unpacked_pk_buffer[BACUOV_PARAM_M*(BACUOV_PARAM_O*(BACUOV_PARAM_O+1)/2)*DEGREE_OF_CIRCULANCY];
    unsigned int i, j, k, l, m, n;

    // init public key
    bacuov_public_key_init(&pk);

    // deserialize the public key, with packing
    bacuov_unpack(unpacked_pk_buffer, pk_buffer, sizeof(unpacked_pk_buffer));
    m = 0;
    for( i = 0 ; i < BACUOV_PARAM_M ; ++i )
    {
        for( j = 0 ; j < BACUOV_PARAM_O ; ++j )
        {
            for( k = j ; k < BACUOV_PARAM_O ; ++k )
            {
                for( l = 0 ; l < DEGREE_OF_CIRCULANCY ; ++l )
                {
                    pk.PPoo[i].data[j*BACUOV_PARAM_O + k].data[l] = unpacked_pk_buffer[m];
                    pk.PPoo[i].data[k*BACUOV_PARAM_O + j].data[l] = unpacked_pk_buffer[m];
                    m = m + 1;
                }
            }
        }
    }
    for( i = 0 ; i < SECURITY_LEVEL/4 ; ++i )
    {
        pk.seed_PCT[i] = pk_buffer[bacuov_pack_length(sizeof(unpacked_pk_buffer))+i];
    }

    // deserialize signature, with packing
    bacuov_unpack((unsigned char *)signature.vector, sm + smlen - CRYPTO_BYTES, (BACUOV_PARAM_V+BACUOV_PARAM_O)*DEGREE_OF_CIRCULANCY*EXTENSION_DEGREE);

    // test validity
    signature_is_valid = bacuov_verify(pk, sm, smlen - CRYPTO_BYTES, &signature);

    // copy message, if necessary
    if( signature_is_valid == 1 )
    {
        *msglen = smlen - bacuov_pack_length(sizeof(signature));
        for( i = 0 ; i < *msglen ; ++i )
        {
            msg[i] = sm[i];
        }
    }

    // free memory
    bacuov_public_key_destroy(pk);

    return !signature_is_valid;
}

