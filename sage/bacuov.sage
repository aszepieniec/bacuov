import binascii
from CompactFIPS202 import SHAKE256, SHA3_256

def hash_digest( m ):
    return binascii.hexlify(SHA3_256(bytearray('1'+m)))

def hash_to_element( E, m ):
    hexmap = dict({'0': 0, '1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8, '9': 9, 'a': 10, 'b': 11, 'c': 12, 'd': 13, 'e': 14, 'f': 15})
    if E.is_prime_field():
        digest = hash_digest(m)
        integer = sum(hexmap[digest[i]] * 16^i for i in range(0, len(digest)))
        return integer % E.order()
    else:
        d = E.polynomial().degree()
        z = E.gen()
        F = E.polynomial().base_ring()
        poly = sum(hash_to_element(F, m+"||"+str(i)) * z^i for i in range(0, d))
        return poly

def hash_to_vector( E, d, m ):
    VS = MatrixSpace(E, d, 1)
    vec = copy(VS.zero())
    for i in range(0,d):
        vec[i,0] = hash_to_element(E, m+"||"+str(i))
    return vec

def anticirculant_matrix( vector ):
    MS = MatrixSpace(vector.parent().base_ring(), vector.nrows(), vector.nrows())
    mat = copy(MS.zero())
    mat[:,0] = vector
    for i in range(1, vector.nrows()):
        mat[:-1,i] = mat[1:,i-1]
        mat[-1,i] = mat[0,i-1]
    return mat

class bacuov_secret_key:
    def __init__( self, F, V, O, l ):
        Fx.<x> = PolynomialRing(F, "x")
        divisor = x^l - 1
        QR = QuotientRing(Fx, divisor)
        self.FF = [matrix([[QR.zero() for j in range(0, V+O)] for i in range(0, V+O)]) for k in range(0, O*l)]
        self.S = matrix([[QR.zero() for j in range(0, O)] for i in range(0, V)])

class bacuov_public_key:
    def __init__( self, F, V, O, l ):
        Fx.<x> = PolynomialRing(F, "x")
        divisor = x^l - 1
        QR = QuotientRing(Fx, divisor)
        self.PP = [matrix([[QR.zero() for j in range(0, O)] for i in range(0, O)]) for k in range(0, O*l)]
        self.seed_PCT = bytearray([0])

def bacuov_sample_circulant_ring_element( QR, buff ):
    Fx = QR.modulus().parent()
    x = Fx.gen()
    F = Fx.base_ring()
    p = F.order()
    d = QR.modulus().degree()

    integer = 0
    for i in range(0, len(buff)):
        integer = integer*256 + buff[i]

    expansion = []
    while integer > 0:
        expansion.append(integer % p)
        integer = floor(integer / p)

    poly = sum(x^i * F(expansion[i]) for i in range(0, min(len(expansion), d)))
    return QR(poly)

def bacuov_print_circulant_ring_element( elm ):
    d = elm.parent().modulus().degree()
    F = elm.parent().modulus().parent().base_ring()
    coeffs = elm.lift().coefficients(sparse=False)
    while len(coeffs) != d:
        coeffs.append(F.zero())
    coeffs = [str(coeffs[i]) for i in range(0,len(coeffs))]
    print "(" + ",".join(coeffs) + ")",

def bacuov_print_circulant_ring_matrix( mat ):
    print "[",
    for i in range(0, mat.nrows()):
        if i != 0:
            print " ",
        print "[",
        for j in range(0, mat.ncols()):
            bacuov_print_circulant_ring_element(mat[i,j])
        print "]",
        if i != mat.nrows()-1:
            print ""
    print "]"

def matrixify( qr_element ):
    QR = qr_element.parent()
    modulus = QR.modulus()
    d = modulus.degree()
    Fx = modulus.parent()
    F = Fx.base_ring()
    MS = MatrixSpace(F, d, d)

    coeffs = qr_element.lift().coefficients(sparse=False)

    mat = copy(MS.zero())
    for j in range(0, min(len(coeffs), d)):
        mat[0,j] = coeffs[j]
    for i in range(1, d):
        for j in range(0, d-1):
            mat[i,j] = mat[i-1,j+1];
        mat[i,d-1] = mat[i-1,0]

    return mat

def flip( qr_element ):
    QR = qr_element.parent()
    modulus = QR.modulus()
    d = modulus.degree()
    Fx = modulus.parent()
    F = Fx.base_ring()

    coeffs = qr_element.lift().coefficients(sparse=False)
    while len(coeffs) != d:
        coeffs.append(F.zero())

    coeffs = [coeffs[d-1-i] for i in range(0,d)]
    return sum(QR(coeffs[i]) * QR.gen()^i for i in range(0, len(coeffs)))

def bacuov_generate_S( QR, V, O, seed ):
    l = QR.modulus().degree()
    # append indicator
    for i in range(0, 4):
        seed += bytearray([int(0 >> (8*i))])

    # expand
    buff = SHAKE256(seed, l * V * O)
    print "buffer:", binascii.hexlify(buff)

    S = matrix([[QR.zero() for j in range(0,O)] for j in range(0,V)])
    k = 0;
    for i in range(0, V):
        for j in range(0, O):
            S[i,j] = bacuov_sample_circulant_ring_element(QR, buff[(l*k):(l*(k+1))])
            k = k + 1

    return S

def bacuov_generate_vinegar_coefficients( QR, V, rand, index ):

    l = QR.modulus().degree()
    # append indicator
    randomness = copy(rand)
    for i in range(0, 4):
        randomness += bytearray([int(1 >> (8*i))])
    for i in range(0, 4):
        randomness += bytearray([int((index >> (8*i)) % 256)])

    # expand
    buff = SHAKE256(randomness, l * V * (V+1) / 2)
    #print "buffer:", binascii.hexlify(buff), "<-- got it!"

    Pi = matrix([[QR.zero() for i in range(0,V)] for j in range(0,V)])
    k = 0;
    for i in range(0, V):
        Pi[i,i] = bacuov_sample_circulant_ring_element(QR, buff[(l*k):(l*(k+1))])
        k = k + 1
        for j in range(i+1, V):
            Pi[i,j] = bacuov_sample_circulant_ring_element(QR, buff[(l*k):(l*(k+1))])
            Pi[j,i] = Pi[i,j]
            k = k + 1

    return Pi
   
def bacuov_generate_linear_coefficients( QR, V, O, rand, index ): 

    l = QR.modulus().degree()
    # append indicator
    randomness = copy(rand)
    for i in range(0, 4):
        randomness += bytearray([int(2 >> (8*i))])
    for i in range(0, 4):
        randomness += bytearray([int((index >> (8*i)) % 256)])

    # expand
    buff = SHAKE256(randomness, l * V * O)
    print "buffer:", binascii.hexlify(buff)

    Pi = matrix([[QR.zero() for j in range(0,O)] for i in range(0,V)])
    k = 0;
    for i in range(0, V):
        for j in range(0, O):
            Pi[i,j] = bacuov_sample_circulant_ring_element(QR, buff[(l*k):(l*(k+1))])
            k = k + 1

    return Pi

def bacuov_keygen( SECURITY_LEVEL, F, V, O, l, randomness ):
    o = O * l
    v = V * l
    N = O + V
    n = o + v
    m = o

    Fx = PolynomialRing(F, "x")
    modulus = x^l - 1
    QR = QuotientRing(Fx, modulus)

    pk = bacuov_public_key(F, V, O, l)

    print "inside keygen with randomness", binascii.hexlify(randomness)

    # expand randomness into seeds
    output = SHAKE256(randomness, 2*SECURITY_LEVEL/4)

    # grab seed for S and for PCT
    seed_S = output[0:(SECURITY_LEVEL/4)]
    pk.seed_PCT = output[(SECURITY_LEVEL/4):]

    print "PCT seed:", binascii.hexlify(pk.seed_PCT);

    # sample secret linear transform S
    VS = MatrixSpace(F, l, 1)
    S_top_right = bacuov_generate_S(QR, V, O, seed_S)
    S = block_matrix([[MatrixSpace(F, V, V).identity_matrix(), S_top_right], [MatrixSpace(F, O, V).zero(), MatrixSpace(F, O, O).identity_matrix()]])
    S = block_matrix([[matrixify(S[i,j]) for j in range(0, N)] for i in range(0, N)])

    PP = [copy(MatrixSpace(F, n, n).zero()) for k in range(0,m)]
    PPc = [copy(MatrixSpace(QR, N, N).zero()) for k in range(0,m)]
    FFc = [copy(MatrixSpace(QR, N, N).zero()) for k in range(0,m)]
    # loop over all m polynomials
    for i in range(0,m):
        # generate top left and top right blocks of Pi
        Pi0V0V = bacuov_generate_vinegar_coefficients(QR, V, pk.seed_PCT, i)
        Pi0VVN = bacuov_generate_linear_coefficients(QR, V, O, pk.seed_PCT, i)
        Pi = block_matrix([[Pi0V0V, Pi0VVN], [Pi0VVN.transpose(), MatrixSpace(QR, O, O).zero()]])
        PPc[i] = Pi
        PP[i] = block_matrix([[matrixify(Pi[k,j]) for j in range(0, N)] for k in range(0, N)])
        FFc[i][0:V, 0:V] = PPc[i][0:V, 0:V]
        for j in range(0, V):
            for k in range(0, V):
                FFc[i][j,k] = flip(QR.gen()*FFc[i][j,k])

        FFc[i][0:V, V:N] = PPc[i][0:V, 0:V] * S_top_right
        for j in range(0,V):
            for k in range(V,N):
                FFc[i][j,k] = flip(QR.gen()^2 * FFc[i][j,k])

        mat = PPc[i][0:V, V:N] * QR.gen()
        for j in range(0, V):
            for k in range(0, O):
                mat[j,k] = flip(mat[j,k])

        FFc[i][0:V, V:N] = -FFc[i][0:V, V:N] + mat
        FFc[i][V:N, 0:V] = FFc[i][0:V, V:N].transpose()

        mat = FFc[i][0:N, 0:N]
        for j in range(0, N):
            for k in range(0, N):
                mat[j,k] = flip(mat[j,k]) * QR.gen()

        temp = mat[V:N, 0:V] * S_top_right
        for j in range(0, temp.nrows()):
            for k in range(0, temp.ncols()):
                temp[j,k] = temp[j,k] / QR.gen()

        #PPc[i][V:N, V:N] = S_top_right.transpose() * mat[0:V, 0:V] * S_top_right + mat[V:N, 0:V] * S_top_right + S_top_right.transpose() * mat[0:V, V:N]
        PPc[i][V:N, V:N] = S_top_right.transpose() * mat[0:V, 0:V] * S_top_right + temp.transpose() + temp
    

    ## explicit construction of redundant quadratic forms
    ##
    ### shorthands
    ##S_ = S[0:v, v:n]
    ##J = copy(MatrixSpace(F, l, l).zero())
    ##for i in range(0, l):
    ##    J[i,l-1-i] = 1
    ##Jo = copy(MatrixSpace(F, o, o).zero())
    ##Jv = copy(MatrixSpace(F, v, v).zero())
    ##JvS_Jo = Jv * S_ * Jo
    ##for k in range(0, O):
    ##    Jo[k*l:(k+1)*l, k*l:(k+1)*l] = J
    ##for k in range(0, V):
    ##    Jv[k*l:(k+1)*l, k*l:(k+1)*l] = J

    ### infer secret quadratic forms and remaining part of public
    ##FF = [copy(MatrixSpace(F, n, n).zero()) for i in range(0,m)]
    ##for k in range(0, len(FF)):
    ##    FF[k][0:v, 0:v] = Jv * PP[k][0:v, 0:v] * Jv
    ##    FF[k][0:v, v:n] = -Jv * PP[k][0:v, 0:v] * Jv*S_*Jo + Jv * PP[k][0:v, v:n] * Jo
    ##    FF[k][v:n, 0:v] = -Jo * S_.transpose() * Jv * PP[k][0:v, 0:v] * Jv + Jo * PP[k][0:v, v:n].transpose() * Jv
    ##    PP[k][v:n, v:n] = S_.transpose() * FF[k][0:v, 0:v] * S_ + Jo * FF[k][v:n, 0:v] * S_ + S_.transpose() * FF[k][0:v, v:n] * Jo

    ##for k in range(0, len(FF)):
    ##    mat = block_matrix([[matrixify(FFc[k][i,j]) for j in range(0, N)] for i in range(0, N)])
    ##    if mat != FF[k]:
    ##        print "error at k =", k
    ##        print "FFc[k]:"
    ##        print FFc[k]
    ##        print "whereas FF[k]:"
    ##        print FF[k]
    ##        return

    ##    mat = block_matrix([[matrixify(PPc[k][i,j]) for j in range(0, N)] for i in range(0, N)])
    ##    if mat != PP[k]:
    ##        print "error at k =", k
    ##        print "PPc[k]:"
    ##        print PPc[k]
    ##        print "whereas PP[k]:"
    ##        print PP[k]
    ##        return


    ##print "suucceeeessss \\o/"

    # compress public polynomials -- take out every block's first row
    oiloil = []
    for k in range(0, len(PPc)):
        for i in range(V, N):
            if i % l == 0:
                for j in range(i, N):
                    oiloil.append(PPc[k][i,j])

    return (FFc, S_top_right), (seed, oiloil, l, m, n)

def sign( E, sk, doc ):
    FF, S = sk
    n = S.nrows()
    m = len(FF)
    o = m
    v = n-o

    # get hash
    target = hash_to_vector(E, m, doc)

    # find vinegar variables that leads to invertible coefficient matrix
    is_invertible = False
    while not is_invertible:
    	xv = MatrixSpace(E, v, 1).random_element()
    	coefficient_matrix = copy(MatrixSpace(E, o, o).zero())
    	for i in range(0, o):
    	    coefficient_matrix[i,:] = xv.transpose() * (FF[i][0:v,v:(v+o)] + FF[i][v:(v+o),0:v].transpose())
        if coefficient_matrix.determinant() != 0:
            is_invertible = True

    # solve linear system to obtain matching oil variables
    target -= matrix([[(xv.transpose() * FF[i][0:v, 0:v] * xv)[0,0]] for i in range(0, m)])
    xo = coefficient_matrix.inverse() * target
    x = block_matrix([[xv], [xo]])

    # invert secret linear transform
    s = MatrixSpace(E, S.nrows(), S.ncols())(S.inverse()) * x

    return s

def verify( pk, doc, sig ):
    seed, oiloil, l, m, n = pk
    field = oiloil[0].parent()
    N = n / l
    O = m / l
    V = N - O
    v = V * l
 
    # prg all but lower-right of public quadratic rows
    PP = [copy(MatrixSpace(field, n, n).zero()) for k in range(0,m)]
    for k in range(0, len(PP)):
        for i in range(0, v+l-1):
            if i % l == 0:
                for j in range(i, n):
                    PP[k][i,j] = hash_to_element(field, seed + "||" + "pk" + str(k) + "||" + str(i) + "||" + str(j))
            else:
                for j in range(i - (i%l), n):
                    if j % l == l-1:
                        PP[k][i,j] = PP[k][i-1,j-l+1]
                    else:
                        PP[k][i,j] = PP[k][i-1,j+1]
        for i in range(0, n):
            for j in range(0, i):
                PP[k][i,j] = PP[k][j,i]

    # copy oil-oil first-rows
    z = 0
    for k in range(0, m):
        for i in range(v, n):
            if i % l == 0:
                for j in range(i, n):
                    PP[k][i,j] = oiloil[z]
                    z += 1

    # infer rest of block from first-rows
    for k in range(0, m):
        for I in range(V, N):
            for J in range(I, N):
                for i in range(1, l):
                    PP[k][I*l+i, J*l + l - 1] = PP[k][I*l+i-1, J*l]
                    PP[k][I*l+i, (J*l):(J*l+l-1)] = PP[k][I*l+i-1, (J*l+1):(J*l+l)]

    # make them symmetric
    for k in range(0, m):
        for i in range(v, n):
            for j in range(0, i):
                PP[k][i,j] = PP[k][j,i]

    # evaluate public polynomials
    E = sig[0,0].parent()
    target = hash_to_vector(E, len(PP), doc)
    compare = matrix([[(sig.transpose() * PP[i] * sig)[0,0]] for i in range(0, len(PP))])
    if target == compare:
        return True
    return False

