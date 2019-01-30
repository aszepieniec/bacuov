import binascii
from CompactFIPS202 import SHAKE256, SHA3_256

def hash_digest( m ):
    return binascii.hexlify(SHA3_256(bytearray('1'+m)))

def sample_circulant_ring_element( R, buff ):
    Fx = R.modulus().parent()
    x = Fx.gen()
    F = Fx.base_ring()
    q = F.order()
    d = R.modulus().degree()
    integer = sum(256^i * buff[i] for i in range(0, len(buff)))
    expansion = []
    while integer > 0:
        expansion.append(integer % q)
        integer = floor(integer / q)
    poly = sum(Fx(expansion[i]) * x^i for i in range(0, len(expansion)))
    return R(poly)

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

    integer = sum(256^i * buff[i] for i in range(0, len(buff)))

    expansion = []
    while integer > 0:
        expansion.append(integer % p)
        integer = floor(integer / p)

    poly = sum(x^i * F(expansion[i]) for i in range(0, min(len(expansion), d)))
    return QR(poly)

def bacuov_generate_S( QR, V, O, seed ):
    l = QR.modulus().degree()
    # append indicator
    for i in range(0, 8):
        seed += bytearray([int(0 >> (8*i))])

    # expand
    buff = SHAKE256(seed, l * V * O)

    S = matrix([[bacuov_sample_circulant_ring_element(QR, buff[(l*(i*O+j)):(l*(i*O+j+1))]) for j in range(0, O)] for j in range(0, V)])
    return S

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
    S_top_right = block_matrix([[matrixify(S_top_right[i,j]) for j in range(0, O)] for i in range(0, V)])
    S = block_matrix([[MatrixSpace(F, v, v).zero(), S_top_right], [MatrixSpace(F, o, v).zero(), MatrixSpace(F, o, o).zero()]])
    for i in range(0, N):
	    for j in range(0, l):
                S[i*l + l - 1 - j, i*l + j] = 1

    # sample seed
    seed = str(bytearray([ZZ(Integers(256).random_element()) for i in range(0, 32)]))

    # prg all but lower-right of public quadratic rows
    PP = [copy(MatrixSpace(F, n, n).zero()) for k in range(0,m)]
    for k in range(0, len(PP)):
        for i in range(0, v+l-1):
            if i % l == 0:
                for j in range(i, n):
                    PP[k][i,j] = hash_to_element(F, seed + "||" + "pk" + str(k) + "||" + str(i) + "||" + str(j))
            else:
                for j in range(i - (i%l), n):
                    if j % l == l-1:
                        PP[k][i,j] = PP[k][i-1,j-l+1]
                    else:
                        PP[k][i,j] = PP[k][i-1,j+1]
        for i in range(0, n):
            for j in range(0, i):
                PP[k][i,j] = PP[k][j,i]

    # shorthands
    S_ = S[0:v, v:n]
    J = copy(MatrixSpace(F, l, l).zero())
    for i in range(0, l):
        J[i,l-1-i] = 1
    Jo = copy(MatrixSpace(F, o, o).zero())
    Jv = copy(MatrixSpace(F, v, v).zero())
    JvS_Jo = Jv * S_ * Jo
    for k in range(0, O):
        Jo[k*l:(k+1)*l, k*l:(k+1)*l] = J
    for k in range(0, V):
        Jv[k*l:(k+1)*l, k*l:(k+1)*l] = J

    # infer secret quadratic forms and remaining part of public
    FF = [copy(MatrixSpace(F, n, n).zero()) for i in range(0,m)]
    for k in range(0, len(FF)):
        FF[k][0:v, 0:v] = Jv * PP[k][0:v, 0:v] * Jv
        FF[k][0:v, v:n] = -Jv * PP[k][0:v, 0:v] * Jv*S_*Jo + Jv * PP[k][0:v, v:n] * Jo
        FF[k][v:n, 0:v] = -Jo * S_.transpose() * Jv * PP[k][0:v, 0:v] * Jv + Jo * PP[k][0:v, v:n].transpose() * Jv
        PP[k][v:n, v:n] = S_.transpose() * FF[k][0:v, 0:v] * S_ + Jo * FF[k][v:n, 0:v] * S_ + S_.transpose() * FF[k][0:v, v:n] * Jo

    # compress public polynomials -- take out every block's first row
    oiloil = []
    for k in range(0, len(PP)):
        for i in range(v, n):
            if i % l == 0:
                for j in range(i, n):
                    oiloil.append(PP[k][i,j])

    return (FF, S), (seed, oiloil, l, m, n)

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

