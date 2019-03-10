import binascii
from CompactFIPS202 import SHAKE256, SHA3_256

def hash_digest( m ):
    return binascii.hexlify(SHA3_256(bytearray('1'+m)))

def hash_to_element( E, m ):
    hexmap = dict({'0': 0, '1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8, '9': 9, 'a': 10, 'b': 11, 'c': 12, 'd': 13, 'e': 14, 'f': 15})
    if E.is_prime_F():
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
    r = E.modulus().degree()
    F = E.modulus().parent().base_ring()
    z = E.gen()
    VS = MatrixSpace(E, d, 1)
    vec = copy(VS.zero())
    print "hashing:", binascii.hexlify(m)
    digest = SHAKE256(m, d*r)
    print "hash output:", binascii.hexlify(digest)
    array = [F(digest[i]) for i in range(0, d*r)]
    for i in range(0, d):
        vec[i,0] = sum(array[i*r + j] * z^j for j in range(0, r))
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
        self.F = F
        self.V = V
        self.O = O
        self.l = l
        self.N = V+O
        self.v = V*l
        self.o = O*l
        self.n = self.N*l
        self.m = self.O*l
        Fx.<x> = PolynomialRing(F, "x")
        divisor = x^l - 1
        self.QR = QuotientRing(Fx, divisor)
        self.FF = [matrix([[self.QR.zero() for j in range(0, V+O)] for i in range(0, V+O)]) for k in range(0, O*l)]
        self.S = matrix([[self.QR.zero() for j in range(0, O)] for i in range(0, V)])

class bacuov_public_key:
    def __init__( self, F, V, O, l ):
        self.F = F
        self.V = V
        self.O = O
        self.l = l
        self.N = V+O
        self.v = V*l
        self.o = O*l
        self.n = self.N*l
        self.m = self.O*l
        Fx.<x> = PolynomialRing(F, "x")
        divisor = x^l - 1
        self.QR = QuotientRing(Fx, divisor)
        self.PP = [matrix([[self.QR.zero() for j in range(0, O)] for i in range(0, O)]) for k in range(0, O*l)]
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

def matrixify_anticirculant( qr_element ):
    QR = qr_element.parent()
    modulus = QR.modulus()
    d = modulus.degree()
    Fx = modulus.parent()
    F = Fx.base_ring()
    MS = MatrixSpace(F, d, d)

    coeffs = qr_element.lift().coefficients(sparse=False)

    while len(coeffs) != d:
        coeffs.append(F(0))

    coeffs = [coeffs[d-1-i] for i in range(0,d)]

    mat = copy(MS.zero())
    for j in range(0, d):
        mat[0,j] = coeffs[j]
    for i in range(1, d):
        for j in range(0, d-1):
            mat[i,j] = mat[i-1,j+1];
        mat[i,d-1] = mat[i-1,0]

    return mat

def matrixify_circulant( qr_element ):
    QR = qr_element.parent()
    modulus = QR.modulus()
    d = modulus.degree()
    Fx = modulus.parent()
    F = Fx.base_ring()
    MS = MatrixSpace(F, d, d)

    coeffs = qr_element.lift().coefficients(sparse=False)

    while len(coeffs) != d:
        coeffs.append(F(0))

    mat = copy(MS.zero())
    for j in range(0, d):
        mat[0,j] = coeffs[j]
    for i in range(1, d):
        for j in range(1, d):
            mat[i,j] = mat[i-1,j-1];
        mat[i,0] = mat[i-1,d-1]

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
    F = QR.modulus().parent().base_ring()
    z = QR.gen()
    # append indicator
    for i in range(0, 4):
        seed += bytearray([int(0 >> (8*i))])

    # expand
    buff = SHAKE256(seed, l * V * O)
    #print "buffer for S:", binascii.hexlify(buff)

    S = matrix([[QR.zero() for j in range(0,O)] for j in range(0,V)])
    for i in range(0, V):
        for j in range(0, O):
            S[i,j] = sum(z^k * F(buff[(i*O + j) * l + k]) for k in range(0, l))

    return S

def bacuov_generate_vinegar_coefficients( QR, V, rand, index ):
    l = QR.modulus().degree()
    F = QR.modulus().parent().base_ring()
    z = QR.gen()
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
    K = 0;
    for i in range(0, V):
        Pi[i,i] = sum(F(buff[K + k]) * z^k for k in range(0, l))
        K += l
        for j in range(i+1, V):
            Pi[i,j] = sum(F(buff[K + k]) * z^k for k in range(0, l))
            Pi[j,i] = Pi[i,j]
            K += l

    return Pi
   
def bacuov_generate_linear_coefficients( QR, V, O, rand, index ): 
    l = QR.modulus().degree()
    F = QR.modulus().parent().base_ring()
    z = QR.gen()
    # append indicator
    randomness = copy(rand)
    for i in range(0, 4):
        randomness += bytearray([int(2 >> (8*i))])
    for i in range(0, 4):
        randomness += bytearray([int((index >> (8*i)) % 256)])

    # expand
    buff = SHAKE256(randomness, l * V * O)
    #print "buffer:", binascii.hexlify(buff)

    Pi = matrix([[QR.zero() for j in range(0,O)] for i in range(0,V)])
    for i in range(0, V):
        for j in range(0, O):
            Pi[i,j] = sum(F(buff[(i*O+j)*l + k]) * z^k for k in range(0, l))

    return Pi

def bacuov_circulant_smart_reverse( elm ):
    QR = elm.parent()
    coeffs = elm.lift().coefficients(sparse=False)
    while len(coeffs) != QR.modulus().degree():
        coeffs.append(elm.parent().zero())
    coeffs = [coeffs[0]] + [coeffs[len(coeffs)-1-i] for i in range(0, len(coeffs)-1)]
    return sum(QR.gen()^i * coeffs[i] for i in range(0, len(coeffs)))

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

    #print "inside keygen with randomness", binascii.hexlify(randomness)

    # expand randomness into seeds
    output = SHAKE256(randomness, 2*SECURITY_LEVEL/4)

    # grab seed for S and for PCT
    seed_S = output[0:(SECURITY_LEVEL/4)]
    pk.seed_PCT = output[(SECURITY_LEVEL/4):]

    # sample secret linear transform S
    VS = MatrixSpace(F, l, 1)
    S_top_right = bacuov_generate_S(QR, V, O, seed_S)
    #print "S:"
    #bacuov_print_circulant_ring_matrix(S_top_right)
    #input(1)
    S = block_matrix([[MatrixSpace(F, V, V).identity_matrix(), S_top_right], [MatrixSpace(F, O, V).zero(), MatrixSpace(F, O, O).identity_matrix()]])
    S = block_matrix([[matrixify_circulant(S[i,j]) for j in range(0, N)] for i in range(0, N)])

    PP = [copy(MatrixSpace(F, n, n).zero()) for k in range(0,m)]
    PPc = [copy(MatrixSpace(QR, N, N).zero()) for k in range(0,m)]
    FFc = [copy(MatrixSpace(QR, N, N).zero()) for k in range(0,m)]
    #S_top_right = matrix([[bacuov_circulant_smart_reverse(S_top_right[i,j]) for j in range(0, S_top_right.ncols())] for i in range(0, S_top_right.nrows())])
    # loop over all m polynomials
    for i in range(0,m):
        # generate top left and top right blocks of Pi
        Pi0V0V = bacuov_generate_vinegar_coefficients(QR, V, pk.seed_PCT, i)
        Pi0VVN = bacuov_generate_linear_coefficients(QR, V, O, pk.seed_PCT, i)
        Pi = block_matrix([[Pi0V0V, Pi0VVN], [Pi0VVN.transpose(), MatrixSpace(QR, O, O).zero()]])
        PPc[i] = Pi

        # get FF[i][0:V, 0:V]
        FFc[i][0:V, 0:V] = PPc[i][0:V, 0:V]

        # get FF[i][0:V, V:N] (and its counterpart)
        FFc[i][0:V, V:N] = PPc[i][0:V, V:N] - PPc[i][0:V, 0:V] * S_top_right
        FFc[i][V:N, 0:V] = FFc[i][0:V, V:N].transpose()

        # get PP[i][V:N, V:N]
        temp = FFc[i][V:N, 0:V] * S_top_right
        PPc[i][V:N, V:N] = S_top_right.transpose() * FFc[i][0:V, 0:V] * S_top_right + temp.transpose() + temp
   
    pk.PP = [PPc[i][V:N, V:N] for i in range(0,m)]

    sk = bacuov_secret_key(F, V, O, l) 
    sk.FF = [FFc[i] for i in range(0,m)]
    sk.S_top_right = S_top_right
    sk.randomness = randomness
    
    return sk, pk

def bacuov_str_gfpem( mat ):
    E = mat[0,0].parent()
    r = E.modulus().degree()
    F = E.modulus().parent().base_ring()
    string = "["
    for i in range(0, mat.nrows()):
        strings = []
        for j in range(0, mat.ncols()):
            coeffs = mat[i,j].polynomial().coefficients(sparse=False)
            while len(coeffs) != r:
                coeffs.append(F(0))
            strings.append("".join(str(coeffs[k]) for k in range(0, r)))
        string += "["
        string += ", ".join(strings)
        string += "]"
        if i != mat.nrows()-1:
            string += "\n"
    string += "]"
    return string

def gfpe_to_str( elm ):
    E = elm.parent()
    r = E.modulus().degree()
    F = E.modulus().parent().base_ring()
    coeffs = elm.polynomial().coefficients(sparse=False)
    while len(coeffs) != r:
        coeffs.append(F(0))
    return "".join(str(c) for c in coeffs)

def bacuov_sign( SECURITY_LEVEL, E, sk, doc ):
    FF = sk.FF
    S_top_right = matrix([[bacuov_circulant_smart_reverse(sk.S_top_right[i,j]) for j in range(0, sk.S_top_right.ncols())] for i in range(0, sk.S_top_right.nrows())])

    V = sk.V
    l = sk.l
    v = V*l
    m = sk.m
    o = m
    O = o/l
    N = V + O
    n = v+o
    r = E.modulus().degree()
    z = E.gen()

    # get randomness
    randomness = sk.randomness + SHA3_256(doc)

    # determine target from hash-of-document
    array = SHAKE256(doc, m*r)
    target = matrix([[sum(E.gen()^j * array[i*r+j] for j in range(0, r))] for i in range(0, m)])
    # target = hash_to_vector(E, m, doc)

    # find vinegar variables that lead to invertible coefficient matrix
    is_invertible = False
    for trial_index in range(0, 256):
        vinegar_array = SHAKE256(randomness + chr(trial_index), v*r)
    	xv = copy(MatrixSpace(E, v, 1).zero())
        for i in range(0, v):
            xv[i,0] = sum(E(vinegar_array[i*r + j]) * z^j for j in range(0, r))
        
    	coefficient_matrix = copy(MatrixSpace(E, o, o).zero())
    	for i in range(0, o):
            FFi = block_matrix([[matrixify_anticirculant(FF[i][j,k]) for k in range(0, FF[i].ncols())] for j in range(0, FF[i].nrows())])
    	    coefficient_matrix[i,:] = xv.transpose() * (FFi[0:v,v:(v+o)] + FFi[v:(v+o),0:v].transpose())

        if coefficient_matrix.determinant() != 0:
            break

    # solve linear system to obtain matching oil variables
    b = target
    for i in range(0, m):
        FFi = block_matrix([[matrixify_anticirculant(FF[i][j,k]) for k in range(0, FF[i].ncols())] for j in range(0, FF[i].nrows())])
        b[i,0] = b[i,0] - (xv.transpose() * FFi[0:v, 0:v] * xv)[0,0]

    xo = coefficient_matrix.inverse() * b
    x = block_matrix([[xv], [xo]])
    
    # invert secret linear transform
    QR = S_top_right[0,0].parent()
    Sinv_complete = block_matrix([[MatrixSpace(QR, V, V).identity_matrix(), -S_top_right], [MatrixSpace(QR, O, V).zero(), MatrixSpace(QR, O, O).identity_matrix()]])
                                                                          # ^ inverse
    Sinv_pressed = block_matrix([[matrixify_circulant(Sinv_complete[i,j]) for j in range(0, Sinv_complete.ncols())] for i in range(0, Sinv_complete.nrows())])
    s = Sinv_pressed * x

    return s

def bacuov_verify( pk, doc, sig ):
    F = pk.F
    n = pk.n
    m = pk.m
    l = pk.l
    QR = pk.QR
    E = sig[0,0].parent()
    r = E.modulus().degree()
    N = pk.N
    O = pk.O
    V = pk.V
    
    # get P matrices from PRG and public key, evaluate while at it
    evaluation = matrix([[E.zero()] for i in range(0,m)])
    for i in range(0,m):
        # generate top left and top right blocks of Pi
        Pi0V0V = bacuov_generate_vinegar_coefficients(QR, V, pk.seed_PCT, i)
        Pi0VVN = bacuov_generate_linear_coefficients(QR, V, O, pk.seed_PCT, i)
        Pi = block_matrix([[Pi0V0V, Pi0VVN], [Pi0VVN.transpose(), pk.PP[i]]])
        PPi = block_matrix([[matrixify_anticirculant(Pi[j,k]) for k in range(0, N)] for j in range(0, N)])
        evaluation[i,0] = (sig.transpose() * PPi * sig)[0,0]

    # get hash
    # target = hash_to_vector(E, m, doc)
    array = SHAKE256(doc, m*r)
    target = matrix([[sum(E.gen()^j * array[i*r+j] for j in range(0, r))] for i in range(0, m)])


    # compare evaluation to hash
    if target == evaluation:
        return True

    else:
        return False

