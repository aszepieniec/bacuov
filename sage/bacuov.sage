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
    r = E.modulus().degree()
    F = E.modulus().parent().base_ring()
    z = E.gen()
    VS = MatrixSpace(E, d, 1)
    vec = copy(VS.zero())
    print "hashing:", binascii.hexlify(m)
    digest = SHAKE256(m, d*r)
    #digest = SHA3_256(m)
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

def bacuov_keygen( SECURITY_LEVEL, F, V, O, l, randomness ):
    o = O * l
    v = V * l
    N = O + V
    n = o + v
    m = o
    print "got m:", m

    Fx = PolynomialRing(F, "x")
    modulus = x^l - 1
    QR = QuotientRing(Fx, modulus)

    pk = bacuov_public_key(F, V, O, l)

    #print "inside keygen with randomness", binascii.hexlify(randomness)

    # expand randomness into seeds
    print "SHAKE input:", binascii.hexlify(randomness)
    output = SHAKE256(randomness, 2*SECURITY_LEVEL/4)
    print "SHAKE output:", binascii.hexlify(output)

    # grab seed for S and for PCT
    seed_S = output[0:(SECURITY_LEVEL/4)]
    pk.seed_PCT = output[(SECURITY_LEVEL/4):]

    print "PCT seed:", binascii.hexlify(pk.seed_PCT);

    # sample secret linear transform S
    VS = MatrixSpace(F, l, 1)
    S_top_right = bacuov_generate_S(QR, V, O, seed_S)
    #print "S:"
    #bacuov_print_circulant_ring_matrix(S_top_right)
    #input(1)
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
        #PP[i] = block_matrix([[matrixify(Pi[k,j]) for j in range(0, N)] for k in range(0, N)])

        #print "Pi0VVN:"
        #bacuov_print_circulant_ring_matrix(Pi0VVN)
        #input(1)

        # get FF[i][0:V, 0:V]
        FFc[i][0:V, 0:V] = PPc[i][0:V, 0:V]
        for j in range(0, V):
            for k in range(0, V):
                FFc[i][j,k] = flip(QR.gen()*FFc[i][j,k])

        # get FF[i][0:V, V:N] (and its counterpart)
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

        #print "FFov:"
        #bacuov_print_circulant_ring_matrix(FFc[i][0:V, V:N])
        #input(1)

        # get PP[i][V:N, V:N]
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

    #print "FF[0][0:V, 0:V]:"
    #bacuov_print_circulant_ring_matrix(FFc[0][0:V,0:V])
    #print "FF[0][0:V, V:N]:"
    #bacuov_print_circulant_ring_matrix(FFc[0][0:V,V:N])
    #print "PP[0][V:N, V:N]:"
    #bacuov_print_circulant_ring_matrix(PPc[0][V:N,V:N])

    # compress public polynomials -- take out every block's first row
    oiloil = []
    for k in range(0, len(PPc)):
        for i in range(V, N):
            if i % l == 0:
                for j in range(i, N):
                    oiloil.append(PPc[k][i,j])

    return (randomness, FFc, S_top_right), (seed, oiloil, l, m, n)

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
    randomness, FF, S_top_right = sk
    V = S_top_right.nrows()
    l = S_top_right[0,0].parent().modulus().degree()
    v = V*l
    m = len(FF)
    o = m
    O = o/l
    n = v+o
    r = E.modulus().degree()
    z = E.gen()

    # get hash
    target = hash_to_vector(E, m, doc)

    # find vinegar variables that lead to invertible coefficient matrix
    is_invertible = False
    for trial_index in range(0, 256):
        array = SHAKE256(doc, m*r)
        array += randomness
        array.append(trial_index)
        vinegar_array = SHAKE256(array, v*r)
    	xv = MatrixSpace(E, v, 1).random_element() # no.
        for i in range(0, v):
            xv[i,0] = sum(E(vinegar_array[i*r + j]) * z^j for j in range(0, r))
        #print "vinegar vector:"
        #print bacuov_str_gfpem(xv)
        
    	coefficient_matrix = copy(MatrixSpace(E, o, o).zero())
    	for i in range(0, o):
            FFi = block_matrix([[matrixify(FF[i][j,k]) for k in range(0, FF[i].ncols())] for j in range(0, FF[i].nrows())])
            #if i >= 1:
            #    print "FFvoi:"
            #    #print FF[i][0:V, V:(V+O)]
            #    print bacuov_str_gfpem(FFi[0:v, v:(v+o)])
            #    print "xv: "
            #    print bacuov_str_gfpem(xv.transpose())
            #    input(1)
    	    coefficient_matrix[i,:] = xv.transpose() * (FFi[0:v,v:(v+o)] + FFi[v:(v+o),0:v].transpose())

            #if i >= 1:
            #    print "line of coefficient matrix:"
            #    #print bacuov_str_gfpem(coefficient_matrix[i,:])
            #    print bacuov_str_gfpem(xv.transpose() * FFi[0:v, v:(v+o)])
            #    input(1)

        if coefficient_matrix.determinant() != 0:
            break

    # solve linear system to obtain matching oil variables
    b = target
    for i in range(0, m):
        FFi = block_matrix([[matrixify(FF[i][j,k]) for k in range(0, FF[i].ncols())] for j in range(0, FF[i].nrows())])
        b[i,0] = b[i,0] - (xv.transpose() * FFi[0:v, 0:v] * xv)[0,0]

    xo = coefficient_matrix.inverse() * b
    x = block_matrix([[xv], [xo]])

    # invert secret linear transform
    #s = MatrixSpace(E, S.nrows(), S.ncols())(S.inverse()) * x
    #S_ = block_matrix([[matrixify(S_top_right[i,j]) for j in range(0, S_top_right.ncols())] for i in range(0, S_top_right.nrows())])
    #S_xo = S_ * xo
    #s = x
    #s[0:v,0] = s[0:v,0] - S_xo
    S_full = block_matrix([[MatrixSpace(E, V, V).identity_matrix(), S_top_right], [MatrixSpace(E, O, V).zero(), MatrixSpace(E, O, O).identity_matrix()]])
    S_ = block_matrix([[matrixify(S_full[i,j]) for j in range(0, S_full.ncols())] for i in range(0, S_full.nrows())])
    Sinv = S_.inverse()
    s = Sinv * x

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

