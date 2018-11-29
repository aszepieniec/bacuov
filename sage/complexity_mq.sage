###
# HS
# Computes the Hilbert series of a semi-regular homogeneous quadratic
# system of n variables and m equations.
#
def HS( m, n ):
    Fz.<z> = PowerSeriesRing(QQ, "z", 100)
    acc = Fz(1)
    for i in range(0, m):
        acc = acc * (1-z^2)

    return PolynomialRing(QQ, "z")(acc / (1-z)^n)#

# HS2
# Computes the Hilbert series of a semi-regular homogeneous quadratic
# system of n variables and m equations over GF(2).
#
def HS2( m, n ):
    Fz.<z> = PowerSeriesRing(QQ, "z", 100)
    acc = Fz(1)
    for i in range(0, n):
        acc = acc * (1+z)

    return PolynomialRing(QQ, "z")(acc / (1+z^2)^m)

##
# dreg
# Computes the degree of regularity of a semi-regular homogeneous
# quadratic system of n variables and m equations.
#
def dreg( q, m, n ):
    if q == 2:
        Snm = HS2(m, n)
    else:
        Snm = HS(m,n)
    coeffs = Snm.coefficients(sparse=False)
    for i in range(0, len(coeffs)):
        if coeffs[i] <= 0:
            return i

    # if we get here, the next coefficient is zero
    return Snm.degree()+1

##
# tau
# Computes the number of coefficients of a quadratic system.
#
def tau_( n ):
    return n * (n+1) / 2

##
# Tee
# Computes the number of monomials at a certain degree.
# 
def Tee_( n, Dreg ):
    return binomial(n+Dreg, Dreg)

##
# solving_complexity
# Computes the log2( ) of the complexity of solving a semi-regular
# homogeneous quadratic system using sparse linear algebra
# techniques.
#
def solving_complexity( q, m, n ):
    Dreg = dreg(q, m, n)
    t = tau_(n)
    T = Tee_(n, Dreg)
    c = log(1.0*t*T*T, 2.0)
    return c

##
# argmin
# please come built into the next version of python
#
def argmin( array ):
    minarray = min(array)
    for k in range(0, len(array)):
        if array[k] == minarray:
            return k

##
# classical_hybrid
# Computes the complexity of the classical hybrid attack for a range
# values for k
#
def classical_hybrid( q, m, n ):
    N = n
    M = m
    if N > m:
        alpha = 1.0 * N / M
        N = m - floor(alpha) + 1
        M = N
    
    complexity = []
    for k in range(0, N):
        complexity.append(solving_complexity(q, M, N-k) + 1.0*k*log(1.0*q,2.0))

    k = argmin(complexity)

    return complexity[k], k

##
# quantum_hybrid
# Computes the complexity of the quantum hybrid attack for a range
# values for k
#
def quantum_hybrid( q, m, n ):
    N = n
    M = m
    if N > m:
        alpha = 1.0 * N / M
        N = m - floor(alpha) + 1
        M = N

    if M <= 1:
        return 0, 0
    
    complexity = []
    for k in range(0, N):
        complexity.append(solving_complexity(q, M, N-k) + 0.5*k*log(1.0*q,2.0))

    k = argmin(complexity)

    return complexity[k], k
 
