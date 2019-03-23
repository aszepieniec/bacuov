load("complexity_mq.sage")

def params_plain( q, v, o ):
    m = o
    n = o + v
    pk = ceil(log(1.0*q, 2.0)) * (n*(n+1)/2) * m
    sig = ceil(log(1.0*q, 2.0)) * n

    print "|pk| = ", pk, "bits =", (1.0*pk/8), "bytes =", (1.0*pk/8/1024), "kilobytes"
    print "|sig| = ", sig, "bits =", (1.0*sig/8), "bytes =", (1.0*sig/8/1024), "kilobytes"

    print "Kipnis-Shamir:", log(1.0 * q^(v-o), 2.0) / 2
    c, k = quantum_hybrid(q, m, n)
    print "Direct Algebraic:", c, "with k =", k
    c, k = quantum_hybrid(q, m, v)
    print "UOV Reconciliation:", c, "with k =", k

def params_pct( q, v, o ):
    m = o
    n = o + v
    pk = ceil(log(1.0*q, 2.0)) * (o*(o+1)/2) * m
    sig = ceil(log(1.0*q, 2.0)) * n

    print "|pk| = ", pk, "bits =", (1.0*pk/8), "bytes =", (1.0*pk/8/1024), "kilobytes"
    print "|sig| = ", sig, "bits =", (1.0*sig/8), "bytes =", (1.0*sig/8/1024), "kilobytes"

    print "Kipnis-Shamir:", log(1.0 * q^(v-o), 2.0) / 2
    c, k = quantum_hybrid(q, m, n)
    print "Direct Algebraic:", c, "with k =", k
    c, k = quantum_hybrid(q, m, v)
    print "UOV Reconciliation:", c, "with k =", k

def params_luov( q, v, o, r ):
    m = o
    n = o + v
    pk = ceil(log(1.0*q, 2.0)) * (o*(o+1)/2) * m
    sig = ceil(log(1.0*q, 2.0))*r * n

    print "|pk| = ", pk, "bits =", (1.0*pk/8), "bytes =", (1.0*pk/8/1024), "kilobytes"
    print "|sig| = ", sig, "bits =", (1.0*sig/8), "bytes =", (1.0*sig/8/1024), "kilobytes"

    print "Kipnis-Shamir:", log(1.0 * q^(v-o), 2.0) / 2
    c, k = quantum_hybrid(q^r, m, n)
    print "Direct Algebraic:", c, "with k =", k
    # pessimistic
    #c, k = quantum_hybrid(q, m, v)
    #print "UOV Reconciliation:", c, "with k =", k
    # optimistic
    c, k = quantum_hybrid(q, v, v)
    print "UOV Reconciliation:", c, "with k =", k

def params_bac( q, V, O, l, r ):
    v = V * l
    o = O * l
    m = o
    n = o + v
    pk = ceil(log(1.0*q, 2.0)) * (O*(O+1)/2) * l * m
    sig = ceil(log(1.0*q, 2.0))*r * n


    c, k = quantum_hybrid(q^r, m, n)
    direct_algebraic = c

    c, k = quantum_hybrid(q, max(V,m), V)

    Fx.<x> = PolynomialRing(FiniteField(q), "x")
    factorization = (x^l - 1).factor()
    maxdeg = max(f.degree() for (f, m) in factorization)

    
    c, k = quantum_hybrid(q^maxdeg, max(V,m), V)
    uov_reconciliation = c

    kipnis_shamir = log(1.0 * q^(maxdeg*(V-O)), 2.0) / 2

    if verbose == True:
        print "|pk| = ", pk, "bits =", (1.0*pk/8), "bytes =", (1.0*pk/8/1024), "kilobytes"
        print "|sig| = ", sig, "bits =", (1.0*sig/8), "bytes =", (1.0*sig/8/1024), "kilobytes"
        print "Direct Algebraic:", direct_algebraic
        print "Kipnis-Shamir (conservative):", log(1.0 * q^(V-O), 2.0) / 2
        c, k = quantum_hybrid(q, max(V,m), V)
        print "UOV Reconciliation (aggressive):", c
        print "factorization:", factorization
        print "max degree:", maxdeg
        print "UOV Reconciliation (aggressive):", uov_reconciliation
        print "Kipnis-Shamir (aggressive):", kipnis_shamir

    return pk, sig, direct_algebraic, uov_reconciliation, kipnis_shamir

def sweep_bac(qq, VV, OO, ll, rr, sec_lvl):
    params_list = []
    
    for q in qq:
        print "q:", q
        for V in VV:
            print "V:", V
            for O in OO:
                print "O:", O
                for l in ll:
                    print "l:", l
                    for r in rr:
                        params = params_bac(q, V, O, l, r)
                        if params[2] > sec_lvl and params[3] > sec_lvl and params[4] > sec_lvl:
                            params_list.append((q, V, O, l, r, params))

    sorted_list = sorted(params_list, key=lambda item : (item[5][0] + item[5][1]))

    for elm in sorted_list[0:10]:
        print elm, "(", 1.0*(elm[5][0] + elm[5][1])/8/1024, "kB)"


