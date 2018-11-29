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

    print "|pk| = ", pk, "bits =", (1.0*pk/8), "bytes =", (1.0*pk/8/1024), "kilobytes"
    print "|sig| = ", sig, "bits =", (1.0*sig/8), "bytes =", (1.0*sig/8/1024), "kilobytes"

    print "Kipnis-Shamir:", log(1.0 * q^(V-O), 2.0) / 2
    c, k = quantum_hybrid(q^r, m, n)
    print "Direct Algebraic:", c, "with k =", k
    c, k = quantum_hybrid(q, max(V,m), V)
    print "UOV Reconciliation:", c, "with k =", k


