import sys
load("bacuov.sage")
from CompactFIPS202 import SHAKE256
import binascii

def test( num_trials, seed ):

    SECURITY_LEVEL = 256
    randomness = SHAKE256(seed, SECURITY_LEVEL/4)
    print "set key seed:", binascii.hexlify(randomness)

    F = FiniteField(7)
    V = 7
    O = 2
    l = 11
    r = 5

    Fx.<x> = PolynomialRing(F, "x")
    #poly = x^2
    #while not poly.is_irreducible():
    #    poly = x^5 + sum(x^i * F.random_element() for i in range(0, 5))
    poly = F(1) + F(5)*x + F(1)*x^2 + F(2)*x^3 + F(1)*x^4 + x^5
    E = FiniteField(7^5, "x", modulus=poly)

    num_successes = 0
    for trial_index in range(0, num_trials):
        randomness = SHAKE256(randomness[0:(SECURITY_LEVEL/4)], 2*SECURITY_LEVEL/4)
        print "random buffer:", binascii.hexlify(randomness)
        key_seed = randomness[0:(SECURITY_LEVEL/4)]
        document = randomness[(SECURITY_LEVEL/4):]
        sk, pk = bacuov_keygen(SECURITY_LEVEL, F, V, O, l, key_seed)
        sig = bacuov_sign(SECURITY_LEVEL, E, sk, document)

        if bacuov_verify(pk, document, sig) == True:
            num_successes += 1

    print "Ran", num_trials, "trials with", num_successes, "successes and", (num_trials - num_successes), "failures."
    print "Failures:"
    print " *", (num_trials-num_successes)
    print "Successes:"
    print " *", num_successes, "total successes"

if len(sys.argv) != 3 or len(sys.argv[2]) % 2 != 0:
    print "usage: sage test_bacuov [num trials, eg 13] [random seed in hex, eg d13d13deadbeef]"
else:
    arg2 = bytearray(sys.argv[2].decode('hex'))
    test(int(sys.argv[1]), bytearray(arg2))

