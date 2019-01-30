import sys
load("bacuov.sage")
from CompactFIPS202 import SHAKE256
import binascii

def test( num_trials, seed ):

    SECURITY_LEVEL = 256
    key_seed = SHAKE256(seed, SECURITY_LEVEL/4)
    print "set key seed:", binascii.hexlify(key_seed)

    F = FiniteField(7)
    V = 7
    O = 2
    l = 3

    num_successes = 0
    num_integrity_failures = 0
    num_decoding_failures = 0
    for trial_index in range(0, num_trials):
        key_seed = SHAKE256(key_seed, SECURITY_LEVEL/4)
        sk, pk = bacuov_keygen(SECURITY_LEVEL, F, V, O, l, key_seed)

    print "Ran", num_trials, "trials with", num_successes, "successes and", (num_integrity_failures + num_decoding_failures), "failures."
    print "Failures:"
    print " *", num_decoding_failures, "decoding errors"
    print " *", num_integrity_failures, "integrity errors"
    print "Successes:"
    print " *", num_successes, "total successes"

if len(sys.argv) != 3 or len(sys.argv[2]) % 2 != 0:
    print "usage: sage test_bacuov [num trials, eg 13] [random seed in hex, eg d13d13deadbeef]"
else:
    arg2 = bytearray(sys.argv[2].decode('hex'))
    test(int(sys.argv[1]), bytearray(arg2))

