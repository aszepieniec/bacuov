function AttackAntiCirculant( field, V, O, l )
    sk, pk := AntiCirculantKeyGen( field, V, O, l );
    Fx := PolynomialRing(field, #pk);
    xx := [Fx.i : i in [1..#pk]];

    n := (V+O)*l;

    // insert solution
    s := Matrix(field, [[Random(field)] : i in [1..n]]);
    target := Matrix(field, [[(Transpose(s) * pk[k] * s)[1,1]] : k in [1..#pk]]);

    fx := Matrix(Fx, [[Fx!s[i,1]] : i in [1..n]]);
    for i in [1..#pk] do
        fx[i,1] := xx[i];
    end for;

    equations := [];
    for k in [1..#pk] do
        equations := equations cat [((Transpose(fx) * Matrix(Fx, pk[k]) * fx)[1,1] - target[k,1])];
    end for;

    // add field equations
    for k in [1..#pk] do
        equations := equations cat [xx[k]^Characteristic(field) - xx[k]];
    end for;

    tick := Realtime();
    gb := GroebnerBasis(equations);
    tock := Realtime();

    return tock-tick;
end function;

function AttackRegular( field, v, o )
    sk, pk := KeyGen( field, v, o );
    Fx := PolynomialRing(field, #pk);
    xx := [Fx.i : i in [1..#pk]];

    n := v + o;

    // insert solution
    s := Matrix(field, [[Random(field)] : i in [1..n]]);
    target := Matrix(field, [[(Transpose(s) * pk[k] * s)[1,1]] : k in [1..#pk]]);

    fx := Matrix(Fx, [[Fx!s[i,1]] : i in [1..n]]);
    for i in [1..#pk] do
        fx[i,1] := xx[i];
    end for;

    equations := [];
    for k in [1..#pk] do
        equations := equations cat [(Transpose(fx) * Matrix(Fx, pk[k]) * fx)[1,1] - target[k,1]];
    end for;

    // add field equations
    for k in [1..#pk] do
        equations := equations cat [xx[k]^Characteristic(field) - xx[k]];
    end for;

    tick := Realtime();
    gb := GroebnerBasis(equations);
    tock := Realtime();

    return tock-tick;
end function;

function Attack( q, m )
    o := m;
    v := 2*o;
    tick := Realtime();
    if q lt 256 then
        print "l = 1";
        tim := AttackRegular(FiniteField(q), v, o);
        PrintFile("timings.dat", "q =\t" cat Sprint(q) cat "\tv =\t" cat Sprint(v) cat "\to =\t" cat Sprint(o) cat "\tl =\t1\ttime =\t" cat Sprint(tim));
    end if;
    for l in [2..m] do
        if o mod l eq 0 then
            print "l = ", l;
            tim := AttackAntiCirculant(FiniteField(q), v/l, o/l, l);
            PrintFile("timings.dat", "q =\t" cat Sprint(q) cat "\tv =\t" cat Sprint(v) cat "\to =\t" cat Sprint(o) cat "\tl =\t" cat Sprint(l) cat "\ttime =\t" cat Sprint(tim));
        end if;
    end for;
    tock := Realtime();

    return tock-tick;
end function;

tuples := <<2, 21>, <2, 18>, <2, 12>, <2, 8>, <7, 21>, <7, 18>, <7, 12>, <7, 8>, <7, 6>, <251, 21>, <251, 18>, <251, 12>, <251, 8>, <251, 6>>;

