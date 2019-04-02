function AttackExtensionAntiCirculant( field, V, O, l, r )
    sk, pk := AntiCirculantKeyGen( field, V, O, l );

    extension_field := FiniteField((#field)^r);
    Fx := PolynomialRing(extension_field, #pk);
    xx := [Fx.i : i in [1..#pk]];

    n := (V+O)*l;

    // insert solution
    s := Matrix(extension_field, [[Random(extension_field)] : i in [1..n]]);
    target := Matrix(extension_field, [[(Transpose(s) * Matrix(extension_field, [[extension_field!pk[k][i,j] : j in [1..n]] : i in [1..n]]) * s)[1,1]] : k in [1..#pk]]);

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
        equations := equations cat [xx[k]^(#extension_field) - xx[k]];
    end for;

    tick := Realtime();
    gb := GroebnerBasis(equations);
    tock := Realtime();

    return tock-tick;
end function;

function AttackExtensionRegular( field, v, o, r )
    sk, pk := KeyGen( field, v, o );

    extension_field := FiniteField((#field)^r);
    Fx := PolynomialRing(extension_field, #pk);
    xx := [Fx.i : i in [1..#pk]];

    n := v + o;

    // insert solution
    s := Matrix(extension_field, [[Random(extension_field)] : i in [1..n]]);
    target := Matrix(extension_field, [[(Transpose(s) * Matrix(extension_field, [[extension_field!pk[k][i,j] : j in [1..n]] : i in [1..n]]) * s)[1,1]] : k in [1..#pk]]);

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
        equations := equations cat [xx[k]^(#extension_field) - xx[k]];
    end for;

    tick := Realtime();
    gb := GroebnerBasis(equations);
    tock := Realtime();

    return tock-tick;
end function;

function AttackExtension( q, m, r )
    o := m;
    v := 2*o;
    tick := Realtime();
    if q lt 256 then
        print "l = 1";
        tim := AttackExtensionRegular(FiniteField(q), v, o, r);
        PrintFile("timings_extension.dat", "q =\t" cat Sprint(q) cat "\tr = \t" cat Sprint(r) cat "\tv =\t" cat Sprint(v) cat "\to =\t" cat Sprint(o) cat "\tl =\t1\ttime =\t" cat Sprint(tim));
    end if;
    for l in [2..m] do
        if o mod l eq 0 then
            print "l =", l;
            tim := AttackExtensionAntiCirculant(FiniteField(q), v/l, o/l, l, r);
            PrintFile("timings_extension.dat", "q =\t" cat Sprint(q) cat "\tr = \t" cat Sprint(r) cat "\tv =\t" cat Sprint(v) cat "\to =\t" cat Sprint(o) cat "\tl =\t" cat Sprint(l) cat "\ttime =\t" cat Sprint(tim));
        end if;
    end for;
    tock := Realtime();

    return tock-tick;
end function;


