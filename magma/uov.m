function AntiCirculantKeyGen( field, V, O, l )
    // infer arguments
    v := V * l;
    o := O * l;
    n := v + o;
    N := V + O;

    // sample S
    S := Matrix(field, [[Zero(field) : i in [1..n]] : j in [1..n]]);
    for i in [1..v] do
        for j in [(v+1)..n] do
            S[i,j] := Random(field);
        end for;
    end for;
    for I in [0..(N-1)] do
        for j in [1..l] do
            S[I*l+j, I*l+1+l-j] := 1;
        end for;
    end for;

    // make anti-circulant
    if l ne 1 then
        for i in [2..v] do
            if i mod l ne 1 then
                for j in [(v+1)..n] do
                    if j mod l eq 0 then
                        S[i,j] := S[i-1, j-l+1];
                    else
                        S[i,j] := S[i-1, j+1];
                    end if;
                end for;
            end if;
        end for;
    end if;

    // sample FF
    FF := [Matrix(field, [[Random(field) : j in [1..n]] : i in [1..n]]) : k in [1..o]];
    for k in [1..o] do
        // set oil-oil coefficients to zero
        for i in [(v+1)..n] do
            for j in [(v+1)..n] do
                FF[k][i,j] := 0;
            end for;
        end for;
        // make the rest anti-circulant
        if l ne 1 then
            for i in [2..n] do
                if i mod l ne 1 then
                    for j in [1..n] do
                        if j mod l eq 0 then
                            FF[k][i,j] := FF[k][i-1, j-l+1];
                        else
                            FF[k][i,j] := FF[k][i-1, j+1];
                        end if;
                    end for;
                end if;
            end for;
        end if;
    end for;
    // make upper triangular ...
    if Characteristic(field) eq 2 then
        for k in [1..o] do
            for i in [1..n] do
                for j in [1..i] do
                    FF[k][j,i] := FF[k][j,i] + FF[k][i,j];
                    FF[k][i,j] := 0;
                end for;
            end for;
        end for;
    // or symmetric
    else
        for k in [1..o] do
            for i in [1..n] do
                for j in [1..i] do
                    FF[k][i,j] := FF[k][j,i];
                end for;
            end for;
        end for;
    end if;

    // compute PP
    PP := [Transpose(S) * FF[k] * S : k in [1..o]];

    // if necessary, zero bottom halfs
    if Characteristic(field) eq 2 then
        for k in [1..o] do
            for i in [1..n] do
                for j in [1..(i-1)] do
                    PP[k][j,i] := PP[k][j,i] + PP[k][i,j];
                    PP[k][i,j] := 0;
                end for;
            end for;
        end for;
    end if;

    return <S, FF>, PP;

end function;

function KeyGen( field, v, o )
    n := v + o;
    S := IdentityMatrix(field, n);
    for i in [1..v] do
        for j in [(v+1)..n] do
            S[i,j] := Random(field);
        end for;
    end for;

    // sample F
    FF := [Matrix(field, [[Random(field) : j in [1..n]] : i in [1..n]]) : k in [1..o]];
    for k in [1..o] do
        for i in [(v+1)..n] do
            for j in [(v+1)..n] do
                FF[k][i,j] := 0;
            end for;
        end for;
    end for;
    if Characteristic(field) eq 2 then
        for k in [1..o] do
            for i in [1..n] do
                for j in [1..i] do
                    FF[k][j,i] := FF[k][j,i] + FF[k][i,j];
                    FF[k][i,j] := 0;
                end for;
            end for;
        end for;
    else
        for k in [1..o] do
            for i in [1..n] do
                for j in [1..i] do
                    FF[k][i,j] := FF[k][j,i];
                end for;
            end for;
        end for;
    end if;

    // compute PP
    PP := [Transpose(S) * FF[k] * S : k in [1..o]];

    // if necessary, zero bottom halfs
    if Characteristic(field) eq 2 then
        for k in [1..o] do
            for i in [1..n] do
                for j in [1..i] do
                    PP[k][j,i] := PP[k][j,i] + PP[k][i,j];
                    PP[k][i,j] := 0;
                end for;
            end for;
        end for;
    end if;

    return <S, FF>, PP;

end function;

function Sign( sk, msg )
    S := sk[1];
    FF := sk[2];
    n := NumberOfColumns(S);
    o := #FF;
    v := n-o;
    field := Parent(S[1,1]);

   
    is_invertible := false;
    while is_invertible eq false do
        xv := Matrix(field, [[Random(field)] : k in [1..v]]);
        coefficient_matrix := Matrix(field, [Eltseq(Transpose(xv) * (Submatrix(FF[k], 1, v+1, v, o) + Transpose(Submatrix(FF[k], v+1, 1, o, v)))) : k in [1..o]]);
        if Determinant(coefficient_matrix) ne 0 then
            is_invertible := true;
        end if;
    end while;
    target_vector := msg - Matrix(field, [Eltseq(Transpose(xv) * Submatrix(FF[k], 1, 1, v, v) * xv) : k in [1..o]]);
    xo := coefficient_matrix^-1 * target_vector;
    x := Matrix(field, [[xi] : xi in Eltseq(xv) cat Eltseq(xo)]);

    for k in [1..o] do
        Transpose(x) * FF[k] * x;
    end for;

    s := S^-1 * x;
    return s;
end function;

function Verify( pk, msg, sig )
    PP := pk;
    field := Parent(PP[1][1,1]);
    m := #PP;

    for k in [1..m] do
        e := (Transpose(sig) * PP[k] * sig)[1,1];
        e;
        if e ne msg[k,1] then
            return false;
        end if;
    end for;

    return true;
end function;

