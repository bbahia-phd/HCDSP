import TSVD: tsvd
export tsvd

function complex_adjoint(Q::AbstractArray{Quaternion{T}}) where T

    # Obtain symplex and perplex
    Xs = complex.(real.(Q),  imagi.(Q));
    Xp = complex.(imagj.(Q), imagk.(Q));

    return vcat(hcat(Xs,Xp),hcat(-1 .* conj.(Xp), conj(Xs)))
end

function tsvd(Q::AbstractArray{Quaternion{T}},k) where T

    Qc = complex_adjoint(Q);

    return TSVD.tsvd(Qc,2k)
end
