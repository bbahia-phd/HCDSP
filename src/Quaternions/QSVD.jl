import TSVD: tsvd
export tsvd, fwd_complex_adjoint, bck_complex_adjoint

function fwd_complex_adjoint(Q::AbstractArray{Quaternion{T}}) where T

    # Obtain symplex and perplex
    Xs = complex.(real.(Q),  imagi.(Q));
    Xp = complex.(imagj.(Q), imagk.(Q));

    return vcat(hcat(Xs,Xp),hcat(-1 .* conj.(Xp), conj(Xs)))
end

function bck_complex_adjoint(C::AbstractArray{Complex{T}}) where T
    n = Int.(size(C) ./ 2)
    return Quaternion.(C[1:n[1],1:n[2]],C[1:n[1],n[2]+1:end])
end