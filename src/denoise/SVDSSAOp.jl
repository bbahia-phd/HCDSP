#import TSVD: tsvd

# SVD-based rank-reduction
function rank_reduction(IN::AbstractArray{T},k::Int) where {T <: Number}
    U,_,_ = tsvd(IN,k)
    return U*U'*IN
end

# QSVD-based rank-reduction
function rank_reduction(IN::AbstractArray{Quaternion{T}},k::Int) where {T <: Number}
    Q = fwd_complex_adjoint(IN);
    U,_,_ = tsvd(Q,k)
    OUT = U*U'*Q
    return bck_complex_adjoint(OUT) 
end

function SVDSSAOp(IN::AbstractArray{T,N}, k::Int) where {T,N}
# IN is the frequency slice
# k is its expected rank
    
    # Hankelize
    H = HankelOp(IN);

    # Rank reduction
    H .= rank_reduction(H,k);

    # Averaging
    OUT = AveragingOp(H,size(IN))
   
   return OUT
end