## rQR-based rank-reduction
function rqr(IN::AbstractArray{T},k::Int) where {T <: Number}

    # Dimensionality reduction
    n = size(IN);
    Ω = rand(n[2],k);
    P = IN*Ω;

    # QR
    Q,_ = qr(P)
    Q = Q[:,1:k];

    # "Rank" reduction
    OUT = Q*Q'*IN
end

function rQROp(IN::AbstractArray{T}, k::Int) where {N,T}
# IN is the frequency slice
# k  is its expected rank
  
   # Hankelize
   H = HankelOp(IN);
   
   # Rank reduction
   H .= rqr(H,k);
   
   # Averaging
   OUT = AveragingOp(H,size(IN))
   
   return OUT
end



