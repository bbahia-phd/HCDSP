function LANCSSAOp(IN,k)

    # Hankelize
    H = HankelOp(IN);
    
    # Rank reduction
    U, Bk, V = HCDSP.lanbpro(H,k)
    
    # Averaging
    OUT = AveragingOp(U*U'*H,size(IN))

    return OUT
end