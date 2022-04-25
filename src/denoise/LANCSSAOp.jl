function LANCSSAOp(IN,k)

    # Hankelize
    H = HankelOp(IN);
    
    # Rank reduction
    U,_,_ = HCDSP.lanbpro(H,k)
    
    # Averaging
    OUT = AveragingOp(U*U'*H,size(IN))

    return OUT
end