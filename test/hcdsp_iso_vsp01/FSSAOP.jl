function FSSAOp(IN,k)

    # OUT
    out = zero(IN)

    # dimensions of array
    dims = size(IN)    

    # Matrix dimensions
    L = floor.(Int64, dims ./ 2) .+ 1;
    K = dims .- L .+ 1;

    fwd(x) = mbh_multiply(IN,x,flag="fwd");
    adj(x) = mbh_multiply(IN,x,flag="adj");

    A(x,i;kwargs...) = i == 1 ? fwd(x) : adj(x)

    U, Bk, V = HCDSP.lanbpro(A,k,m=prod(L),n=prod(K))

    Ub = U*Bk;

    # do fast anti-diagonal averaging using rank-1 approx
    for i in 1:k
        out += HCDSP.anti_diagonal_summation(Ub[:,i],V[:,i],L,K);
    end

    count = HCDSP.count_copy_times([100,100])

    return out ./ count
end

function QFSSAOp(IN,k)

    # OUT
    out = zero(IN)

    # dimensions of array
    dims = size(IN)    

    # Matrix dimensions
    L = floor.(Int64, dims ./ 2) .+ 1;
    K = dims .- L .+ 1;

    fwd(x) = qmbh_multiply(IN,x,flag="fwd");
    adj(x) = qmbh_multiply(IN,x,flag="adj");

    A(x,i;kwargs...) = i == 1 ? fwd(x) : adj(x)

    U, Bk, V = HCDSP.lanbpro(A,k,m=prod(L),n=prod(K),qflag=true)

    Ub = U*Bk;

    # do fast anti-diagonal averaging using rank-1 approx
    for i in 1:k
        out += anti_diagonal_summation(Ub[:,i],V[:,i],L,K);
    end

    count = HCDSP.count_copy_times([40,40])

    return out ./ count
end

afwd(IN,x) = vcat(
    qmbh_multiply(IN,x,flag="fwd"),
    qmbh_multiply(invi.(IN),x,flag="fwd"),
    qmbh_multiply(invj.(IN),x,flag="fwd"),
    qmbh_multiply(invk.(IN),x,flag="fwd") )


function aadj(IN,x)
    l = floor(Int64,length(x)/4)
    ind1 = 1:l;     xo  = @view x[ind1];
    ind2 = l+1:2l;  xi = @view x[ind2];
    ind3 = 2l+1:3l; xj = @view x[ind3]; 
    ind4 = 3l+1:4l; xk = @view x[ind4];

    out = qmbh_multiply(      IN, xo, flag="adj") .+ 
          qmbh_multiply(invi.(IN),xi, flag="adj") .+
          qmbh_multiply(invj.(IN),xj, flag="adj") .+
          qmbh_multiply(invk.(IN),xk, flag="adj")

    return out     
end

function AQFSSAOp(IN,k)

    # dimensions of array
    dims_in = size(IN);
    dims_pad = nextpow.(2,dims_in);

    # PAD IN
    INP = PadOp(IN,nin=dims_in,npad=dims_pad,flag="fwd");

    # OUT
    OUTP = zero(INP)

    # Matrix dimensions
    L = floor.(Int64, dims_pad ./ 2) .+ 1;
    K = dims_pad .- L .+ 1;

    A(x,i;kwargs...) = i == 1 ? afwd(INP,x) : aadj(INP,x)

    U, Bk, V = HCDSP.lanbpro(A,k,m=4prod(L),n=prod(K),qflag=true)

    Ub = U*Bk;

    # do fast anti-diagonal averaging using rank-1 approx
    for i in 1:k
        OUTP += anti_diagonal_summation(Ub[1:prod(L),i],V[:,i],L,K);
    end

    count = HCDSP.count_copy_times([dims_pad[1],dims_pad[2]])

    OUTP = OUTP ./ count;

    return PadOp(OUTP,nin=dims_in,npad=dims_pad,flag="adj");
end

function QRFSSAOp(IN::AbstractArray{Quaternion{T}},k) where T

    # dimensions of array
    dims = size(IN);

    # OUT
    OUT = zero(IN)

    # Matrix dimensions
    L = floor.(Int64, dims ./ 2) .+ 1;
    K = dims .- L .+ 1;

    # Random matrix for projection
    Ω = Quaternion.(rand(prod(K),k),rand(prod(K),k),rand(prod(K),k));

    # projection
    P = zerosq((prod(L),k))
    for i = 1:k
        P[:,i] .= qmbh_multiply(IN,Ω[:,i],flag="fwd")
    end

    # Quaternion QR
    Qr,_ = qr(P)

    for i in 1:k

        q = copy(Qr[:,i]);

        t = qmbh_multiply(IN,q,flag="adj");

        OUT .+= anti_diagonal_summation(q,t,L,K);        
    end

    count = HCDSP.count_copy_times([dims[1],dims[2]])

    OUT = OUT ./ count;

    return OUT

end


function QRFSSAOp(IN::AbstractArray{T},k) where T

    # dimensions of array
    dims = size(IN);

    # OUT
    OUT = zero(IN)

    # Matrix dimensions
    L = floor.(Int64, dims ./ 2) .+ 1;
    K = dims .- L .+ 1;

    # Random matrix for projection
    Ω = rand(prod(K),k);

    # projection
    P = zeros(T,(prod(L),k))
    for i = 1:k
        P[:,i] .= mbh_multiply(IN,Ω[:,i],flag="fwd")
    end

    # Quaternion QR
    Qr,_ = LinearAlgebra.qr(P)

    for i in 1:k

        q = copy(Qr[:,i]);

        t = mbh_multiply(IN,q,flag="adj");

        OUT .+= anti_diagonal_summation(q,t,L,K);        
    end

    count = HCDSP.count_copy_times([dims[1],dims[2]])

    OUT = OUT ./ count;

    return OUT

end