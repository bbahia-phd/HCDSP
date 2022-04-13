## Forward modeling Hankelization operator
function HankelOp(IN::AbstractArray{T,M};
                  N = size(IN),
                  L = Int.(floor.(N ./ 2) .+ 1),
                  K = Int.(N .- L .+ 1)) where {T <: Number,M}
    
    OUT = similar(IN,prod(L),prod(K));

    # Allocate
    if M == 1
        OUT = Hankelize(IN)       
    else
        AUX = similar(IN,L[1],K[1],N[2]);
        
        # Level one hankel matrices
        @inbounds for i in 1:N[2]
            AUX[:,:,i] .= Hankelize(IN[:,i])
        end
            
        # MBH
        @inbounds for j in 1:L[2], i in 1:K[2]
            n = i+j-1;
            ind1 = 1+(j-1)*L[1]:L[1]+(j-1)*L[1];
            ind2 = 1+(i-1)*K[1]:K[1]+(i-1)*K[1];
            OUT[ind1,ind2] = AUX[:,:,n]
        end 
    end
    
    return OUT
end
######################################################
function Hankelize(IN::AbstractArray{T};
                   L::Int = Int(floor(length(IN)/2))+1,
                   K::Int = length(IN)-L+1) where {T <: Number}

    OUT = similar(IN,L,K);
    
    @inbounds for i in 1:K
        ind = 1+i-1:L+i-1
        OUT[:,i] .= IN[ind]
    end
   
    return OUT
end
######################################################
## Pseudo-inverse of the Hankelization Operator Anti-diagonal averaging
function AveragingOp(IN::AbstractArray{T},
                     N::NTuple{M,Int};
                     L = Int.(floor.(N ./ 2) .+ 1),
                     K = Int.(N .- L .+ 1)) where {T <: Number, M}
    
    # Allocate
    OUT = zeros(T,N);
    
    if M == 1
#        count = zeros(Int,prod(L),prod(K))
          OUT .= AntiDiagAvg(IN,N[1])       
    else
        count = zeros(Int,L[1],K[1],N[2])
          AUX = zeros(T,L[1],K[1],N[2]);
        
        # Level one hankel matrices
        @inbounds for j in 1:L[2], i in 1:K[2]
            n = i+j-1;
            ind1 = 1+(j-1)*L[1]:L[1]+(j-1)*L[1];
            ind2 = 1+(i-1)*K[1]:K[1]+(i-1)*K[1];            

            count[:,:,n] .+= ones(L[1],K[1]);            
              AUX[:,:,n] .+= IN[ind1,ind2]
        end
        AUX ./= count
        
        @inbounds for j in 1:N[2]
            OUT[:,j] .= AntiDiagAvg(AUX[:,:,j],N[1])
        end
        
    end
    
    return OUT
end
######################################################
function AntiDiagAvg(IN::AbstractArray{T},N::Int) where {T <: Number}

    # Loop bounds
    l,k = size(IN);    
   
    # Allocate
      out = zeros(T,N);
    count = zeros(Int,N);
    
    # Antidiagonal averaging
    @inbounds for i in 1:l, j in 1:k
        count[i+j-1] += 1;
          out[i+j-1] += IN[i,j];        
    end
            
    out ./= count;
    return out
end

