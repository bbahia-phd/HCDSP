function imputation_op(inp::AbstractArray{Ti,N},
                       Op::Function, args...;
                       iter = 10,
                       α = range(1,stop=0,length=iter),
                       ε=1e-16) where {Ti,Ts,N}

    # get sampling
    smp = SamplingOp(inp)

    # old or new
    return imputation_op(inp,smp,Op,args...; iter=iter, α=α, ε=ε)

end

function imputation_op(inp::AbstractArray{Ti,N},
                       smp::Union{AbstractArray{Ts,N},Function},
                       Op::Function, args...;
                       iter = 10,
                       α = range(1,stop=0,length=iter),
                       ε=1e-16) where {Ti,Ts,N}


    # allocate
    old,new,tmp = copy(inp),zero(inp),zero(inp);


    # safe-guard scalar α
    α = α .* ones(Ti,iter);

    # iterative reinsertion
    for i in 1:iter

        #@show i

        # apply operator
        tmp .= Op(old,args...)

        # reinsert
        new .= α[i] .* inp .+ (one(Ti) .- α[i] .* smp) .* tmp

        # Convergence
        if norm(new .- old,2)^2/norm(new,2)^2 < ε
            break;
        end
        
        old .= copy(new);
    end

    # old or new
    return new

end

function ImputationOp(IN::Array{Complex{Float64}},PARAM::Dict{String,Any})
    
    # Prelim
    Op = PARAM["Op"]
  iter = PARAM["iter"]
     N = PARAM["size"]
    haskey(PARAM,"alpha") ? a = PARAM["alpha"].*ones(iter,1) : a = range(1, stop = 0, length = iter)
    haskey(PARAM,"samp")  ? S = PARAM["samp"]                : S = vec( GetSamplingOp( reshape(IN,N) ) )
    OLD = copy(IN);
    NEW = zeros(eltype(IN),size(IN))
    OUT = zeros(eltype(IN),size(IN))
    
    # Iterate
    for i in 1:iter
        TMP = Op(OLD);
        
        NEW = a[i].*IN .+ (1 .- a[i].*S).*TMP
        
        # Convergence
        if norm(vec(NEW-OLD),2)^2/norm(vec(NEW),2)^2 < 1e-16
            break;
        end
        
        OLD = copy(NEW);
    end
    
    OUT = NEW;
    
    return vec(OUT)
    
end

function ImputationOp(IN::Array{Float64},PARAM::Dict{String,Any})
    
    # Prelim
    Op = PARAM["Op"]
  iter = PARAM["iter"]
     N = PARAM["size"]
    haskey(PARAM,"alpha") ? a = PARAM["alpha"].*ones(iter,1) : a = range(1, stop = 0, length = iter)
    haskey(PARAM,"samp")  ? S = PARAM["samp"]                : S = vec( GetSamplingOp( reshape(IN,N) ) )
    OLD = copy(IN);
    NEW = zeros(eltype(IN),size(IN))
    OUT = zeros(eltype(IN),size(IN))
    
    # Iterate
    for i in 1:iter
        TMP = Op(OLD);
        
        NEW = a[i].*IN .+ (1 .- a[i].*S).*TMP
        
        # Convergence
        if norm(vec(NEW-OLD),2)^2/norm(vec(NEW),2)^2 < 1e-16
            break;
        end
        
        OLD = copy(NEW);
    end
    
    OUT = NEW;
    
    return vec(OUT)
    
end