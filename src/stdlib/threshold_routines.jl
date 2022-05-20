"""
Returns the maximum FT coefficient of an IN array.
Thanks to Dr Abma for sharing his codes.
"""
function get_max_ft_coeff(IN::AbstractArray)
    ft_coeff = maximum( pmap(get_max_abs,IN) );
    return ft_coeff
end

function get_max_abs(IN::AbstractArray{T}) where {T <: Number}
    coeff = maximum( abs.( fft(IN) ) )
    return coeff
end

"""
Define threshold schedule.
""" 
function thresh_sched(IN::AbstractArray,K::Int,Pi::Real,Pf::Real,flag::String)

    # get max FT coeff among patches
    @show cmax = get_max_ft_coeff(IN)

    # output vector of coefficients
    out = zeros(eltype(cmax),K)

    # Upper and lower limits
    Up  = (Pi/100)*cmax
    Low = (Pf/100)*cmax

    # Define "scheduler" based on flag
    sched = if flag == "linear"

        function _linear(k::Int,K::Int,Up::Real,Low::Real)
            return Up + (Low - Up)*(k-1)/(K-1)
        end
        _linear

        elseif flag == "exp"

        function _expo(k::Int,K::Int,Up::Real,Low::Real)
            return Up*exp( ( log(Low/Up) )*(k-1)/(K-1) )
        end
        _expo

        elseif flag == "abma"

        function _abma(k::Int,K::Int,Up::Real,Low::Real)
            return Up*((K-k)/K)^2
        end
        _abma
    end

    @inbounds begin
        for k in eachindex(out)
            out[k] = sched(k,K,Up,Low)
        end
    end
    return out
end


"""
Define threshold schedule.
""" 
function _sched(IN::Array{T,N},K::Int,Pi::Real,Pf::Real,flag::String) where {T,N}

    # get max FT coeff among patches
    cmax = get_max_abs(IN)

    # output vector of coefficients
    out = Array{eltype(cmax)}(undef,K)

    # Upper and lower limits
    Up  = (Pi/100)*cmax
    Low = (Pf/100)*cmax

    # Define "scheduler" based on flag
    sched = if flag == "linear"

        function _linear(k::Int,K::Int,Up::Real,Low::Real)
            return Up + (Low - Up)*(k-1)/(K-1)
        end
        _linear

        elseif flag == "exp"

        function _expo(k::Int,K::Int,Up::Real,Low::Real)
            return Up*exp( ( log(Low/Up) )*(k-1)/(K-1) )
        end
        _expo

        elseif flag == "abma"

        function _abma(k::Int,K::Int,Up::Real,Low::Real)
            return Up*((K-k)/K)^2
        end
        _abma
    end

    @inbounds begin
        for k in eachindex(out)
            out[k] = sched(k,K,Up,Low)
        end
    end
    return out
end

"""
Define threshold schedule given upper and lower bounds.
""" 
function _schedule(Up::Real,Low::Real,K::Int,flag::String) where {T,N}

    # output vector of coefficients
    out = Array{eltype(Up)}(undef,K)

    # Define "scheduler" based on flag
    sched = if flag == "linear"

        function _linear(k::Int,K::Int,Up::Real,Low::Real)
            return Up + (Low - Up)*(k-1)/(K-1)
        end
        _linear

        elseif flag == "exp"

        function _expo(k::Int,K::Int,Up::Real,Low::Real)
            return Up*exp( ( log(Low/Up) )*(k-1)/(K-1) )
        end
        _expo

        elseif flag == "abma"

        function _abma(k::Int,K::Int,Up::Real,Low::Real)
            return Up*((K-k)/K)^2
        end
        _abma
    end

    @inbounds begin
        for k in eachindex(out)
            out[k] = sched(k,K,Up,Low)
        end
    end
    return out
end
