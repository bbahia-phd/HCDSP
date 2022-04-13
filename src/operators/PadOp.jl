"""
Operator to do N-dimensional padding and truncation.

IN is the time domain data to be padded or truncated
nin is the input dimension
nout is the output dimensions
flag = "fwd" for padding
flag = "adj" for trunaction
Passes the dot product test
"""
function PadOp(IN::AbstractArray{T}; nin::Tuple{Vararg{Int}}=(), npad::Tuple{Vararg{Int}}=(), flag::String) where {T <: Number}

    # Indexes to get
    inx = CartesianIndices(nin)
    
    if flag == "fwd"
        # Zero-padding: IN is of size nin
        OUT = zeros(eltype(IN),npad...)
        for i in inx
            OUT[i] = IN[i]
        end

        elseif flag == "adj"
        
        # Truncation: IN is of size npad
        OUT = zeros(eltype(IN),nin...)
        for i in inx
            OUT[i] = IN[i]
        end
    end

    return OUT
end