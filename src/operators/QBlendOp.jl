"""
    Operator to blend a 3D common-receiver gather (CRG).

    PARAM: Dictionary? I don't think so... maybe a NamedTuple.
    
    flag: "fwd" or "adj" operator 

    if flag == "fwd"        
         IN: (nt,nsx,nsy) array representing a CRG
        OUT: (ntb,1) array representing long blended record
    
    if flag == "adj"
         IN: (ntb,1) array representing long blended record
        OUT: (nt,nsx,nsy) array representing a pseudo-deblended CRG
    
    ```
    julia> b = SeisBlendOp(IN,PARAM,flag)
        
    ```
"""
function SeisBlendOp(IN::AbstractArray{T},PARAM::NamedTuple{(:nt, :nx, :ny, :dt, :tau, :sx, :sy),Tuple{Int, Int, Int, T, Vector{T},Vector{Int},Vector{Int}}},flag::String) where {T <: Number}
    
    # Get parameters
    tau = PARAM[:tau]  # firing times
     sx = PARAM[:sx]   # orderd sources indexes in x
     sy = PARAM[:sy]   # orderd sources indexes in y
     dt = PARAM[:dt]   # sampling in time
     nt = PARAM[:nt]   # time samples
     nx = PARAM[:nx]   # sources in x
     ny = PARAM[:ny]   # sources in y
    
    # Firing times to indexes
    j = floor.(Int64, tau ./ dt) .+ 1
    
    # Number of samples in blended record 
    ntb = maximum(j) + nt - 1;

    # Blending
    if flag == "fwd"     
        
        # Output allocation
        OUT = zeros(eltype(IN), ntb);
        
        @inbounds for is in eachindex(j)
            # Indexes 
            it, i1, i2 = 0, j[is], j[is]+nt-1;
            @inbounds for jj in i1:i2
                it += 1;
                OUT[jj] = OUT[jj] + IN[it,sx[is],sy[is]];
            end
        end

        # Pseudo-deblending
        elseif flag == "adj"
        
        # Output allocation
        OUT = zeros(eltype(IN), nt, nx, ny);
        
        @inbounds for is in eachindex(j)
            # Indexes 
            it, i1, i2 = 0, j[is], j[is]+nt-1;
            @inbounds for jj in i1:i2
                it += 1;            
                OUT[it,sx[is],sy[is]] = OUT[it,sx[is],sy[is]] + IN[jj];
            end            
        end

    end
    
    return OUT

end


function QBlendOp(IN::AbstractArray{Quaternion{T}},PARAM::NamedTuple{(:nt, :nx, :ny, :dt, :tau, :sx, :sy),Tuple{Int, Int, Int, T, Vector{T}, Vector{Int}, Vector{Int}}}, flag::String) where {T <: Number}

    # Get parameters                                                                                                                                                                                                                                                            
    tau = PARAM[:tau]  # firing times                                                                                                                                                                                                                                           
     sx = PARAM[:sx]   # orderd sources indexes in x                                                                                                                                                                                                                            
     sy = PARAM[:sy]   # orderd sources indexes in y                                                                                                                                                                                                                            
     dt = PARAM[:dt]   # sampling in time                                                                                                                                                                                                                                       
     nt = PARAM[:nt]   # time samples                                                                                                                                                                                                                                           
     nx = PARAM[:nx]   # sources in x                                                                                                                                                                                                                                           
     ny = PARAM[:ny]   # sources in y                                                                                                                                                                                                                                           

    # Firing times to indexes                                                                                                                                                                                                                                                   
    j = floor.(Int64, tau ./ dt) .+ 1

    # Number of samples in blended record                                                                                                                                                                                                                                       
    ntb = maximum(j) + nt - 1;

    # Blending                                                                                                                                                                                                                                                                  
    if flag == "fwd"

        # Output allocation                                                                                                                                                                                                                                                     
        OUT = zeros(eltype(IN), ntb);

        @inbounds for is in eachindex(j)
            # Indexes                                                                                                                                                                                                                                                           
            it, i1, i2 = 0, j[is], j[is]+nt-1;
            @inbounds for jj in i1:i2
                it += 1;
                OUT[jj] = OUT[jj] + IN[it,sx[is],sy[is]];
            end
        end

        # Pseudo-deblending                                                                                                                                                                                                                                                     
        elseif flag == "adj"

        # Output allocation                                                                                                                                                                                                                                                     
        OUT = zeros(eltype(IN), nt, nx, ny);

        @inbounds for is in eachindex(j)
            # Indexes                                                                                                                                                                                                                                                           
            it, i1, i2 = 0, j[is], j[is]+nt-1;
            @inbounds for jj in i1:i2
                it += 1;
                OUT[it,sx[is],sy[is]] = OUT[it,sx[is],sy[is]] + IN[jj];
            end
        end

    end

    return OUT

end

function QBlendOp2!(DATA::AbstractArray{Quaternion{T},1},MODEL::AbstractArray{Quaternion{T},3},PARAM::NamedTuple{(:nt, :nx, :ny, :dt, :tau, :sx, :sy),Tuple{Int, Int, Int, T, Vector{T}, Vector{Int}, Vector{Int}}}, flag::String) where {T <: Number}

    # Get parameters                                                                                                                                                                                                                                                            
    tau = PARAM[:tau]  # firing times                                                                                                                                                                                                                                           
        sx = PARAM[:sx]   # orderd sources indexes in x                                                                                                                                                                                                                            
        sy = PARAM[:sy]   # orderd sources indexes in y                                                                                                                                                                                                                            
        dt = PARAM[:dt]   # sampling in time                                                                                                                                                                                                                                       
        nt = PARAM[:nt]   # time samples                                                                                                                                                                                                                                           
        nx = PARAM[:nx]   # sources in x                                                                                                                                                                                                                                           
        ny = PARAM[:ny]   # sources in y                                                                                                                                                                                                                                           

    # Firing times to indexes                                                                                                                                                                                                                                                   
    j = floor.(Int64, tau ./ dt) .+ 1

    # Number of samples in blended record                                                                                                                                                                                                                                       
    ntb = maximum(j) + nt - 1;

    # Blending                                                                                                                                                                                                                                                                  
    if flag == "fwd"
        
        @inbounds for is in eachindex(j)
            # Indexes                                                                                                                                                                                                                                                           
            it, i1, i2 = 0, j[is], j[is]+nt-1;
            @inbounds for jj in i1:i2
                it += 1;
                DATA[jj] = DATA[jj] + MODEL[it,sx[is],sy[is]];
            end
        end

        # Pseudo-deblending                                                                                                                                                                                                                                                     
        elseif flag == "adj"
        
        @inbounds for is in eachindex(j)
            # Indexes                                                                                                                                                                                                                                                           
            it, i1, i2 = 0, j[is], j[is]+nt-1;
            @inbounds for jj in i1:i2
                it += 1;
                MODEL[it,sx[is],sy[is]] = MODEL[it,sx[is],sy[is]] + DATA[jj];
            end
        end

    end

end
        
