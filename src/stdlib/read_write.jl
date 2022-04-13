"""
    read_write(file_name, flag; n=(), input=[])

Read or write binary file `file_name` according to `flag`.

If `flag="w"` then `input` is an array of size `n`.
If `flag="r"` then `input` is empty but `n` still represents array size.
# Arguments
-`file_name::AbstractString`: path to a binary file
-`flag::AbstractString`: "r" for read and "w" for write
-`n::Tuple`: tuple indicating array size
-`input::AbstractArray`: Array of size `n`

# Examples
```julia-repl
julia>
julia>
julia>
```
"""
function read_write(filename::AbstractString, flag::AbstractString ; n::Tuple{Vararg{Int}}=(), input=[])
    
    filename="$filename"
    fid = open(filename,flag)
    
    if flag == "r"
        
        input = zeros(Float32,prod(n));
        read!(fid, input);
        close(fid)

        return input;
        
    else flag == "w"
        
        write(fid,convert.(Float32,input))
        close(fid)

        return nothing
    end
end
