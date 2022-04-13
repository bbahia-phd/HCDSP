function prediction_quality(pred::AbstractArray,ideal::AbstractArray)

    α = dot(ideal,pred) / dot(pred,pred);
    Q = 10 * log( norm(ideal,2) / norm(ideal - α .* pred,2)  )

    return Q
end
