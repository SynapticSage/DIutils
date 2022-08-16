module arr

    import ..Utils
    using Infiltrator

    export permutevec
    function permutevec(V::AbstractVector, d::Int)
        ind = Vector{Union{Vector{<:CartesianIndex},Colon}}(undef, d)
        [setindex!(ind, Utils.na, i) for i in 1:(d-1)]
        ind[d] = Colon()
        V[ind...]
    end

end
