module arr

    import ..Utils: na
    using Infiltrator

    export permutevec
    function permutevec(V::AbstractVector, d::Int)
        ind = Vector{Union{Vector{<:CartesianIndex},Colon}}(undef, d)
        [setindex!(ind, Utils.na, i) for i in 1:(d-1)]
        ind[d] = Colon()
        V[ind...]
    end

    function atleast2d(X::Array)
        if ndims(X) == 1
            X[:, na]
        else
            X
        end
    end

    function atleast3d(X::Array)
        if ndims(X) == 1
            X[:, na, na]
        elseif ndims(X) == 2
            X[:, :, na]
        else
            X
        end
    end

end
