module clean
    
    import ..Utils
    using Infiltrator
    using NaNStatistics, Statistics

    """
    quantile_filter_dims
    """
    function quantile_filter_dims(X, q; dims=2)
        inds = inds_quantile_filter_dims(X,q; dims)
        X[inds, :]
    end
    function inds_quantile_filter_dims(X, q; dims=2)
        if ndims(X) > 2; @warn "careful: not designed for ndims>2 right now"; end

        slices = eachslice(X; dims)
        matches = []
        for slice in slices
            q1 = nanquantile(slice, q[1])
            q2 = nanquantile(slice, q[2])
            push!(matches, Utils.in_range(slice, [q1,q2]))
        end
        accumulate(.*, matches)[end] 
    end

end
