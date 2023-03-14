module arr

    import ..DIutils: na
    import ..DIutils
    using Infiltrator, DataFrames, DimensionalData, LazyGrids

    export permutevec
    function permutevec(V::AbstractVector, d::Int)
        ind = Vector{Union{Vector{<:CartesianIndex},Colon}}(undef, d)
        [setindex!(ind, DIutils.na, i) for i in 1:(d-1)]
        ind[d] = Colon()
        V[ind...]
    end

    function atleast2d(X::AbstractArray)
        if ndims(X) == 1
            X[:, na]
        else
            X
        end
    end

    function atleast3d(X::AbstractArray)
        if ndims(X) == 1
            X[:, na, na]
        elseif ndims(X) == 2
            X[:, :, na]
        else
            X
        end
    end

    """
    convert dimarray into dataframe
    """
    function DataFrames.DataFrame(D::DimArray;
            name::Symbol=:value)::DataFrame
        n = Symbol.(DimensionalData.name.(D.dims))
        dims   = vec.(ndgrid(vec.(Array.(D.dims))...))
        value = vec(Array(D.data))
        DataFrame([dims..., value], [n..., name])
    end

    function DimensionalData.DimArray(D::DataFrame;name::Symbol=:value)::DimArray
        dimnames = setdiff(propertynames(D), [name])
        if !issorted(D, dimnames)
            D = sort(D, dimnames)
        end
        Dg = groupby(D, dimnames)
        dimvals = [unique(D[!,d]) for d in dimnames]
        diminds = [1:length(d) for d  in dimvals]
        @assert dimnames == Dg.cols
        vals = fill(NaN, length.(dimvals)...)
        iter = zip(Iterators.product(diminds...), Iterators.product(dimvals...))
        Threads.@threads for (inds, elem) in collect(iter)
            nt = NamedTuple(zip(dimnames, elem))
            d = Dg[nt]
            if :value in propertynames(d)
                vals[inds...] = d.value[1]
            end
        end
        dims = Tuple([DimensionalData.Dim{dimname}(dim) 
                      for (dim, dimname) in zip(dimvals, dimnames)])
        DimArray(vals, dims)
    end


end
