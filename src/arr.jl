"""
    arr.jl

Module containing shorcuts for working with arrays and DimArrays (dimensioned
arrays)
"""
module arr

    import ..DIutils: na
    import ..DIutils
    using Infiltrator, DataFrames, DimensionalData, LazyGrids, Statistics,
          NaNStatistics

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

    export get_quantile_filtered
    """
        get_quantile_filtered
    filters out the quantiles of the data
    # Arguments
    - `X::Matrix`: the data to filter
    - `q::Real`: the quantile to filter
    """ 
    function get_quantile_filtered(X::Matrix, q::Real)
        if size(X,1) > size(X,2)
            trans = true
            x = X'
        else
            trans = false
        end
        lims = get_lims(x, q)
        inds = fill(true, size(x,2))
        for lim in lims
            inds = inds .& (x[:,lim[1]] .< lim[2])
        end
        @assert any(inds)
        trans ? x[:, findall(inds)] : x[:, findall(inds)]'
    end

    function get_quantile_filtered(X::AbstractVector, q::Real; 
        set_missing=false, set_zero=false, set_lim=false)
        q_lower = quantile(X, q)
        q_upper = quantile(X, 1-q)
        if set_missing
            X[(X .< q_lower) .& (X .> q_upper)] .= missing
            return X
        elseif set_zero
            X[(X .< q_lower) .& (X .> q_upper)] .= 0
            return X
        elseif set_lim
            X[(X .< q_lower)] .= q_lower
            X[(X .> q_upper)] .= q_upper
            return X
        else
            return X[(X .< q_lower) .& (X .> q_upper)]
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

    """
        detect_undef_elements
    detects the elements that are undef in an array
    # Arguments
    - `X::AbstractArray`: the array to detect
    # Returns
    - `inds::BitArray`: the indices of the elements that are undef
    # Examples
    ```julia
    julia> v=Vector{DataFrame}(undef, 3)
    3-element Vector{DataFrame}:
     #undef
     #undef
     #undef
    julia> detect_undef_elements(v)
    3-element BitArray{1}:
     1
     1
     1
    ```
    """
    function detect_undef_elements(X::AbstractArray)
        inds = fill(true, size(X))
        for i in eachindex(X)
           # Detect elements that are undef via UndefInitializer
           inds[i] = !isassigned(X, i)
        end
        inds
    end

    """
        undef_to_missing(X::AbstractArray)
    converts undef elements to missing
    # Arguments
    - `X::AbstractArray`: the array to convert
    # Returns
    - `X::AbstractArray`: the converted array
    """
    function undef_to_missing(X::AbstractArray)
        inds = detect_undef_elements(X)
        X = Missing <: eltype(X) ? X : 
            convert(Array{Union{Missing,eltype(X)}}, X)
        X[inds] .= missing
        X
    end


end
