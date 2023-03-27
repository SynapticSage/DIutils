module statistic

    export  pfunc
    using StatsBase
    using Infiltrator
    import ..DIutils

    function pfunc(x::Real)
        if x < 1e-3
            "***"
        elseif x < 1e-2
            "**"
        elseif x < 0.05
            "*"
        elseif x < 0.1
            "†"
        else
            ""  
        end
    end
    pfunc(x::Vector{<:Real}) = pfunc.(x)

    using Statistics, Bootstrap
    """
    TITLE: MeanBoot
    Purpose: Shortcut for mean bootstrap
    """
    meanboot(data, N) = confint(bootstrap(mean, data, BasicSampling(N)),
                                BasicConfInt(0.95));
    meanboot(data, N, ci) = confint(bootstrap(mean, data, BasicSampling(N)),
                                BasicConfInt(ci));
    boot(data, N, func) = confint(bootstrap(func, data, BasicSampling(N)),
                                BasicConfInt(0.95));

    export stat_quant
    """
    TITLE: stat_quant
    Purpose: Shortcut for quantile bootstrap of a statistic

    # Arguments
    - `stat::Function`: the statistic to bootstrap
    - `x`: the data to bootstrap
    - `q`: the quantile to bootstrap
    - `N`: the number of bootstrap samples
    # Returns
    - the quantile of the statistic
    """
    function stat_quant(stat::Function, x::AbstractArray, 
         q::Union{<:Real,<:AbstractVector{<:Real}}, N=1000)
        b=bootstrap(stat, collect(skipmissing(x)), Bootstrap.BasicSampling(N))
        # quantile(skipmissing(x), 0.975) - nanmean(skipmissing(x))
        quantile(b.t1[1], q)
    end

    """
    TITLE: bootGroups
    Purpose: shorcut for groupby bootstrapping means
    """
    function bootGroups(x, groups; centerByFields=[], field=:occNorm, mb=meanboot)
        if length(centerByFields) > 0
        end
        xg = groupby(x, groups);
        xr = combine(xg, field => (x -> mb(x,100)) => [:mean, :lower, :upper]);
        return xr
    end

    """
    TITLE: bootGroups
    Purpose: shorcut for groupby bootstrapping means
    """
    function bootstrapArray(x, func; N, dims)
    end
    """
    TITLE: bootGroups
    Purpose: shorcut for groupby bootstrapping means
    """
    function ciArray(x, func; N, dims)
    end

    """
        dummycode
    
    (create documentation here)
    """
    function dummycode(X::AbstractMatrix)
        hcat([dummycode(x) for x in eachcol(X)]...)
    end
    function dummycode(X::AbstractVector)
        UX = sort(unique(X))
        X = X .== UX'
    end

end
