module statistic

    export  pfunc

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
end
