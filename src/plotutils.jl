module plotutils

    import ..DIutils as Utils
    using Plots
    using Dates
    using ColorSchemes
    using Infiltrator
    import DIutils.statistic: stat_quant
    using Statistics

    export  plotcolor

    #    _  _      ____      _   
    #  _| || |_   / ___| ___| |_ 
    # |_  ..  _| | |  _ / _ \ __|
    # |_      _| | |_| |  __/ |_ 
    #   |_||_|    \____|\___|\__|
    #                            
    # UTILTIES

    export getplotcolor
    function getplotcolor(propvec::AbstractArray, cmap::Symbol=:vik;
                         alpha=:none, alpha_min=0.1, alpha_max=1.0)
        getplotcolor(propvec, colorschemes[cmap]; alpha, alpha_min, alpha_max)
    end
    function getplotcolor(propvec::AbstractArray, cmap::ColorScheme;
                         alpha=:none, alpha_min=0.1, alpha_max=1.0)
        nans   = isnan.(propvec)
        valids = (!).(nans)
        propvec = Float64.(propvec)
        # propvec[valids] = Utils.norm_extrema(propvec[valids])
        results = fill(RGB(NaN,NaN,NaN), size(propvec))
        results[valids] = get(cmap, propvec[valids], :extrema)
        if alpha == :lowtohigh
            results[valids] = [RGBA(c.r, c.b, c.b, p) for (c,p) 
                in zip(results[valids], propvec[valids])]
        elseif alpha == :hightolow
            results[valids] = [RGBA(c.r, c.b, c.b, 1-p) for (c,p) 
                in zip(results[valids], propvec[valids])]
        elseif alpha == :none
        else
            throw(ArgumentError("alpha must be :lowtohigh, :hightolow, or :none"))
        end

        results
    end

    export lower_stat_quant, upper_stat_quant
    """
        lower_stat_quant(stat::Function, x, q, N=1000)
    obtains the relative lower bound of the quantile of a statistic
    """
    function lower_stat_quant(stat::Function, x, q, N=1000)
        stat(skipmissing(x)) - stat_quant(stat, x, q, N) 
    end
    """
        upper_stat_quant(stat::Function, x, q, N=1000)
    obtains the relative upper bound of the quantile of a statistic
    """
    function upper_stat_quant(stat::Function, x, q, N=1000)
        stat_quant(stat, x, q, N) - stat(skipmissing(x)) 
    end

    """
        get_ylims

    grabs the ylims for each subplot in the entire plot
    and returns rows (subplots) x columns (lower, upper)
    """
    function get_ylims(x::Plots.Plot)::Matrix
        response=[]
        for sub in x.subplots
            if sub isa Plots.Plot
                [push!(response, s) for s in sub.subplots]
            elseif sub isa Plots.Subplot
                push!(response, get_ylims(sub))
            end
        end
        hcat(collect.(response)...)'
    end
    get_ylims(x::Plots.Subplot) = ylims(x)
    get_ylims(x) = @infiltrate # If we don't recognize the type, infiltrate


    export get_lims
    """
        get_lims
    grabs the quantiles for each row in the matrix
    """
    function get_lims(X::Union{Vector,Matrix}, q=0.01)
        X = X isa Vector ? Matrix(X[:, Utils.na]) : X
        lims = []
        for i in axes(X, 1)
            push!(lims, quantile(X[i,:], [q, 1-q]))
        end
        [(q...,) for q in lims]
    end

    #    _  _     ____       _    
    #  _| || |_  / ___|  ___| |_  
    # |_  ..  _| \___ \ / _ \ __| 
    # |_      _|  ___) |  __/ |_  
    #   |_||_|   |____/ \___|\__| 
    #                             
    # UTILTIES

    function set_theme_timebased(time::Union{Real,Nothing}=nothing)
        time = time === nothing ? hour(now()) : time
        if Utils.in_range(time, [0,5]) ||
           Utils.in_range(time, [20, 24])
            Plots.theme(:brigh)
            theme="dark"
        else
            Plots.theme(:bright)
            theme="bright"
        end
    end
    
    function set_subplot_attr!(x::Plots.Plot, value, name)
        for sub in x.subplots
            set_subplot_attr!(x, value, name)
        end
    end
    function set_subplot_attr!(x::Plots.Subplot, value, name)
        x.attr[name] = value
    end

    function set_series_attr!(x::Plots.Plot, value, name)
        for sub in x.subplots
            set_subplot_attr!(x, value, name)
        end
    end
    function set_series_attr!(x::Plots.Subplot, value, name)
        for ser in x.series_list
            set_series_attr!(x, value, name)
        end
    end
    function set_series_attr!(x::Plots.Series, value, name)
        x.plotattributes[name] = value
    end

end
