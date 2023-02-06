module plotutils

    import ..Utils
    using Plots
    using Dates
    using ColorSchemes
    using Infiltrator

    export  plotcolor

    #    _  _      ____      _   
    #  _| || |_   / ___| ___| |_ 
    # |_  ..  _| | |  _ / _ \ __|
    # |_      _| | |_| |  __/ |_ 
    #   |_||_|    \____|\___|\__|
    #                            
    # UTILTIES

    export getplotcolor
    function getplotcolor(propvec::AbstractVector, cmap::Symbol)
        propvec = Utils.norm_extrema(propvec)
        get.([colorschemes[cmap]], propvec)
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
