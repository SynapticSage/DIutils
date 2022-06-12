module utils

    import Random
    using CSV, DataFrames
    using Gadfly
    using Colors, ColorSchemes
    using Pushover
    using Statistics, NaNStatistics
    using Plots, DrWatson
    using ProgressMeter
    using Infiltrator
    using ThreadsX

    import SearchSortedNearest
    searchsortednearest = SearchSortedNearest.searchsortednearest
    searchsortednext    = SearchSortedNearest.searchsortednext

    export skipnan
    export itsizeof, piso
    export squeeze
    export searchsortednext, searchsortednearest
    export remove_key_item

    skipnan(x) = Iterators.filter(!isnan, x)
    na = [CartesianIndex()]

    """
    mkdifne

    like mkdir, except makes if not exist
    """
    function mkifne(path)
        if !(isdir(path)); mkdir(path); end
    end

    function itsizeof(X)
        [size(x) for x in X]
    end
    function piso(X::Union{Vector,Tuple})
        println([size(x) for x in X])
    end
    function piso(X::T) where T <: Dict
        println(Dict(key=>size(x) for (key,x) in X))
    end
    function pnf(X::Union{Vector,Tuple})
        println([mean(isnan.(x)) for x in X])
    end
    function pnf(X::T) where T <: Dict
        println(Dict(key=>mean(isnan.(x)) for (key,x) in X))
    end

    function norm_extrema(x::Vector{T1}, minmax::Union{Vector{T2},Tuple{T2, T2}}) where
        T1 <: Real where T2 <: Real
        if minmax isa Tuple
            minmax = [minmax...]
        end
        if minmax[2] == minmax[1]
            x = minmax[1] * ones(size(x))
        else
            x = (x .- minimum(x))./(maximum(x) - minimum(x))
            x = x .* diff(minmax) .+ minmax[1]
        end
    end

    function in_range(X::AbstractArray, range::Union{Tuple, Vector})
        X .≥ range[1] .&& X .< range[2]
    end

    function squeeze(A::AbstractArray)  
        s = size(A)
        A = dropdims(A, dims = tuple(findall(size(A) .== 1)...))
        return A
    end  

    function dextrema(A::AbstractArray; kws...)
        diff([nanminimum(A), nanmaximum(A)], kws...)
    end
    range_extrema = dextrema

    function randomize_int(X)
        Xmin = minimum(X);
        Xmax = maximum(X);
        initial = collect(Xmin:Xmax);
        final   = Random.shuffle(initial);
        mapping(x) = Dict(initial .=> final)[x]
        map(mapping, X)
    end

    """
    savef

    saves in multiple formats
    """
    function savef(args...;formats=["png","svg","pdf"])
        for format in formats
            @debug "format=$format"
            Plots.savefig(plotsdir(args[1:end-1]..., args[end]*".$format"))
        end
    end

    guides = Dict(); # shortcut for guides
    guides[(:x, :stopWell)] = Guide.xlabel("Goal")
    guides[(:y, :stopWell)] = Guide.ylabel("Goal")
    guides[(:x, :neuron)] = Guide.xlabel("Neuron")
    guides[(:y, :neuron)] = Guide.ylabel("Neuron")
    guides[(:x, :rayleighZ)] = Guide.xlabel("Rayleigh Ζ")
    guides[(:x, :rayleighZ_diff)] = Guide.xlabel("Rayleigh Z\nDifferences")
    guides[(:x, :gt_shuffle)] = Guide.xlabel("Percent\nreal > shuff") 

    """
    TITLE: goalVectorTheme
    Purpose: theme for goal vector shit
    """
    function goalVectorTheme()
        theme = Theme(major_label_color=colorant"white", major_label_font_size=14pt,
                      minor_label_color=colorant"white",
                      key_label_color=colorant"white",
                      key_title_color=colorant"white",
                      panel_fill=colorant"black",
                      background_color=colorant"black")
        Gadfly.push_theme(theme)
    end

    function getPushoverClient()
        token = open(expanduser("~/.pushover.token"),"r") do f
            token = read(f, String)
        end
        user = open(expanduser("~/.pushover.user"),"r") do f
            user = read(f, String)
        end
        return PushoverClient(user, token)
    end

    function pushover(pos...; kws...)
        send(getPushoverClient(), pos...; kws...)
    end


    function findgroups(pos...)
        X = Matrix(hcat(pos...))
        uX = unique(X, dims=1)
        g = zeros(Int,size(X,1))
        #P = Progress(size(X,1))
        for (i,row) in enumerate(eachrow(X))
            answer = findfirst(utils.squeeze(all(uX .== row[na, :], dims=2)))
            if answer != nothing
                g[i] = answer
            end
            #next!(P)
        end
        g
    end

    function ismember(X::AbstractVector, Y::AbstractVector)
        Z = ThreadsX.collect(any(Y .== x) for x in eachrow(X))
        replace(Z .!= nothing, missing=>false)
    end

    include("dict.jl")
    include("namedtup.jl")

end
