module DIutils

    using DrWatson, ProgressMeter, Infiltrator, Reexport

    import Random
    using CSV, DataFrames, Colors, ColorSchemes
    using Statistics, NaNStatistics, Plots, ThreadsX

    import SearchSortedNearest
    searchsortednearest  = SearchSortedNearest.searchsortednearest
    searchsortednext     = SearchSortedNearest.searchsortednext
    searchsortedprevious = SearchSortedNearest.searchsortedprevious

    export skipnan
    export itsizeof, piso
    export searchsortednext, searchsortednearest
    export remove_key_item
    export thropt

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

    #function norm_extrema(x::AbstractArray{T1}, minmax::Union{Vector{T2},Tuple{T2, T2}}) where
    #    T1 <: Real where T2 <: Real
    #    if minmax isa Tuple
    #        minmax = [minmax...]
    #    end
    #    if minmax[2] == minmax[1]
    #        x = minmax[1] * ones(size(x))
    #    else
    #        x = (x .- minimum(x))./(maximum(x) - minimum(x))
    #        x = x .* diff(minmax) .+ minmax[1]
    #    end
    #end

    export norm_extrema
    function norm_extrema(x::AbstractArray{T1}, 
             minmax::Union{Vector{T2},Tuple{T2, T2}}=[0,1] ) where
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

    export nannorm_extrema
    function nannorm_extrema(x::AbstractArray{T1}, minmax::Union{Vector{T2},Tuple{T2, T2}}=(0,1)) where
        T1 <: Real where T2 <: Real
        if minmax isa Tuple
            minmax = [minmax...]
        end
        if minmax[2] == minmax[1]
            x = minmax[1] * ones(size(x))
        else
            x = (x .- nanminimum(x))./(nanmaximum(x) - nanminimum(x))
            x = x .* diff(minmax) .+ minmax[1]
        end
    end

    function norm_percent(x::AbstractArray{T1}, quant::Real) where T1 <: Real 
        Q = any(isnan.(x)) ? nanquantile(x, quant) : quantile(x,quant)
        (x.-Q)./maximum(x) * 100
    end
    function norm_percent(x::AbstractArray{T1}) where T1 <: Real 
        x./maximum(x) * 100
    end

    function in_range(X::AbstractArray, range::Union{Tuple, Vector})
        X .≥ range[1] .&& X .< range[2]
    end
    function in_range(X::T where T<:Real, range::S where S<:Union{Tuple, Vector, SubArray{1}})
        X ≥ range[1] .&& X < range[2]
    end
    function not_in_range(X::Union{Real,AbstractArray}, range::Union{Tuple, Vector})
        (!).(in_range(X, range))
    end
    function in_rangeq(X::AbstractArray, qlim::Union{Tuple,Vector})
        range = [nanquantile(X,q) for q in qlim]
        in_range(X, range)
    end

    function vals_in_rangeq(X::AbstractArray, qlim::Union{Tuple,Vector})
        X[in_rangeq(X, qlim)]
    end

    export squeeze
    function squeeze(A::AbstractArray)  
        sz = size(A)
        A = dropdims(A, dims = tuple(findall(sz .== 1)...))
        return A
    end  
    function squeeze_zeros(A::AbstractArray, dims)
        inds = any(i->i != 0, A, dims=dims) 
        A[inds]
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

    # Move to Plots module
    # guides = Dict(); # shortcut for guides
    # guides[(:x, :stopWell)] = Guide.xlabel("Goal")
    # guides[(:y, :stopWell)] = Guide.ylabel("Goal")
    # guides[(:x, :neuron)] = Guide.xlabel("Neuron")
    # guides[(:y, :neuron)] = Guide.ylabel("Neuron")
    # guides[(:x, :rayleighZ)] = Guide.xlabel("Rayleigh Ζ")
    # guides[(:x, :rayleighZ_diff)] = Guide.xlabel("Rayleigh Z\nDifferences")
    # guides[(:x, :gt_shuffle)] = Guide.xlabel("Percent\nreal > shuff") 
    # """
    # TITLE: goalVectorTheme
    # Purpose: theme for goal vector shit
    # """
    # function goalVectorTheme()
    #     theme = Theme(major_label_color=colorant"white", major_label_font_size=14pt,
    #                   minor_label_color=colorant"white",
    #                   key_label_color=colorant"white",
    #                   key_title_color=colorant"white",
    #                   panel_fill=colorant"black",
    #                   background_color=colorant"black")
    #     Gadfly.push_theme(theme)
    # end

    # function getPushoverClient()
    #     token = open(expanduser("~/.pushover.token"),"r") do f
    #         token = read(f, String)
    #     end
    #     user = open(expanduser("~/.pushover.user"),"r") do f
    #         user = read(f, String)
    #     end
    #     return PushoverClient(user, token)
    # end

    export pushover
    function pushover(message; title=nothing)
        # send(getPushoverClient(), pos...; kws...)
        cmd = if title === nothing
            `pushover-cli "$message"`
        else
            `pushover-cli "$message" "$title"`
        end
        run(cmd)
    end


    function findgroups(pos...)
        X = Matrix(hcat(pos...))
        uX = unique(X, dims=1)
        g = zeros(Int,size(X,1))
        #P = Progress(size(X,1))
        for (i,row) in enumerate(eachrow(X))
            answer = findfirst(DIutils.squeeze(all(uX .== row[na, :], dims=2)))
            if answer !== nothing
                g[i] = answer
            end
            #next!(P)
        end
        g
    end

    """
    ∀ x ∈ X, x ∈ Y
    """
    function ismember(X::AbstractVector, Y::AbstractVector)
        Z = ThreadsX.collect(any(Y .== x) for x in eachrow(X))
        #replace(Z .!= nothing, missing=>false)
    end
    function ismember_trans(X::AbstractVector, Y::AbstractVector)
        any(X .== Y', dims=2)
    end

    function notisnan(X::Vector)
        (!).(isnan.(X))
    end

    function getind(val, pos, nd)
        ans = [Colon() for _ in 1:nd]
        ans[pos] = val
        ans
    end

    function indicesMatrixForm(A::AbstractArray)
        @time C = collect.(getproperty.(CartesianIndices(A),:I))
        Matrix(hcat(C...)')
    end
    

    include("dict.jl")
    include("namedtup.jl")
    include("macros.jl")
    include("binning.jl")
    include("filtreg.jl")
    include("arr.jl")
    include("statistic.jl")
    include("mlj.jl")
    include("plotutils.jl")
    include("clean.jl")
    include("Table.jl")
    # plotutils.plotutils
    export dict, namedtup, macros, binning, filtreg, arr, statistic, mlj, plotutils, clean, Table
    # @reexport using .dict, .namedtup, .macros, .binning, .filtreg, .arr, .statistic, .mlj, .plotutils, .clean, .Table
    #Reexport.reexport(plotutils)

end
