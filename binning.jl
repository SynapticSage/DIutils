module binning
    
    using Statistics
    import ..Utils
    using Entropies
    import DataStructures: OrderedDict
    using DataFrames
    using LazyGrids: ndgrid
    using ProgressMeter
    using Missings
    using LazySets
    using Infiltrator

    function center_to_edge(grid::AbstractVector)
        grid = collect(grid)
        Δ = median(diff(grid))
        δ = Δ/2
        grid = collect(minimum(grid)-δ:Δ:maximum(grid)+δ)
    end
    function edge_to_center(grid::AbstractArray)
        grid = collect(grid)
        grid = dropdims(mean([vec(grid[1:end-1]) vec(grid[2:end])], dims=2), dims=2)
    end

    digitize(X::AbstractArray, nbins) = 
        Int16.(floor.(Utils.norm_extrema(X, [0, nbins-1])) .+ 1)

    export Grid
    abstract type Grid end
    export Occupancy
    abstract type Occupancy end

    #=
      ____      _     _ 
     / ___|_ __(_) __| |
    | |  _| '__| |/ _` |
    | |_| | |  | | (_| |
     \____|_|  |_|\__,_|
    =#

    export GridAdaptive
    struct GridAdaptive <: Grid

        props::Array{String}
        centers::Tuple
        edges::Tuple
        grid::Array{Array{Float32}}
        radii::Union{Array{Float32},
               Array{Vector{Float32}}}

        function GridAdaptive(props::Vector, centers::Union{Array,Tuple}; 
                radiidefault=:default)
            if eltype(props) == Symbol
                props = String.(props)
            end
            if radiidefault == :single_halfwidth
                radii = get_default_radii(centers)
            elseif radiidefault == :vector_halfwidth
                radii = get_default_vector_radii(centers)
            else
                radii = radiidefault
            end
            centers = centers isa Array ? Tuple(centers) : centers
            GridAdaptive(props, centers, radii)
        end
        function GridAdaptive(props::Vector, centers::Union{<:AbstractArray,Tuple},
            radii::Union{Float32,Vector{Float32}})
            if eltype(props) == Symbol
                props = String.(props)
            end
            centers = centers isa Array ? Tuple(centers) : centers
            centers = Tuple(Float32.(c) for c in centers)
            grid = sparse_to_full(centers)
            radii = fill(radii, size(grid))
            edges = center_to_edge.([[c...] for c in centers])
            edges = Tuple((e...,) for e in edges)
            new(props, centers, edges, grid, radii)
        end
        function GridAdaptive(props::Vector, centers::Tuple, grid::Array, radii::Array)
            if eltype(props) == Symbol
                props = String.(props)
            end
            @assert(size(grid) == size(radii))
            edges = center_to_edge.([[c...] for c in centers])
            edges = Tuple((e...,) for e in edges)
            new(props, centers, edges, grid, radii)
        end
        function GridAdaptive(props::Vector; width::Vector, boundary::Vector,
                radiidefault::Union{Symbol, Float32, Vector{Float32}}=:default)
            centers = Tuple(Tuple(collect(s:w:e))
                            for (w, (s, e)) in zip(width, boundary))
            GridAdaptive(props, centers; radiidefault)
        end
    end


    # Setup iteration
    Base.length(g::GridAdaptive)  = length(g.grid)
    Base.size(g::GridAdaptive)    = size(g.grid)
    Base.iterate(g::GridAdaptive) = Base.iterate(zip(g.grid, g.radii))
    #Base.done(g::GridAdaptive, state::Int) = length(g.centers) == state
    function Base.iterate(g::GridAdaptive, state::Tuple{Int,Int})
        iterate(zip(g.grid, g.radii), state)
    end
    cenumerate(g::GridAdaptive) = zip(CartesianIndices(g.grid), g)
    export content
    """
        nanradii(g::GridAdaptive)

    gets the location of naned radii
    """
    function content(g::GridAdaptive; nantemplate=false)  
        res = map(x->any(isnan.(x)), g.radii)
        if nantemplate
            replace(res, 1=>NaN, 0=>1)
        else
            (!).(res)
        end
    end

    Base.ndims(grd::GridAdaptive) = length(grd.props)

    function Base.sum(grd::GridAdaptive; dims)
        dims = Tuple(collect(dims))
        newdims = setdiff(1:ndims(grd), dims)
        grid, props, centers = getindex.(grd.grid,[newdims]), 
                               grd.props[newdims],
                               grd.centers[newdims]
        mean_sq(a,b) = sqrt.( a.^2 .+ b.^2 )
        radii   = accumulate(mean_sq, grd.radii, dims=dims)
        GridAdaptive(props, centers, grd.edges, grid, radii)
    end

    export get_grid_bounded
    """
        get_grid_bounded(behavior::DataFrame, props::Vector;
            widths::OrderedDict, boundary::OrderedDict,
            thresh::Float32=1.25f0, # Threshold in seconds
            adaptive::Bool=true,
            dt::Union{Nothing,Float32}=nothing, # Total time of sample
            radiusinc::Union{Float32,Vector{Float32}}=0.2f0, # Spatial unit of RF
            ϵ::Float32=0.1f0,
            maxrad::Union{Vector{Float32},Float32,Nothing}=nothing,
            radiidefault::Union{Symbol, Float32, Vector{Float32}}=:default,
            method::Symbol=:converge_to_radius,
            info::Bool=false,
            kws...)::GridAdaptive

    obtains the dynamic sampling grid from only the animals behavior

    # Parameters
    radiidefault -- (Optional) 
        :default|:single_halfwidth|:vector_halfwidth|Float32|Vector{Float32}

    ### Notes
    if radiusinc is vector, then it will instead expand a hyperbox instead 
    of a hypersphere
    """
    function get_grid_bounded(behavior::DataFrame, props::Vector;
        widths::OrderedDict, boundary::OrderedDict,
        thresh::Float32=1.25f0, # Threshold in seconds
        adaptive::Bool=true,
        dt::Union{Nothing,Float32}=nothing, # Total time of sample
        radiusinc::Union{Float32,Vector{Float32}}=0.2f0, # Spatial unit of RF
        ϵ::Float32=0.1f0,
        maxrad::Union{Vector{Float32},Float32,Nothing}=nothing,
        method::Symbol=:converge_to_radius,
        radiidefault::Union{Symbol, Float32, Vector{Float32}}=:default,
        info::Bool=false,
        steplimit=nothing,
        kws...)::GridAdaptive

        #ensureTimescale!(behavior)
        dt = dt === nothing ?
             median(diff(behavior.time)) : dt
        widths = OrderedDict(prop => widths[prop] for prop in props)
        maxrad = maxrad === nothing ? 3 * maximum(values(widths)) : maxrad
        vals = return_vals(behavior, props)
        cv(x) = collect(values(x))
        G = GridAdaptive(props; width=cv(widths), boundary=cv(boundary), 
                        radiidefault)
        radiusinc = ((valtype(G.radii) <: Vector) && !(typeof(radiusinc) <: Vector)) ?
                    fill(radiusinc, length(G.centers)) : radiusinc
        R = Vector{valtype(G.radii)}(undef, length(G))
        if method == :converge_to_radius_w_inertia
            thresh += ϵ
        end
        method = adaptive ? eval(method) : fixed_radius
        radiusinc = radiusinc .+ 1
        if info || !(isdefined(Main, :PlutoRunner))
            @info "grid" thresh dt maxrad radiusinc 
            @info boundary 
            @info widths
            prog = P = Progress(length(G); desc="Grid")
        end
        @debug "Beginning threaded get_grid_bounded loop"
        #Threads.@threads for (index, (this_center, radius)) in collect(enumerate(G))
        Threads.@threads for (index, (this_center, radius)) in collect(enumerate(G))
            #i+=1
            @debug "pre $(i)" radius G.radii
            dt = Float32.(dt)
            newradius = method(vals, this_center, radius; dt, thresh, maxrad,
                radiusinc, ϵ, steplimit)
            @debug "post $(i)" newradius G.radii
            R[index] = newradius
            if !(isdefined(Main, :PlutoRunner))
                next!(P)
            end
        end
        G.radii .= reshape(R, size(G))
        return G
    end

    export get_grid
    """
            get_grid(behavior::DataFrame, props::Vector;
        widths::Union{<:Float32,Vector{<:Float32},OrderedDict},
        boundary=nothing, other_kws...)::GridAdaptive

    get an adaptive grid object
    """
    function get_grid(df::DataFrame, props::Vector;
        widths::Union{<:Float32,Vector{<:Float32},OrderedDict},
        boundary=nothing, other_kws...)::GridAdaptive
        if typeof(widths) <: Float32
            widths = OrderedDict{}(prop => widths for prop in props)
        elseif typeof(widths) <: AbstractVector
            widths = OrderedDict{}(prop => width for (width, prop) in zip(widths, props))
        else
            @assert(widths isa OrderedDict)
        end
        bd = get_boundary(df, props)
        boundary = boundary === nothing ? bd :
                   fill_missing_boundary(boundary, bd)
        get_grid_bounded(df, props; widths=widths, boundary, other_kws...)
    end
    """
        get_grid(grid::GridAdaptive, task::DataFrames)

    Cuts out points that are out of bounds from the task grid
    """
    function get_grid(grid::GridAdaptive, tsk::DataFrame; scale=1)
        boundarypoints = :epoch in propertynames(tsk) ? 
            first(groupby(tsk[tsk.name .== "boundary", :],:epoch)) : 
                 tsk[tsk.name .== "boundary", :]
        xb, yb = boundarypoints.x, boundarypoints.y
        append!(xb, xb[1]); append!(yb, yb[1]);
        taskbound = VPolygon(hcat(xy,yb)'; apply_convex_hull=false)
        if scale != 1
            # Get a linear scale
            q = taskbound .* scale
            # Get the new boundary
            V′ = VPolygon(q.M * hcat(q.X.vertices...))
            # Shift it to the old
            vdiff = mean(V′.vertices) .- mean(taskbound.vertices)
            # Get the new polygon
            taskbound = VPolygon(V′.vertices .- [vdiff])
        end

        # NOw let's identify the grid inside the task boundary
        filt = Singleton.(grid.centers) .∈ taskbound
        centers =  copy(grid.centers)
        for c ∈ centers[filt]
            c .= NaN
        end

        GridAdaptive(grid.props, centers, )

    end

    function fill_missing_boundary(bd_userinput::AbstractDict,
                                   bd_template::OrderedDict)::OrderedDict
        for (key, value) in bd_userinput
            if value !== nothing
                bd_template[key] = bd_userinput[key]
            end
        end
        return bd_template
    end

    """
        converge_to_radius

    simple converge to radius. it expands a circle/sphere until a certain
    amount of sample time is acheived. it stops when threshold is reached.
    if it reaches maxrad radius before that threshold, this radius becomes
    NaN, nullifying this sample.
    """
    function converge_to_radius(vals::Matrix{Float32}, center::Vector{Float32},
            radius::Float32; dt::Float32, thresh::Float32, 
            steplimit::Union{Nothing,Int}=nothing,
            maxrad::Float32, radiusinc::Float32, kws...)::Float32
        #list = []
        #push!(list, radius)
        steplimit = steplimit === nothing ? Inf : steplimit
        @inbounds while get_samptime(vals, center, radius; dt) < thresh
            @fastmath radius = radius * radiusinc
            #push!(list, radius)
            if radius > maxrad
                radius = NaN32
                break
            end
            if (steplimit -= 1) <= 0
                break
            end
        end
        radius
    end
    function converge_to_radius(vals::Matrix{Float32}, center::Vector{Float32},
        radius::Vector{Float32}; dt::Float32, thresh::Float32,
        maxrad::Union{Float32,Vector{Float32}},
        steplimit=nothing,
        radiusinc::Union{Float32,Vector{Float32}}, kws...)::Vector{Float32}
        #list = []
        #push!(list, radius)
        #radius = copy(radius)
        steplimit = steplimit === nothing ? Inf : steplimit
        @inbounds while get_samptime(vals, center, radius; dt) < thresh
            @debug radius
            @fastmath radius = radius .* radiusinc
            #push!(list, radius)
            if any(radius .> maxrad) || any(radius .=== NaN32)
                radius .= NaN32
                break
            end
            if (steplimit -= 1) <= 0
                break
            end
        end
        radius
    end

    function fixed_radius(vals::Matrix{Float32}, center::Vector{Float32},
            radius::Union{Float32,Vector{Float32}}; 
            dt::Float32, kws...)::Vector{Float32}

        samptime = get_samptime(vals, center, radius; dt)
        if samptime == 0
            radius = fill(NaN32, size(radius))
        end
        radius
    end
    function fixed_radius(vals::Matrix{Float32}, center::Vector{Float32},
            radius::Float32; dt::Float32, kws...)::Float32
        samptime = get_samptime(vals, center, radius; dt)
        if samptime == 0
            radius = NaN32
        end
        radius
    end


    function converge_to_radius_w_inertia(vals::Matrix{Float32},
        center::Vector{Float32}, radius::Float32; dt::Float32,
        thresh::Float32, maxrad::Float32, radiusinc::Float32,
        ϵ::Float32, kws...)::Float32
        Δ = []
        tₛ = get_samptime(vals, center, radius; dt)
        δ = tₛ - thresh
        push!(Δ, δ)
        s = sign(δ)
        tolerance_acheived = abs(δ) < ϵ
        if s == 1
            tolerance_acheived = true
        end
        increase_phase = true
        reversal = 1
        while !(tolerance_acheived)
            while s != reversal && !(tolerance_acheived)
                radius += radiusinc
                if increase_phase
                    radiusinc *= 2
                else
                    radiusinc /= 2
                end
                tₛ = get_samptime(vals, center, radius; dt)
                δ = tₛ - thresh
                push!(Δ, δ)
                s = sign(δ)
                tolerance_acheived = abs(δ) < ϵ
            end
            reversal *= -1
            radiusinc *= -1
            increase_phase = false
            if abs(radiusinc) < 0.01
                break
            end
        end
        @debug "tolernace acheived = $tolerance_acheived"
        if radius > maxrad
            radius = NaN
        end
        radius
    end


    #=
    ,---.    |          |    o           
    |---|,---|,---.,---.|--- ..    ,,---.
    |   ||   |,---||   ||    | \  / |---'
    `   '`---'`---^|---'`---'`  `'  `---'
                   |                     
    =#

    """
        default_radii

    Determines how much to expand a hyperphere to fit the corners of the grid
    tagent to the sphere
    """
    function get_default_radii(centers::Tuple)::Union{Float32,Vector{Float32}}
        C = Vector{Float32}()
        # Push the max diff of each dimension
        for center in centers
            c = maximum(diff([center...]))
            push!(C, c)
        end
        # Turn diameter into radius
        C ./= 2
        # Return the euclidean distance to a corner of the hypercube from the
        # center
        if length(unique(C)) == 1
            sqrt(sum(C .^ 2))
        else
            C
        end
    end

    """
        get_default_vector_radii

    Determines how much to expand a hyperphere to fit the corners of the grid
    tagent to the sphere
    """
    function get_default_vector_radii(centers::Tuple)::Union{Float32,Vector{Float32}}
        C = Vector{Float32}()
        # Push the max diff of each dimension
        C = vcat(centers[:]...)
        [median(diff(sort(unique(c)))./2) for c in C]
    end


    #=
      ___                                              
     / _ \  ___ ___ _   _ _ __   __ _ _ __   ___ _   _ 
    | | | |/ __/ __| | | | '_ \ / _` | '_ \ / __| | | |
    | |_| | (_| (__| |_| | |_) | (_| | | | | (__| |_| |
     \___/ \___\___|\__,_| .__/ \__,_|_| |_|\___|\__, |
                         |_|                     |___/ 
    =#
    export IndexedAdaptiveOcc
    struct IndexedAdaptiveOcc <: Occupancy
        grid::GridAdaptive
        count::Array{Int32}
        prob::Probabilities
        camerarate::Float32
        inds::Array{Vector{Int32}}
    end
    Base.length(g::IndexedAdaptiveOcc) = length(g.grid)
    Base.size(g::IndexedAdaptiveOcc)   = size(g.grid)
    Base.iterate(g::IndexedAdaptiveOcc) = Base.iterate(V) do V
        zip(g.grid.grid, inds)
    end
    export get_occupancy_indexed
    function get_occupancy_indexed(data::DataFrame, 
                           grid::GridAdaptive)::IndexedAdaptiveOcc
        vals = return_vals(data, grid.props)
        count = zeros(Int32, size(grid))
        inds = Array{Union{Missing,Vector{Int}}}(missing, size(grid)...)
        for (index, (center, radius)) in collect(enumerate(grid))
            binary_locs = inside(vals, center, radius)
            @inbounds inds[index] = findall(binary_locs)
            @inbounds count[index] = sum(binary_locs)
        end
        count = reshape(count, size(grid))
        inds  = reshape(inds, size(grid))
        prob = Probabilities(Float32.(vec(count)))
        if :time in propertynames(data)
            camerarate = Float32(1 / median(diff(data.time)))
            if camerarate > 300
                camerarate /= 60
            end
        else
            camerarate = 1
        end
        IndexedAdaptiveOcc(grid, count, prob, camerarate, inds)
    end

    export AdaptiveOcc
    struct AdaptiveOcc <: Occupancy
        grid::GridAdaptive
        count::Array{Int32}
        prob::Probabilities
        camerarate::Float32
    end
    function get_occupancy(data::DataFrame, 
                           grid::GridAdaptive)::AdaptiveOcc
        vals = return_vals(data, grid.props)
        count = zeros(Int32, size(grid))
        Threads.@threads for (index, (center, radius)) in collect(enumerate(grid))
            @inbounds count[index] = sum(inside(vals, center, radius))
        end
        count = reshape(count, size(grid))
        prob = Probabilities(Float32.(vec(count)))
        if :time in propertynames(data)
            camerarate = Float32(1 / median(diff(data.time)))
            if camerarate > 300
                camerarate /= 60
            end
        else
            camerarate = 1
        end
        AdaptiveOcc(grid, count, prob, camerarate)
    end

    # ----------------
    # Helper functions
    # ----------------
    # Radius measurements
    export vector_dist
    function vector_dist(vals::Array, center::Array)::Vector{Float32}
        @inbounds @fastmath sqrt.(sum((vals .- center[Utils.na, :]) .^ 2, dims=2)[:, 1]) # HUGE SAVINGS FAST MATH HERE
    end
    function indiv_dist(vals::Array, center::Array)::Matrix{Float32}
        abs.(vals .- center[Utils.na, :])
    end
    # Vector radius measurments
    export inside
    function inside(vals::Array, center::Array, radius::Float32)::BitVector
        vector_dist(vals, center) .< radius
    end
    function inside(vals::Array, center::Array, radius::Vector{Float32})::BitVector
        Utils.squeeze(all(indiv_dist(vals, center) .< radius[Utils.na, :], dims=2))
    end
    export get_samptime
    function get_samptime(vals::Array, center::Array, radius::Float32;
        dt::Float32=1 / 30)::Float32
        @fastmath sum(inside(vals, center, radius)) * dt
    end
    function get_samptime(vals::Array, center::Array, radius::Vector{Float32};
        dt::Float32=1 / 30)::Float32
        @fastmath sum(inside(vals, center, radius)) * dt
    end


    """
        get_boundary

    gets the boundary of each property
    """
    function get_boundary(df::DataFrame, props::Vector)::OrderedDict
        boundary   = OrderedDict{eltype(props),Any}()
        for prop in props
            boundary[prop]   = extrema(Utils.skipnan(collect(skipmissing(df[!,prop]))))
        end
        boundary
    end

    """
        resolution_to_width

    converts resolution to width, ie the number of points spread over a
    range of points into the inter-sample width
    """
    function resolution_to_width(resolution::OrderedDict,
            boundary::OrderedDict)::OrderedDict
        width = OrderedDict{String,Any}()
        for prop in keys(resolution)
            width[prop] = boundary[prop] / resolution[prop]
        end
        width
    end

    """
        return_vals

    return a value matrix/dataframe for the requested
    properties in the DataFrame X
    """
    function return_vals(X::DataFrame, props::Vector)::Union{Vector{Float32},
                                                             Matrix{Float32}}
        Y = dropmissing(X[!, props])
        Float32.(hcat([Y[!,prop] for prop in props]...))
    end

    sparse_to_full(G::Vector...) = [g for g in zip(ndgrid(G...)...)]
    sparse_to_full(G::Tuple{Vector}) = [g for g in zip(ndgrid(G...)...)]
    function sparse_to_full(sparse::Tuple...)::Array{Tuple}
        C = [[c...] for c in sparse]
        C = [c for c in zip(ndgrid(C...)...)]
        [collect(c) for c in C]
    end
    function sparse_to_full(sparse::Tuple)::Array{Array}
        C = [[c...] for c in sparse]
        C = [c for c in zip(ndgrid(C...)...)]
        [collect(x) for x in C]
    end
    function full_to_sparse(full::Array)::Array{Array}
        accessor = Vector{Union{Colon,<:Int}}(undef,ndims(full))
        accessor .= 1
        sparse = []
        for i in 1:ndims(full)
            access = copy(accessor)
            access[i] = Colon()
            g = full[access...]
            g = Tuple(g[i] for g in g)
            push!(sparse,g)
        end
        sparse
    end
end
