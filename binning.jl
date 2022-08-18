module binning
    
    export inside, get_samptime
    using Statistics
    import ..Utils

    # Radius measurements
    function vector_dist(vals::Array, center::Array)::Vector{Float32}
        @inbounds @fastmath sqrt.(sum((vals .- center[Utils.na, :]).^2,dims=2)[:,1]) # HUGE SAVINGS FAST MATH HERE
    end
    function inside(vals::Array, center::Array, radius::Float32)::BitVector
        vector_dist(vals, center) .< radius
    end
    function get_samptime(vals::Array, center::Array, radius::Float32;
            dt::Float32=1/30)::Float32
        @fastmath sum(inside(vals, center, radius)) * dt
    end
    # Vector radius measurments
    function indiv_dist(vals::Array, center::Array)::Matrix{Float32}
        abs.(vals .- center[Utils.na, :])
    end
    function inside(vals::Array, center::Array, radius::Vector{Float32})::BitVector
        Utils.squeeze(all(indiv_dist(vals, center) .< radius[Utils.na, :], dims=2))
    end
    function get_samptime(vals::Array, center::Array, radius::Vector{Float32};
            dt::Float32=1/30)::Float32
        @fastmath sum(inside(vals, center, radius)) * dt
    end

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
    


end
