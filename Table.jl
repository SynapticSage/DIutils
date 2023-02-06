"""
Module Contents:

* Time period functionality
** get_periods :: acquire start and stop ranges for distinct categorical
** contrain_otherperiods_by_period :: constrain a second range of periods by a first one
** contrain_range :: contrain data by a start and stop time
** binary on time :: return a set of time ranges where a property is a value
** select_group :: return a set of samples where a categorical property is  value
"""
module Table

    using DrWatson
    using DataFrames
    using DataFrames: ColumnIndex
    using ProgressMeter
    using Statistics
    using Blink
    using TableView
    using DataStructures
    using Infiltrator
    import Utils
    Utils = Utils
    using Reexport: @reexport
    export to_dataframe

    __revise_mode__ = :evalassign
    ∞ = Inf

    CItype = Union{ColumnIndex, Vector{<:ColumnIndex}, InvertedIndex}
    CItype_plusNull = Union{ColumnIndex, Vector{<:ColumnIndex}, Nothing}
    ColumnSelector = Union{Nothing,Vector{String},InvertedIndex,Cols,All,Between}
    DFColVars = Union{DataFrames.ColumnIndex, DataFrames.MultiColumnIndex}
                                                  
    # --.--o                             o         |
    #   |  .,-.-.,---.    ,---.,---.,---..,---.,---|
    #   |  || | ||---'    |   ||---'|    ||   ||   |
    #   `  `` ' '`---'    |---'`---'`    ``---'`---'
    #                     |                         
    #                                                          
    # ,---.               |    o               |    o|         
    # |__. .   .,---.,---.|--- .,---.,---.,---.|    .|--- ,   .
    # |    |   ||   ||    |    ||   ||   |,---||    ||    |   |
    # `    `---'`   '`---'`---'``---'`   '`---^`---'``---'`---|
    #                                                     `---'
    #
    function get_periods(df::DataFrame, property::CItype, pos...;
            removeMissing=false, end_period=:stop,
            timefract::Union{Function,Nothing,Pair}=nothing,
            kws...)
        if removeMissing
            df = dropmissing(df);
        end
        if timefract !== nothing
            timefract = if timefract isa Pair 
                sym, func = timefract
                sym => (x -> mean(func.(x))) => :frac
            else
                (x -> mean(timefract.(x))) => :frac
            end
            pos = [pos..., timefract]
        end
        period = groupby(df, property)
        period = combine(period, pos..., :time => (x->minimum(x)) => :start,
                                         :time => (x->maximum(x)) => end_period)
        period.δ = period[!,end_period] .- period.start;
        period.prop .= property isa Vector ? join(String.(property)," ") : String(property)
        return period
    end
    function isinperiod(x, start, stop)
        (x >= start) & (x < stop)
    end
    isin = isinperiod

    function remove_interperiod_time!(raster::DataFrame, period::DataFrame)
        prop = period.prop[1]
        dropmissing!(raster, prop)
        dropmissing!(period, prop)
        #naninds = isnan.(raster[!,prop])
        #if any(naninds)
        #    print(findall(naninds))
        #    delete!(raster, findall(naninds))
        #end
        raster_select(e) = raster[:, prop] .>= period[e, prop]
        period_select(e) = period[:, prop] .>= period[e, prop]
        for i in size(period,2):-1:2
            Δ = (period[i,"start"] - period[i-1,"end"])
            raster[raster_select(i),"time"] .-= Δ
            period[period_select(i),"start"] .-= Δ
            period[period_select(i),"end"] .-= Δ
        end
        return raster, prop
    end

    function remove_inactive_time!(inactive_time::Float64=1.0, data::DataFrame...)
        sorted_time = sort(data[1].time)
        diffs = diff(sorted_time)
        gap_starts = findall(diffs .> inactive_time)
        gap_stops = []
        for gap in gap_starts
            push!(gap_stops, sorted_time[gap+1])
        end
        gap_starts = sorted_time[gap_starts]
        for (gap_start, gap_stop) in collect(zip(gap_starts, gap_stops))[end:-1:begin]
            Δ = gap_stop - gap_start
            for d in 1:length(data)
                inds = data[d].time .> gap_stop
                data[d][inds,"time"] .-= Δ
                data[d][inds,"time"] .-= Δ
            end
        end
    end

    function constrain_otherperiods_by_period(period::DataFrame,
                                             τ_start, τ_end)
        α = findmax(period.end .≥ τ_start)[2];
        β = findmin(period.start .≤ τ_end)[2];
        constrained_period = period[ α:β ,:]
        constrained_period = trajs[begin:(end-1),:];
        τ_start = minimum(constrained_period.start);
        τ_end   = maximum(constrained_period.end);
        return (constrained_period, τ_start, τ_end)
    end

    function binary_on_times(beh::DataFrame, property::String, value;
                             τ_start=-∞, τ_stop=∞)
        # Determine memory periods
        
        δ = diff(beh[:,property] .== value);
        starts = beh.time[findall(δ == 1)];
        stops = beh.time[findall(δ == -1)];

        periods = ones(length(starts), 2) * NaN;
        candidates = [];
        for i = 1:length(starts)
            candidates = starts[i] .≤ stops;
            if any(candidates)
                besttime = argmax(candidates);
                periods[i, 1] = starts[i];
                periods[i, 2] = stops[besttime];
                starts[i] = ∞;
            end
        end
        constrain = (periods[!,1] .> τ_start) .& (periods[!,2] .≤ τ_stop);
        periods = periods[constrain,:];
    end

    function constrain_range(X::DataFrame, prop::String, τ_start, τ_stop)
        timeset = (X[!,prop] .≥ τ_start) .& (X[!,prop] .< τ_stop);
        X = X[timeset, :];
        return X
    end
    constrain_into_period = constrain_range


    """
    vec_of_matching_colnames

    grabs a vector of matching columns
    """
    function vec_of_matching_colnames(df, matchstr)
        if matchstr isa Symbol
            matchstr = String(matchstr)
        end
        if matchstr === nothing
            matching_cols = VectorString()
        else
            matching_cols = sort(filter(name -> occursin(matchstr, name),
                                        names(df)))
        end
    end

    function expand_colnames(df, colnames)
        newcolnames = Vector{String}([])
        for colname in colnames
            N = vec_of_matching_colnames(df, colname)
            for n in N
                print(n)
                push!(newcolnames, n)
            end
        end
        return newcolnames
    end

    """
    vec_of_matching_columns

    grabs a vector of matching columns
    """
    function vec_of_matching_columns(df, matchstr)
        matching_cols = sort(filter(name -> occursin(matchstr, name),
                                    names(df)))
        [x for x in eachcol(df[:, matching_cols])]
    end

    """
    select_row_group

    convenience function for selecting a categorical label
    """
    function select_row_group(raster, prop, value)
        raster[findall(skipmissing(raster[:,prop].==value)),:];
    end

    """
    `add_columns_from_othertable`

    adds from another table

    this method unfortunately requires there exist an index column
    """
    function add_columns_from_othertable(receiver::DataFrame, sender::DataFrame,
            columns::Vector{String}; index_column=nothing, time_match_column=nothing)
        #=
        index_column , optional
            which column in receiver used to index into send

            if time_match_column and index_column are Nothing, then this is assumed
            to be "index"

        time_match_column, optional
            the other method of adding is to use a time column to register the two
            tables
        =#
        if time_match_column ≠ nothing
            throw(ArgumentError("Not a supported feature yet"))
        end

        if index_column === nothing
            index_column = "index";
        end
        indices = receiver[:,index_column];
        valid = indices .> 0;
        column_expanded = []
        for col in columns
            if col in names(sender)
                push!(column_expanded, col)
            else
                push!(column_expanded, vec_of_matching_colnames(sender, col)...)
            end
        end

        for col in column_expanded
            elem = sender[1,col]
            if elem isa Union{Float16,Float32,Float64,ComplexF32,ComplexF64} && ~ismissing(elem)
                receiver[!,col] .= typeof(elem)(NaN)
            elseif sender[1,col] isa Int
                receiver[!,col] .= -1
            elseif sender[1,col] isa String
                receiver[!,col] .= missing
            else
                if any(ismissing.(sender[!,col]))
                    sender[!,col] = replace(sender[!,col], missing=>NaN)
                end
                receiver[!,col] .= NaN
            end
            try
                receiver[valid, col] = sender[indices[valid], col];
            catch
                println(typeof(receiver[!,col]), typeof(sender[!,col]))
                throw(ErrorException("Bad"))
            end
        end
        receiver
    end

    """
    `add_sort_properties`

    Adds cell centroid properties that we can sort plots by
    """
    function add_sort_properties(spikes, beh::DataFrame, 
            props::Vector{String}; modifier="", skipcenter=false)
        spike_cellgroups = groupby(spikes, "unit");
        @showprogress "Adding $modifier sort properties..." for prop in props
            if !(prop in names(beh))
                colnames = table.vec_of_matching_colnames(beh, prop);
                if isempty(colnames)
                    throw(ArgumentError("prop=$prop does not exist in behavior"))
                end
                prop_nums = String.(Symbol.(1:length(colnames)));
                prop_nums = "_" .* prop_nums;
            else
                prop_nums = [""];
            end
            for prop_num in prop_nums
                spike_cellgroups = table._add_sort_property(beh, spike_cellgroups, prop * prop_num;
                                                     modifier=modifier, skipcenter=skipcenter);
            end
        end
        spikes = combine(spike_cellgroups, x->x)
    end

    function _add_sort_property(beh, spike_cellgroups, prop; modifier="", skipcenter=false)
        if beh[!,prop] isa Vector{Complex}
            throw(ArgumentError("prop=$prop is complex"))
        end
        q = combine(spike_cellgroups, :index => (index->median(skipmissing(beh[index,prop]))) 
                    => :prop)
        if q.prop[1] isa Complex
            throw(ArgumentError("prop=$prop is complex"))
        end
        q[!,"prop_index"] = sortperm(q.prop);
        q[!,"inv_prop_index"] = invperm(q.prop_index);
        #sort(q,"prop_index")
        qfunc(cell)  = Dict(q.unit.=>q.prop)[cell]
        qifunc(cell) = Dict(q.unit.=>q.inv_prop_index)[cell]
        for c = 1:length(spike_cellgroups)
            if !skipcenter
                spike_cellgroups[c][!,modifier*prop*"_center"] .= qfunc.(spike_cellgroups[c].unit);
            end
            spike_cellgroups[c][!,modifier*prop*"_argcenter"] .= qifunc.(spike_cellgroups[c].unit);
        end
        return spike_cellgroups
    end

    include(srcdir("Table", "clean.jl"))
    @reexport using .clean
    include(srcdir("Table", "convert_types.jl"))
    @reexport using .convert_types
    include(srcdir("Table", "group.jl"))
    @reexport using .group
    include(srcdir("Table", "columntype.jl"))
    @reexport using .columntype

end

    #"""
    #`occupancy_normalize`

    ## Arguments
    #`data`
    #dataframe to normalize
    #`beh`
    #used to compute occupancy
    #`props`
    #properties to split by and normalize by counts within
    #`normalize_cols`
    #columns to normalize by 
    #"""
    #function occupancy_normalize(data::DataFrame, beh::DataFrame, 
    #        props::Union{Vector{String},String}; 
    #        normalize_cols::ColumnSelector=nothing)
    #    data = _occupancy_normalize(data, beh, props, normalize_cols=normalize_cols)
    #end
    #function occupancy_normalize!(data::DataFrame, beh::DataFrame, 
    #        props::Union{Vector{String},String}; 
    #        normalize_cols::ColumnSelector=nothing)
    #    _occupancy_normalize(data, beh, props, normalize_cols=normalize_cols, inplace=true);
    #    return nothing
    #end
    #function _occupancy_normalize(data::DataFrame, beh::DataFrame, 
    #        props::Union{Vector{String}, String}; 
    #        normalize_cols::ColumnSelector=nothing, inplace::Bool=false)
    #    if !inplace
    #        data = copy(data)
    #    end
    #    if props isa String
    #        props = [props]
    #    end
    #    if normalize_cols === nothing
    #        normalize_cols = Not(props)
    #    end
    #    uProps = [unique(beh[!,prop]) for prop in props]
    #    for state in Iterators.product(uProps...)
    #        bools = [beh[!,prop] .== state[i] for (i, prop) in enumerate(props)]
    #        bools = accumulate(.&, bools)
    #        bools = bools[end]
    #        count = accumulate(.+, bools)[end]
    #        bools = [data[!,prop] .== state[i] for (i, prop) in enumerate(props)]
    #        bools = accumulate(.&, bools)
    #        bools = bools[end]
    #        try
    #            data[bools, normalize_cols] ./= count
    #        catch InexactError
    #            data[!, normalize_cols] = convert.(Float64, data[!, normalize_cols]);
    #            data[bools, normalize_cols] ./= count
    #        end
    #    end
    #    return data
    #end

