module filtreg

    import ..DIutils 


    findnearest = DIutils.searchsortednearest
    findnext    = DIutils.searchsortednext
    findprev    = DIutils.searchsortedprevious
    using Statistics, NaNStatistics
    using DataFrames
    import DataFrames: ColumnIndex
    using Infiltrator
    using ProgressMeter

    CItype = Union{ColumnIndex, Vector{<:ColumnIndex}}

    export downsample
    export register, registerEventsToContinuous, filterAndRegister
    export filter, filterTables

    """
        downsample(animal, day)

    Downsamples the data

    Input
    =====
    raster : DataFrame
    beh : DataFrame


    Output
    ======
    (raster, behavior)

    """
    function downsample(x; dfactor=10)
        δ = dfactor;
        downsamp = x[begin:δ:end, :];
        downsamp
    end

    """
        register

    function register(source::DataFrame, target::DataFrame; 
            transfer, on::String="time")::Vector{DataFrame}

    register columns in a `source` to a `target` dataframe `on` a certain

    ### Input
    `source` -- the source dataframe to register
    `target` -- the recipient dataframe
    `on` -- the column to register on
    `transfer` -- the columns to transfer from source to target
    `convert_type` -- the type to convert the columns to
    `tolerance` -- the tolerance for the registration
    `tolerance_violation` -- the value to use if the tolerance is violated
    `match` -- the type of match to use
    `checksort` -- whether to check if the data is sorted

    ### Output

    column
    """
    function register(source::AbstractDataFrame, target::AbstractDataFrame; 
            transfer::Vector{<:CItype},
            on::CItype="time",
            tolerance::Union{Float64, Nothing}=0.9999,
            tolerance_violation=missing,
            convert_type::Type=Float64,
            match=:nearest, checksort::Bool=true
        )::Tuple{DataFrame, DataFrame} 

        if checksort
            issorted(source[:, on]) || error("source not sorted")
            issorted(target[:, on]) || error("target not sorted")
        end

        match = if match == :nearest
            findnearest
        elseif match == :next
            findnext
        elseif match in [:prev, :previous]
            findprev
        end
        if tolerance === nothing
            @warn "No given tolerance"
        end
        # Get our columns into target
        @debug "→ → → → → → → → → → → → "
        @debug "Registration"
        @debug "→ → → → → → → → → → → → "

        if transfer == All()
            return source, target
        end
        @debug "columns=$transfer from source->target on $on"

        match_on_source = source[:, on]
        match_on_target = target[:, on]

        target_missings = ismissing.(match_on_target)
        if all(target_missings)
            @error "All target values are missing"
        elseif any(target_missings)
            ct = Union{convert_type, Missing}
        else
            ct = convert_type
        end
        match_on_target = typeof(match_on_target) == ct ? match_on_target :
                          convert(Vector{ct}, match_on_target)
        match_on_source = typeof(match_on_target) == ct ? match_on_source :
                          convert(Vector{ct}, match_on_source)

        # FindNearest
        match_on_source = (match_on_source,)
        indices_of_source_samples_in_target = match.(match_on_source,
                                                           match_on_target)

        # Tolerance
        if tolerance !== nothing
            δ = @inbounds source[indices_of_source_samples_in_target, on] -
                target[!, on]
            out_of_tolerance = abs.(δ) .> tolerance
        else
            out_of_tolerance = zeros(Bool, size(target,1))
        end
        nonmissing = (!).(ismissing.(out_of_tolerance))
        any_out_of_tol = any(out_of_tolerance[nonmissing])

        if tolerance_violation === nothing
            target = @inbounds target[Not(out_of_tolerance), :]
            @inbounds indices_of_source_samples_in_target = 
                indices_of_source_samples_in_target[Not(out_of_tolerance)]
        end

        for col ∈ transfer

            @debug any(out_of_tolerance) ? 
            "mean out of tolerance=>$(mean(out_of_tolerance))" :
            "no out of tolerance"

            rows = size(target,1)
            #if (col ∉ propertynames(target) && col ∉ names(target)) ||
            #    typeof(target[!,col]) != Vector{eltype(source[!,col])}
                target[!,col] = Vector{eltype(source[!,col])}(undef, rows)
            #end

            # Move columns
            target[!, col] = @inbounds source[indices_of_source_samples_in_target, col]

            # If tolerance vio is nothing, drop out of tolerance entries
            if tolerance_violation === missing
                target[!, col] = allowmissing(target[!,col])
            end

            # Write in tolerance violations
            if  (tolerance !== nothing &&
                tolerance_violation !== nothing && 
                any_out_of_tol)
                    @inbounds target[nonmissing .&& out_of_tolerance , col] .= tolerance_violation
            end
        end

        @debug "← ← ← ← ← ← ← ← ← ← ← ← "
        return source, target

    end


    """
        register(source::GroupedDataFrame, target::GroupedDataFrame, 
                      pos...; kwargs...)

    register columns in a `source` set of grouped dataframes to a `target` set
    of grouped dataframes
    # Arguments
    - `source` -- the source grouped dataframe to register
    - `target` -- the target grouped dataframe to register
    - `pos` -- arguments to pass to `register`
    - `kwargs` -- keyword arguments to pass to `register`
    """
    function register(source::GroupedDataFrame, target::GroupedDataFrame, 
                      pos...; kwargs...)
        source_keys = keys(source) |> collect .|> NamedTuple
        target_keys = keys(target) |> collect .|> NamedTuple
        K = intersect(source_keys, target_keys)
        for k ∈ K
            s = source[k]
            t = target[k]
            DIutils.filtreg.register(s, t, pos...; kwargs...)
        end
    end


    function registerEventsToContinuous(events::DataFrame, target::DataFrame;
            transfer::Union{Vector{String},String}, on::String="time",
            targetEltype::Dict{String,<:Type}=Dict{String,Any}(),
            targetCast::Dict{String,<:Function}=Dict{String,Function}(),
            ifNonMissingAppend::Bool=false,
            eventStart::String="start",
            eventStop::String="stop")::DataFrame 

        # Get our columns into target
        @debug "→ → → → → → → → → → → → "
        @debug "Registration"
        @debug "→ → → → → → → → → → → → "
        
        for col ∈ transfer
            if col ∉ names(target)
                @debug "$col not in target"
                T = col ∈ keys(targetEltype) ? targetEltype[col] : eltype(events[!,col])
                @debug "col=$col => type=$T"
                @debug "creating new cool with type = $T"
                target[!,col] = Array{Union{Missing,T}}(missing, size(target,1))
                # Cast to user specified type?
                if col in keys(targetCast)
                    target[!,col] = targetCast[col](target[!,col])
                end
            end
        end

        match = target[!, on]
        p = Progress(size(events,1), desc="Registering")
        Threads.@threads for event in eachrow(events)
            start = event[eventStart]
            stop  = event[eventStop]
            indices_of_source_samples_in_target = (match .>= start) .&&
                                                  (match .< stop)
            #@debug "sum=$(sum(indices_of_source_samples_in_target))"
            for col ∈ transfer
                if col == All()
                    continue
                end
                specialStringBehavior = ifNonMissingAppend && 
                         eltype(target[!, col]) == Union{Missing,String} # has to be String, not InlineString
                #@debug "$col is eltype=$(eltype(target[!,col])) and specialStringBehavior=$specialStringBehavior"
                if specialStringBehavior
                    #@debug "special behavior"
                    missingVals    = indices_of_source_samples_in_target .&& ismissing.(target[!,col])
                    nonMissingVals = indices_of_source_samples_in_target .&& (!).(ismissing.(target[!,col]))
                    existing       = target[nonMissingVals, col]
                    target[nonMissingVals, col] .= existing .* event[col]
                    target[missingVals, col]    .= event[col]
                else
                    target[indices_of_source_samples_in_target, col] .= event[col]
                end
            end
            next!(p)
        end
        @debug "← ← ← ← ← ← ← ← ← ← ← ← "
        return target
    end
    """
        filter(data::DataFrame...; filters::AbstractDict=Dict(),
                filter_skipmissingcols::Bool=false)

    instructions to query/filter values in a set of dataframes

    # Note
    this method works by copy, as opposed to the other methods which even
    without the ! symbol operate by reference. I didn't understand this !
    convention early on well when I started these libraries.
    """
    function filter(data::DataFrame...; filters::AbstractDict=Dict(),
            filter_skipmissingcols::Bool=false,
        )::Vector{DataFrame}
        data = [data...]
        @debug "→ → → → → → → → → → → → "
        @debug "Filtration"
        @debug "→ → → → → → → → → → → → "
        @inbounds for (cols, filt_for_cols) ∈ filters
            for i ∈ 1:length(data)

                @assert !(cols isa Bool) &&
                !(filt_for_cols isa Bool) &&
                !(filt_for_cols isa Vector{Bool}) &&
                !(cols isa Vector{Bool})

                list_cols = cols isa Vector ? String.(cols) : String.([cols])

                column_count = count(col->col∈names(data[i]), list_cols) 
                column_count_not_equal = column_count != length(list_cols)

                if filter_skipmissingcols && 
                    (any(c ∉ names(data[i]) for c ∈ list_cols) || column_count_not_equal)
                    @debug "data[$i] missing col ∈ cols=$cols, skip filter"
                    continue
                end

                if filt_for_cols isa Function
                    #println("Filter is a function")
                    @assert column_count == length(list_cols) "Uh oh, missing col, consider filter_skipmissingcols=true"
                    inds = @inbounds filt_for_cols(data[i][!, cols])
                elseif typeof(filt_for_cols) <: Vector
                    #println("Filter is a set of functions")
                    inds = @inbounds @fastmath accumulate(.&, [ff(data[i][!, cols]) for ff
                                           in filt_for_cols])[end]
                else
                    @debug "cols = $cols"
                    @debug "Typeof(cols) = $(typeof(cols))"
                    @debug "filt_for_cols = $filt_for_cols"
                    throw(TypeError(:filt_for_cols, "",Vector,typeof(filt_for_cols)))
                end

                @assert inds isa Vector || inds isa BitVector &&
                        !(inds isa Vector{Missing})

                if Missing <: eltype(inds)
                    replace!(inds, missing=>false)
                    inds = disallowmissing(inds)
                end
                @debug begin
                    percent = mean(inds)*100
                    "data_$i filtration: $percent percent pass filter on $cols"
                end
                try
                    @inbounds data[i] = data[i][inds, :];
                catch
                    @infiltrate
                end
            end
        end
        for i in 1:length(data)
            @debug "final_size(data($i)) = $(size(data[i]))"
        end
        @debug "← ← ← ← ← ← ← ← ← ← ← ← "
        return data
    end

    """
        filterAndRegister

    combination of  `raw.filter()` and `raw.register()` steps

    # Inputs
    filters
    source::DataFrame 
    target::DataFrame
    filters::Union{Nothing,AbstractDict}=nothing
    transfer=nothing,
    on="time"
    filter_skipmissingcols::Bool=false

    # Output
    source, target

    but these are also modified by reference
    """
    function filterAndRegister(source::DataFrame, target::DataFrame;
            filters::Union{Nothing,AbstractDict}=nothing, transfer=nothing,
            on="time", filter_skipmissingcols::Bool=false
            )::Vector{DataFrame}
        if filters !== nothing
            required_cols  = get_filter_req(filters)
            missing_fields = setdiff(required_cols, names(target))
            if !(isempty(missing_fields))
                @debug """Adding missing_fields=$missing_fields
                          to transfer=$transfer"""
                request = transfer === nothing ? [] : transfer
                transfer = string.(unique(collect(Iterators.flatten([request, missing_fields]))))
            end
        end
        @debug "register"
        if transfer !== nothing
            source, target = register(source, target ; transfer=transfer, on=on)
            #utils.piso(data)
        end
        @debug "filter"
        if filters !== nothing
        @info filter_skipmissingcols
            source, target = DIutils.filtreg.filter(source, target; filters=filters, filter_skipmissingcols)
        end
        return [source, target]
    end

    """
    `filterTables`

    alias for `filterAndRegister`
    """
    filterTables(data::DataFrame...; filters::Union{Nothing,AbstractDict}=nothing,
                 lookupcols=nothing, lookupon="time") =
    filterAndRegister(data...;transfer=lookupcols, on=lookupon, filters=filters)


    """
        get_filter_req

    Gets the fields required by the set of filters to operate
    """
    function get_filter_req(filts::AbstractDict)
        sets = [String.(f) for f ∈ keys(filts)]
        sets = [s isa Vector ? s : [s] for s in sets]
        unique(collect(Iterators.flatten(sets)))
    end


end
