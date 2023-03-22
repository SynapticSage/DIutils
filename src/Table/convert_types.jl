"""
    convert_types

This module contains functions for creating dataframes from various
types of data

# `to_` methods

- to_dataframe(field::AbstractArray; other_labels=Dict())
- to_dataframe(fields::AbstractDict; other_labels=Dict(), 
            key_name::Union{Nothing,Symbol,String,Vector}=nothing,
            exit_hiccup::Bool=false,
            level::Int=0, level_names::Dict=Dict(), kws...)
- to_dataframe(fields::T where T<:NamedTuple; kws...)
- to_dataframe(X::DataFrame; kws...)
- to_dataframe(F::Union{AbstractArray, Real}; props::Vector{String},
                grid::Tuple, gridtype::String, other_labels=Dict(),
                name::String="value", explode::Bool=false, kws...)

Of these dispatches, the most important is AbstractDict. There is a tendency
for the leaves of these calls to hit one of the other types of dispatches,
in particular the AbstractArray, DataFrame, and Real dispatches.

# `from_` methods

rather than being a part of the dispatch herarchy of to_dataframe, these
methods are one off implementations to get a dataframe from some type
with special conditions

- from_dimarrary(X::DimArray, registrant::Tuple{DataFrame, String,
    Vector{String}}...)::DataFrame

"""
module convert_types

    using DataFrames
    export to_dataframe
    using LazyGrids: ndgrid
    using DataStructures
    using Infiltrator
    using DimensionalData
    import DIutils 
    using DIutils.binning: edge_to_center

    export from_dimarrary
    """
        rate_todataframe

    converts a rate dimarray (or count dimarray) to a dataframe, and populates
    columns labeling the times from other dataframes

    # Arguments
    - `X`: a rate dimarray
    - `registrate`: a tuple of (source, on, transfer) where
        - `source` is a dataframe
        - `on` is the column name in `source` to match to `X`
        - `transfer` is a vector of column names in `source` to transfer to `X`

    # Returns
    - `X`: a dataframe with the same dimensions as `X` and the columns from `registrant`
    """
    function from_dimarrary(X::DimArray, registrant::Tuple{DataFrame, String,
    Vector{String}}...)::DataFrame
        dims = X.dims
        time = collect(dims[1])
        X = DataFrame([time, Vector.(collect(eachrow(X)))], [:time, :data])
        for (source, on, transfer) in registrant
            DIutils.filtreg.register(source, X; on, transfer)
        end
        X
    end

    export report_recursion
    """
        report_recursion(dict::AbstractDict)

    # Purpose
    reports the levels that can be recursed to in a dictionary,
    showing the first key and value at each level, as well as their types
    and then calls itself on the value if it is a dictionary
    """
    function report_recursion(dict::AbstractDict; printval=true, 
            separatorcount=2)
        key = first((keys(dict)))
        println("key: $key, typeof(key): $(typeof(key))")
        println("")
        printval ? println("value: $(dict[key])") : nothing
        println("typeof(value): $(typeof(dict[key]))")
        for _ in 1:separatorcount; println("-"^100); end
        if typeof(dict[key]) <: AbstractDict
            report_recursion(dict[key])
        end
    end

    export to_dataframe
    """
    to_dataframe

    # Purpose
    converts a field dictionary (multiple units) into a field dataframe

    # Arguments
    - `fields`: a dictionary of fields, where each field is a dictionary of units
    - `other_labels`: a dictionary of labels to add to the dataframe
    - `key_name`: the name of the key to use for the dataframe
            if nothing, will use the key name from the field dictionary 
            if a vector, recursively peels off key_name for the vector
            if a symbol or string, will use that as the key name
    - `exit_hiccup`: if true, will exit the function if the field is a hiccup
    - `level`: the level of recursion :: used to keep track of the recursion
    - `level_names`: a dictionary of level names :: used to keep track of the
                    recursive level names
    - `_level_hiccup_tracker`: a dictionary to track the first instance where
                        a level is not named ("hiccup") so to show the user
                        the only the first instance
    - `kws`: keyword arguments
        see to_dataframe(field::AbstractArray; other_labels=Dict())
    """
    function to_dataframe(fields::AbstractDict; other_labels=Dict(), 
            key_name::Union{Nothing,Symbol,String,Vector}=nothing,
            exit_hiccup::Bool=false,
            level::Int=0, level_names::Dict=Dict(),
            _level_hiccup_tracker::Dict=Dict(),
            kws...)::DataFrame

        # increment the level tracker for recursion
        level += 1
        
        #@info "function start" key_name level level_names

        D = DataFrame()
        for key in keys(fields)

            #@info "nested key:" level level_names fields key key_name
            kn = "unnamed"

            if (key isa NamedTuple) || (key isa Dict)
                key_dict = Dict(string(k)=>v for (k, v) in pairs(key))
                other_labels = merge(other_labels, key_dict)
                kn = nothing

            elseif key_name !== nothing #&& (key isa NamedTuple || key isa Dict)

                if (typeof(key_name) <: Vector)
                    if level ∉ keys(level_names)
                        if !(isempty(key_name))
                            kn = popfirst!(key_name)
                            level_names[level] = kn
                        end
                    else
                        kn = level_names[level]
                    end
                elseif (key_name == :keyboard)
                    if level ∉ keys(level_names)
                        println("Name key for key=$key, typeof(key)=$key")
                        kn = readline()
                        level_names[level] = kn
                    else
                        kn = level_names[level]
                    end
                else
                    if !(level-1 ∈ keys(level_names) && 
                         level_names[level-1] == "__key_level")
                        level_names[level] = "__key_level"
                        kn = key_name
                    end
                end
            end

            if kn == "unnamed"
                @debug "unhandled key_name" level key key_name other_labels
            end

            if kn !== nothing
                other_labels[kn] = key
            end
            
            if fields[key] !== nothing
                try
                    if key_name !== nothing || (key_name isa Vector &&
                                                !(isempty(key_name)))
                        kws = (;kws..., key_name)
                    end
                    df = to_dataframe(fields[key]; other_labels, level,
                                      level_names, _level_hiccup_tracker,
                                      kws...)
                    append!(D, df , cols=:union)
                catch exception
                    if !haskey(_level_hiccup_tracker, level)
                        @warn "hiccup" level key key_name other_labels
                        _level_hiccup_tracker[level] = 
                            (;key, fields, key_name, level, exception)
                    end
                    if exit_hiccup
                        return nothing
                    end
                end
            end
        end
        return D
    end

    """
        to_dataframe(field::T where T <: NamedTuple)
    
    if we hit a named tuple, then convert it to a dictionary and dispatch
    """
    function to_dataframe(fields::T where T <: NamedTuple; kws...)::DataFrame
        D = to_dataframe(DIutils.namedtuple_to_dict(fields); kws...) 
        return D
    end

    """
        to_dataframe(X::DataFrame; other_labels=nothing, kws...)

    if parent levels of to_dataframe encounter a dataframe while dispatching,
    then return a dataframe with the other_labels added to it
    """
    function to_dataframe(X::DataFrame; other_labels=nothing, kws...)::DataFrame
        # if other_labels !== nothing, then add them to the dataframe
        if other_labels !== nothing
            for (k,v) in other_labels
                if typeof(v) <: Array
                    v = string(v)
                end
                X[!,k] .= v
            end
        end
        X
    end

    """
        to_dataframe(F::Union{AbstractArray,Real};
    props::Vector{String}=Vector{String}([]), grid::Tuple=(), gridtype="center",
    other_labels=Dict(), name::String="value", explode::Bool=true, kws...)

    converts a single field MATRIX into a field DATAFRAME

    # Arguments
    - `F`: a field matrix
    - `props`: a vector of properties to add to the dataframe
    - `grid`: a tuple of grid vectors, namely the aliases for the dimension
                    indices
    - `gridtype`: the type of grid, either "center" or "edge"
    - `other_labels`: a dictionary of labels to add to the dataframe
    - `name`: the name of the value column
    - `explode`: if true, will explode multi-dimension fields into 
                dataframe rows and label dimensions as a column
    
    # Returns
    - `D`: a dataframe of the field
    """
    function to_dataframe(F::Union{AbstractArray,Real};
            props::Vector{String}=Vector{String}([]), grid::Tuple=(),
            gridtype="center", other_labels=Dict(), name::String="value",
            explode::Bool=true, kws...)::DataFrame

        # exploding vector fields into rows requires us to describe the
        # indices of their dimensions
        if explode
            D = ndgrid((1:size(F,i) for i in 1:ndims(F))...)
        else
            D = [1:size(F,i) for i in 1:ndims(F)]
        end
        # let's create a dictinoary that will expand these indices into
        # a dataframe, where keys are the dimension names 
        D = OrderedDict{String,Any}("dim_$d"=>vec(D[d])
                              for d in 1:ndims(F))
        # if the field is a real number, then we need to wrap it in a vector
        if typeof(F) <: Real
            F = [F]
        end
        # if props are given, then let us potentially hand the aliases of
        # the dimensions to the F array already given
        if ~isempty(props)
            if grid[1] isa Tuple
                grid = Tuple([g...] for g in grid)
            end
            if gridtype == "edge"            
                grid = edge_to_center.(grid) 
            end
            if explode
                grid = ndgrid(grid...)
            end
        end
        # let's clean the other labels
        _clean_label_values(other_labels)
        # let's add the other labels to the dictionary which will become
        # the dataframe
        for (label, value) in other_labels
            if label isa Symbol
                label = String(label)
            end
            D[label] = value
        end
        # if we are exploding the field, then let's add the field to the
        # dictionary and convert the dictionary to a dataframe
        if explode
            for (prop, G) in zip(props, grid)
                D[prop] = vec(G)
            end
            D[name] = vec(F);
        # if we are not exploding the field, then let's add the field to
        # the dictionary wrapped in a vector and convert the dictionary
        else
            for (key, value) in D
                D[key] = [value]
            end
            D[name] = [F,];
        end
        D = DataFrame(D)
    end

    """
        _clean_label_values(L::AbstractDict)

    cleans the label values in a dictionary by converting them to strings and
    joining them if they are vectors of strings

    # Arguments
    - `L`: a dictionary of labels
    """
    function _clean_label_values(L::AbstractDict)
        for (k,v) in L
            typev = typeof(v)
            if typev <: AbstractVector
                if eltype(v) == String
                    v = join(v, "-")
                else
                    v = string(v)
                end
                L[k] = v
            end
        end
    end

end
