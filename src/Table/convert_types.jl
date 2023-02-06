module convert_types

    using DataFrames
    export to_dataframe
    using LazyGrids: ndgrid
    using DataStructures
    using Infiltrator
    using DimensionalData
    import Utils

    export from_dimarrary
    """
        rate_todataframe

    converts a rate dimarray (or count dimarray) to a dataframe, and populates columns labeling the
    times from other dataframes
    """
    function from_dimarrary(X::DimArray, registrant::Tuple{DataFrame, String, Vector{String}}...)::DataFrame
        dims = X.dims
        time = collect(dims[1])
        X = DataFrame([time, Vector.(collect(eachrow(X)))], [:time, :data])
        for (source, on, transfer) in registrant
            Utils.filtreg.register(source, X; on, transfer)
        end
        X
    end

    export to_dataframe
    """
    to_dataframe

    #purpose
    converts a field dictionary (multiple units) into a field dataframe
    """
    function to_dataframe(fields::AbstractDict; other_labels=Dict(), 
            key_name::Union{Nothing,Symbol,String,Vector}=nothing,
            exit_hiccup::Bool=false,
            level::Int=0, level_names::Dict=Dict(), kws...)::DataFrame

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
                @warn "unhandled key_name" level key key_name other_labels
                @infiltrate
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
                                      level_names, kws...)
                    append!(D, df , cols=:union)
                catch
                    @warn "hiccup" level key key_name other_labels
                    @infiltrate
                    if exit_hiccup
                        return nothing
                    end
                end
            end
        end
        return D
    end
    function to_dataframe(fields::T where T <: NamedTuple; kws...)::DataFrame
        D = to_dataframe(Utils.namedtuple_to_dict(fields); kws...) 
        return D
    end
    function to_dataframe(X::DataFrame; other_labels=nothing, kws...)::DataFrame
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
    field_to_dataframe(field::AbstractArray; other_labels=Dict())

    +purpose: converts a single field matrix into a field dataframe
    """
    function to_dataframe(F::Union{AbstractArray,Real};
            props::Vector{String}=Vector{String}([]), grid::Tuple=(),
            gridtype="center", other_labels=Dict(), name::String="value",
            explode::Bool=true, kws...)::DataFrame
        if explode
            D = ndgrid((1:size(F,i) for i in 1:ndims(F))...)
        else
            D = [1:size(F,i) for i in 1:ndims(F)]
        end
        D = OrderedDict{String,Any}("dim_$d"=>vec(D[d])
                              for d in 1:ndims(F))
        if typeof(F) <: Real
            F = [F]
        end
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
        _clean_label_values(other_labels)
        for (label, value) in other_labels
            if label isa Symbol
                label = String(label)
            end
            D[label] = value
        end
        if explode
            for (prop, G) in zip(props, grid)
                D[prop] = vec(G)
            end
            D[name] = vec(F);
        else
            for (key, value) in D
                D[key] = [value]
            end
            D[name] = [F,];
        end
        D = DataFrame(D)
    end

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
