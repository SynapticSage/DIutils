module dict

    export remove_key_item
    export lambda_keys
    export to_dict
    using Infiltrator


    function flatten!(X::AbstractDict;keyfilt=x->true,valfilt=x->true)
        while any(typeof.(values(X)) .<: AbstractDict)
            for (K,V) in X
                if typeof(V) <: AbstractDict
                    v = Dict(k=>v for (k,v) in pop!(X, K)
                             if keyfilt(k) && valfilt(v))
                    try
                        merge!(X,v)
                    catch
                    end
                end
            end
        end
        X
    end
    flatten(X::AbstractDict; kws...) = flatten!(copy(X); kws...)

    function remove_key_item(k::Dict, item)
        if item ∈ keys(k)
            pop!(k, item)
        end
        k
    end

    """
    filters a dict by its keys
    """
    function lambda_keys(d::Dict, lambda::Function; nested::Bool=true)
        d = copy(d)
        lambda_keys!(d, lambda)
    end

    """
    filters a dict by its keys, with modification
    """
    function lambda_keys!(d::Dict, lambda::Function; nested::Bool=true)
        for key ∈ keys(d)
            v = pop!(d, key)
            key = lambda(key)
            if nested && typeof(v) <: Union{Dict,NamedTuple}
                d[key] = filter_keys(v, lambda; nested)
            else
                d[key] = v
            end
        end
        return d
    end

    function to_dict()
    end

    """
    filterchange_keys

    filter a dict by its keys and change those filter hits with a changee function
    """
    function filterchange_keys!(D::AbstractDict, filter::Function, change::Union{Function,Nothing}=identity)
        for k in keys(D)
            if filter(k)
                v = pop!(D, k)
                if change !== nothing
                    D[change(k)] = v
                end
            end
        end
        D
    end
    function filterchange_keys(D::AbstractDict, filter::Function, change::Union{Function,Nothing}=nothing)
        D = copy(D)
        filterchange_keys!(D, filter, change)
    end

    function load_dict_to_module!(mod::Module, dict::Dict)::Nothing
        load_keysvals_to_module!(mod, keys(dict), values(dict))
    end

    function load_keysvals_to_module!(mod::Module, keys, vals)::Nothing
        skipped = []
        for (k,v) in zip(keys, vals)
            println("Loading $k")
            try
            Core.eval(mod, :($(Symbol(k)) = $v))
            catch
                push!(skipped, k)
            end
        end
        if !isempty(skipped)
            @warn "skipped" skipped
        end
        nothing
    end

    function load_jld2_to_module!(fn::String)::Nothing
        
    end


    """
    Base.inv

    invert a dict keys=>values to values=>keys
    """
    Base.inv(dict::Dict) = Dict(value => key for (key, value) in dict)


end
