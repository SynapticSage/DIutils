module dict

    export remove_key_item
    export lambda_keys
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


end
