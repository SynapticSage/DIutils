module dict
    export remove_key_item
    export filter_keys
    function remove_key_item(k::Dict, item)
        if item ∈ keys(k)
            pop!(k, item)
        end
        k
    end

    """
    filters a dict by its keys
    """
    function filter_keys(d::Dict, lambda::Function)
        d = copy(d)
        lambda_keys!(d, lambda)
    end

    """
    filters a dict by its keys, with modification
    """
    function filter_keys!(d::Dict, lambda::Function)
        for key ∈ keys(d)
            v = pop!(d, key)
            key = lambda(key)
            d[key] = v
        end
        return d
    end
end
