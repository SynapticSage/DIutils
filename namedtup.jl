module namedtup
    #=
    This sub library is about functions for named
    tuples and for dicts who have named tuple
    keys.
    =#

    using Reexport
    using DataFrames
    using DataStructures: OrderedDict

    export namedtupkeys_to_df, namedtuple_to_dict, remove_key_item
    export bestpartialmatch, argbestpartialmatch, countmatch
    export applyvalues, lambda_key, removeprops, reorderprops
    @reexport using NamedTupleTools

    removeprops(keyset::Base.KeySet, key::Vector{Symbol}) =
                                            removeprops(collect(keyset), key)
    function removeprops(keyset::Vector{NamedTuple}, key::Vector{Symbol})
        for i in 1:length(keyset)
            keyset[i] = pop(keyset[i], key)
        end
        keyset
    end
    function reorderprops(keyset::Vector{NamedTuple}, 
                          pref_order::Vector{Symbol})
        for i in 1:length(keyset)
            keyset[i] = reorderprops(keyset[i], pref_order)
        end
        keyset
    end
    function reorderprops(X::NamedTuple, keys::Vector{Symbol})
        pn     = propertynames(X)
        keys   = intersect(keys, pn)
        unspec = setdiff(pn, keys)
        out = merge(OrderedDict(k=>getproperty(X, k) for k in keys),
                    OrderedDict(u=>getproperty(X, u) for u in unspec))
        NamedTuple(out)
    end
    function pop(X::NamedTuple, key::Symbol)
        X = namedtuple_to_dict(X)
        key ∈ keys(X) ? pop!(X, key) : nothing
        NamedTuple(X)
    end
    function pop(X::NamedTuple, keyz::Vector{Symbol})
        X = namedtuple_to_dict(X)
        for key in keyz
            key ∈ keys(X) ? pop!(X, key) : nothing
        end
        NamedTuple(X)
    end

    function namedtuple_to_dict(X::NamedTuple)
        Dict(zip(keys(X), values(X)))
    end
    function namedtupkeys_to_df(K::Base.KeySet)
        df = DataFrame()
        for k in K
            k = Dict(zip(keys(k), values(k)))
            k = DataFrame(k)
            append!(df, k)
        end
        df
    end

    function remove_key_item(k::NamedTuple, item)
        k = Dict(zip(keys(k), values(k)))
        k = remove_key_item(k, item)
        NamedTuple(k)
    end

    # Dicts with named tuple keys
    # May eventually split this off to its own sublib

    function namedtupkeys_to_df(D::AbstractDict)
        K = keys(D)
        namedtupkeys_to_df(K)
    end

    #                                                                             
    #.    ,          |               ,-,   .                   |--.--          -. 
    #|    |,---.,---.|--- ,---.,---. | |\  |,---.,-.-.,---.,---|  |  .   .,---. | 
    # \  / |---'|    |    |   ||    -: | \ |,---|| | ||---'|   |  |  |   ||   | :-
    #  `'  `---'`---'`---'`---'`     | `  `'`---^` ' '`---'`---'  `  `---'|---' | 
    #                                `-                                   |    -' 
    """
    Counts the number of matches in vector of namedTuples to a search namedtuple
    """
    function countmatch(N::Vector{<:NamedTuple}, search::NamedTuple;
        tolerance_kws...)
        results = []
        for n in N
            key_results = []
            keys_intersect = intersect(keys(n), keys(search))
            for k in keys_intersect
                if typeof(n[k]) <: Real
                    match = isapprox(n[k], search[k]; tolerance_kws...)
                else
                    match = n[k] == search[k]
                end
                push!(key_results, match)
            end
            push!(results,sum(key_results))
        end
        return results
    end
    """
    finds the key that is the best partial match to a search term
    """
    function argbestpartialmatch(K::Base.KeySet, search::NamedTuple)
        bestpartialmatch(collect(K), search)
    end
    function argbestpartialmatch(K::Vector{<:NamedTuple}, search::NamedTuple)
        argmax(countmatch(K, search))
    end
    function bestpartialmatch(K::Base.KeySet, search::NamedTuple)
        bestpartialmatch(collect(K), search)
    end
    function bestpartialmatch(K::Vector{<:NamedTuple}, search::NamedTuple)
        K[argmax(countmatch(K, search))]
    end

    #                                                                
    # |   /                        ,---.                 |o     |    
    # |__/ ,---.,   .,---.    ,---.|__.     ,---.    ,---|.,---.|--- 
    # |  \ |---'|   |`---.    |   ||        ,---|    |   |||    |    
    # `   ``---'`---|`---'    `---'`        `---^    `---'``---'`---'
    #           `---'                                                
    """
    Counts the number of matches a dict with named tuple keys to a search
    namedtuple
    """
    function countmatch(D::AbstractDict{<:NamedTuple,Any}, 
                                 search::NamedTuple; kws...)
        countmatch(keys(D), search; kws...)
    end

    """
    finds the key of a dict{namedtuple,any} that is the best partial match to a
    search term
    """
    function bestpartialmatch(D::AbstractDict{<:NamedTuple, Any},
                              search::NamedTuple)
        argmax(countmatch(D, search); kws...)
    end
    function bestentry(D::AbstractDict{<:NamedTuple, Any}, search::NamedTuple)
        best_key = argmax(countmatch(D, search); kws...)
        D[collect(keys(D))[best_key]]
    end
    function new_func()
    end

    #                                                                           
    # .    ,     |                            ,---.                 |o     |    
    # |    |,---.|    .   .,---.,---.    ,---.|__.     ,---.    ,---|.,---.|--- 
    #  \  / ,---||    |   ||---'`---.    |   ||        ,---|    |   |||    |    
    #   `'  `---^`---'`---'`---'`---'    `---'`        `---^    `---'``---'`---'
    function lambda_values(D::T where T <: AbstractDict{<:NamedTuple,<:Any},
            lambda::T where T <: Function) 
        Dict(key=>lambda_values(value, lambda) for (key,value) in D)
    end
    function lambda_values(D::T where T <: AbstractDict{<:Any,<:AbstractDict},
            lambda::T where T <: Function) 
        Dict(key=>lambda_values(value, lambda) for (key,value) in D)
    end
    function lambda_values(D::T where T <: AbstractDict{<:Any,<:NamedTuple},
            lambda::T where T <: Function) 
        Dict(key=>lambda(value) for (key,value) in D)
    end
    applyvalues = lambda_values


    """
    return a namedtuple changed
    """
    function lambda_key(d::NamedTuple, lambda::Function; nested::Bool=false)
        d = Dict(zip(d))
        filter_keys(d, lambda; nested)
        d = NamedTuple(d)
        return d
    end

end
