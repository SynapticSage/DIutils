module namedtup
    #=
    This sub library is about functions for named
    tuples and for dicts who have named tuple
    keys.
    =#

    using Reexport
    export namedtupkeys_to_df, namedtuple_to_dict, remove_key_item
    export bestpartialmatch, countmatch
    export applyvalues
    export lambda_keys
    @reexport using NamedTupleTools

    function namedtuple_to_dict(X::NamedTuple)
        Dict(zip(keys(X), values(X)))
    end
    function pop(X::NamedTuple, key)
        X = namedtuple_to_dict(X)
        pop!(X, key)
        NamedTuple(X)
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
    function bestpartialmatch(K::Vector{<:NamedTuple},
                              search::NamedTuple)
        argmax(countmatch(K, search); kws...)
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
    function lambda_values(D::T where T <: AbstractDict{<:Any,<:AbstractDict},
            lambda::T where T <: Function) 
        Dict(key=>lambda_values(value, lambda) for (key,value) in D)
    end
    function lambda_values(D::T where T <: AbstractDict{<:Any,<:NamedTuple},
            lambda::T where T <: Function) 
        Dict(key=>lambda(value) for (key,value) in D)
    end
    #applyvalues = lambda_values


    """
    filters a dict by its keys, with modification
    """
    function lambda_keys(d::NamedTuple, lambda::Function; nested::Bool=false)
        d = Dict(zip(d))
        filter_keys(d, lambda; nested)
        d = NamedTuple(d)
        return d
    end

end
