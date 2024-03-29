module dict

export remove_key_item
export lambda_keys
export to_dict
using Infiltrator
import ..DIutils
using DataFrames

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

function load_exfiltrated!(mod::Module, sess::Infiltrator.Session)
    k = propertynames(sess)        
    v = [getproperty(sess, kk) for kk in k]
    load_keysvals_to_module!(mod, k, v)
end

# TODO: UN commment after jld2 in package
# function save_exfiltrated(filename::String, sess::Infiltrator.Session)
#     if !endswith(filename, ".jld2")
#         filename = "$filename.jld2"
#     end
#     k = propertynames(sess)        
#     v = [getproperty(sess, kk) for kk in k]
#     JLD2.jldopen(filename, "w") do file
#         for (kk,vv) in zip(k,v)
#             print("Saving $kk")
#             file[String(kk)] = vv
#         end
#     end
# end
# #
# function load_jld2_to_module!(filename::String, mod::Module)::Nothing
#     if !endswith(filename, ".jld2")
#         filename = "$filename.jld2"
#     end
#     JLD2.open(filename, "r") do file
#         load_keysvals_to_module!(mod, keys(file), values(file))
#     end
# end
#
"""
Base.inv

invert a dict keys=>values to values=>keys
"""
Base.inv(dict::Dict) = Dict(value => key for (key, value) in dict)

export setkeybyfirst
"""
    setkeybyfirst(d::Dict)

set a dict key type by the type of the first key
"""
function setkeybyfirst(d::Dict)
    k = first(keys(d))
    Dict{typeof(k),valtype(d)}(d)
end

export setvalbyfirst
"""
    setvalbyfirst(d::Dict)

set a dict value type by the type of the first value
"""
function setvalbyfirst(d::Dict)
    v = first(values(d))
    Dict{keytype(d),typeof(v)}(d)
end

export setbyfirst
"""
    setbyfirst(d::Dict)

set a dict key and value type by the type of the first key and value
"""
function setbyfirst(d::Dict)
    k = first(keys(d))
    v = first(values(d))
    Dict{typeof(k),typeof(v)}(d)
end

function tostring(d::AbstractDict)
    DIutils.namedtup.tostring(NamedTuple(pairs(d)))
end

# ---- Dicts with namedtuple keys ----
"""
    getkeyprop(d::AbstractDict, prop::Symbol)

get a property of a key in a dict
"""
function getntprop(d::AbstractDict, prop::Symbol)
    map(keys(d)|>collect) do k
        getindex(k, prop)
    end 
end


"""
    getpropuniq(d::AbstractDict, prop::Symbol)

get a property of a key in a dict

# Arguments
-----------
d::AbstractDict
    dict to get property from
prop::Symbol
# Returns
--------
Array
    unique values of the property
"""
function getntpropuniq(d::AbstractDict, prop::Symbol)
    map(keys(d)|>collect) do k
        getindex(k, prop)
    end |> unique
end

function ntkeyframe(d::AbstractDict)::DataFrame
    keymat = hcat([(key|>collect) for key in keys(d)]...)
    keymat = permutedims(keymat, (2,1))
    keyframe = DataFrame(keymat,
        keys(d)|>first|>propertynames.|>string|>collect)
end

"""
getbinaryindex(d::AbstractDict, inds::AbstractVector{Bool})

get a dict with keys filtered by a binary index
"""
function getbinaryindex(d::AbstractDict, 
                        inds::AbstractVector{Bool})::AbstractDict
    typeof(d)(collect(keys(d))[inds], collect(values(d))[inds])
end

end
