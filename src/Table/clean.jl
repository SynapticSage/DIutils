module clean

    export clean_duplicate_cols
    export nan_to_missing, nan_to_missing!
    export naninf_to_missing, naninf_to_missing!
    using DataFrames
    #function combine(data::Dict{Any, Any})
    #    print("Dict combine")
    #    return vcat((x[2] for x in data)...)
    #end

    function nan_to_missing(df::DataFrame, cols::Union{Vector{String}, Vector{Symbol}})
        for col in cols
            df[:, col] = replace(df[!, col], NaN=>missing)
        end
        return df
    end
    function naninf_to_missing(df::DataFrame, cols::Union{Vector{String}, Vector{Symbol}})
        for col in cols
            df[:, col] = replace(df[!, col], NaN=>missing, Inf=>missing,
                                 -Inf=>missing)
        end
        return df
    end

    function nan_to_missing!(df::DataFrame, cols::Union{Vector{String}, Vector{Symbol}})
        for col in cols
            df[!, col] = replace(df[!, col], NaN=>missing)
        end
    end
    function naninf_to_missing!(df::DataFrame, cols::Union{Vector{String}, Vector{Symbol}})
        for col in cols
            df[!, col] = replace(df[!, col], NaN=>missing, Inf=>missing,
                                 -Inf=>missing)
        end
    end

    function clean_duplicate_cols(df::DataFrame)
        for detect in "_" .* string.(1:4)
            splits = split.(names(df), detect)
            check_these = splits[length.(splits) .> 1]
            for check in check_these
                col1, col2 = join(check,detect), check[1]
                dat1, dat2 = eachcol(dropmissing(df[!,[col1,col2]]))
                if all(dat1 .== dat2)
                    @info "Removing $col1"
                    df = df[!, Not(col1)]
                end
            end
        end
        return df
    end
end
