module group

    import ..Table
    import Utils
    using DataFrames
    using DataStructures
    using Infiltrator
    using StatsBase
    using DataFrames: ColumnIndex
    using Random
    using ProgressMeter

    coords(groups) = collect(zip(sort(collect(groups.keymap),by=x->x[2])...))[1]
    pairs(groups) = (collect(zip(sort(collect(groups.keymap))))[i][1] 
                     for i in 1:length(G.keymap))
    
    CItype = Union{ColumnIndex, Vector{<:ColumnIndex}}

    export annotate_periods!
    """
        annotate_periods!

    annotate dataframe with period data (see Table.get_period)
    """
    function annotate_periods!(data::DataFrame, periods::DataFrame;
                name=nothing,period_col=nothing,type=Int32)::DataFrame
        name = name === nothing ? periods[1,:prop] : name
        period_col = period_col === nothing ? name : period_col
        data[!,period_col] = Vector{Union{type,Missing}}(missing, size(data,1))
        period_start = allowmissing(Utils.searchsortedprevious.([periods.start], data.time))
        period_stop  = Utils.searchsortednext.([periods.stop],  data.time)
        period_start[period_start .!= period_stop] .= missing
        data[!,period_col] = period_start
        data
    end

    function named_coords(groups)
    end

    function find_common_mapping(groups::Union{Vector, Symbol, String},
            args::DataFrame...)
        find_common_mapping((groupby(arg, groups) for arg in args)...)
    end
    """
    finds a common index system to move through multiple dataframes
    """
    function find_common_mapping(args::GroupedDataFrame...)
        # -- Obtain a hash of the keys of each keymap --
        hashmaps = []
        for i ∈ 1:length(args)
            hashmap = OrderedDict()
            for (key, value) in args[i].keymap
                hashmap[hash(key)] = value
            end
            push!(hashmaps, hashmap)
        end
        # Intersect the hashmaps
        superset = [keys(hashmap) for hashmap in hashmaps]
        superset = Tuple(intersect(superset...))
        # Throw out keys that fail to match the superset
        for i ∈ 1:length(args)
            hashmap = OrderedDict(key=>hashmaps[i][key] for key in superset)
            hashmaps[i] = hashmap
        end
        # Now that all of the orderings match, I can read out the matrix of
        # group numbers per matching key
        [hashmaps[i][superset[j]] for j ∈ 1:length(superset), i ∈ 1:length(args)]
    end

    function mtg_via_commonmap(groups::Union{Vector, Symbol, String},
            args::DataFrame...)
        G = [groupby(arg, groups) for arg in args]
        M = find_common_mapping(G...)
        G = [[G[1][i],G[2][j]] for (i,j) in eachrow(M)]
    end

    """
    Combine two tables and groupby some props, while distinguishing the tables with a new key
    """
    function multitable_groupby(groups::Union{Vector, Symbol, String}, args::AbstractDataFrame...; dropmissingrows::Bool=false)
        #@info length(args)
        args = vcat(args...; source="source", cols=:union)
        G = [groupby(g, :source) for g in groupby(args, groups)]
        source_count, group_count = maximum([length(G[i]) for i in 1:length(G)]),
                                    length(G)
        #@info countmap([length(G[i]) for i in 1:length(G)])
        GG = Array{Union{Missing, SubDataFrame}}(undef, group_count, source_count)
        for i in 1:group_count
            for j in 1:length(G[i])
                GG[i,j] = G[i][j]
            end
            if length(G[i]) < source_count
                for j in (length(G[i])+1):source_count
                    GG[i,j] = missing
                end
            end
        end
        if dropmissingrows
            GG = GG[sum(ismissing.(G),dims=2).==0,:]
        end
        GG
    end

    function nt_keys(GroupDF::GroupedDataFrame)::Vector{NamedTuple}
        result = Vector{NamedTuple}()
        by = collect(values(GroupDF.keymap))
        sortme = collect(keys(GroupDF.keymap))
        K = sortme[by]
        for kk in K
            symbols = Tuple(Symbol(col) for col in GroupDF.cols)
            valuez = Tuple(val for val in kk)
            key_new = NamedTuple(
                                 Dict(a=>b for (a,b) in zip(symbols,valuez))
                                )
            push!(result, key_new)
        end
        return result
    end

    export equalize
    function equalize(X::DataFrame, group::CItype, balance::ColumnIndex; kws...)
        equalize(groupby(X,group), balance; kws...)
    end
    function equalize(X::GroupedDataFrame, balance::ColumnIndex; 
            relabel::Bool=true, thresh::Int = 1, qthresh = nothing)
        counts= combine(X, balance => (x->length(unique(x))) => :count)
        thresh = qthresh === nothing ? thresh : Int(round(quantile(counts.count, qthresh)))
        m = Int(max(minimum(counts.count), thresh))
        @info "thresh" m
        G=[]
        for g in X
            
            subg = groupby(g, balance)
            subg = filter((x->nrow(x) > m), subg)
            if length(subg) < m
                continue
            end

            # Balance our sets
            inds = collect(1:size(subg,1))
            shuffle!(inds)
            subg = subg[inds[1:m]]
            #Relabel
            if relabel
                subg = reindex(subg, balance)
            end
            push!(G,combine(subg,identity))
        end
        vcat(G...)
    end


    function reindex(X::GroupedDataFrame, dims)
        dims = dims isa Vector ? dims : [dims]
        uDimD = Dict()
        Xc = combine(X,identity)
        for dim in dims
            uDim = sort(unique(Xc[!,dim]))
            uDimD[dim] = Dict(uDim[i]=>i for i in 1:length(uDim))
        end
        for (x,dim) in Iterators.product(X,dims)
            x[!,dim] .= getindex.([uDimD[dim]], x[!,dim])
        end
        X
    end
    

end
