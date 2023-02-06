
module mlj

    using Infiltrator
    using DataStructures: OrderedDict
    using DataFrames
    using MLJBase: round3, _standard_errors, PerformanceEvaluation

    function DataFrames.DataFrame(e::PerformanceEvaluation)
        _measure = [repr(MIME("text/plain"), m) for m in e.measure]
        _measurement = round3.(e.measurement)
        _per_fold = [round3.(v) for v in e.per_fold]
        _sterr = round3.(_standard_errors(e))

        # Only show the standard error if the number of folds is higher than 1.
        show_sterr = any(!isnothing, _sterr)
        data = show_sterr ?
            hcat(_measure, e.operation, _measurement, _sterr, _per_fold) :
            hcat(_measure, e.operation, _measurement, _per_fold)
        header = show_sterr ?
            ["measure", "operation", "measurement", "1.96*SE", "per_fold"] :
            ["measure", "operation", "measurement", "per_fold"]

        D = DataFrame(OrderedDict(h=>d for (h,d) in zip(header, eachcol(data))))
        D[!,:measure] = lowercase.([x[1] for x in split.(String.(Symbol.(D.measure)),'(')])
        D[!,:idxmeasure] = sortperm(D[!,:measure])
        hcat(D[!,[:measure]], D[!,Not(:measure)])
    end

    function measureidxDict(df::DataFrame)
        OrderedDict(unique(df.measure) .=> unique(df.idxmeasure))
    end

end
