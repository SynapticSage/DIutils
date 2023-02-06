module columntype

    using DataFrames
    export handle_complex_columns
    export vec_arrayofarrays!

    """
        handle_complex_columns

    Converts complex columns into two real columns, optionally dropping
    the original columns
    """
    function handle_complex_columns(df::DataFrame; amp_label="amp",
            angle_label="ang", drop=true, type=Float32, props_to_mod=nothing)
        df_orig = df
        df = dropmissing(df)
        nametypes = zip(names(df),eltype.(eachcol(df)))
        iscomplex(type) = ((type == ComplexF32) || (type == ComplexF64) || (type == Complex))
        nametypes = filter(x->iscomplex(x[2]), collect(nametypes))
        complex_types = [name for (name,type) in nametypes]
        df = nothing
        for name in copy(complex_types)
            df_orig[:, name] = replace(df_orig[!, name], missing=>NaN)
            amp_name = amp_label*titlecase(name,strict=false)
            df_orig[:, amp_name]   = type.(abs.(df_orig[!, name]))
            phase_name = angle_label*titlecase(name,strict=false)
            df_orig[:, phase_name] = type.(angle.(df_orig[!, name]))
            println(props_to_mod, name)
            if props_to_mod !== nothing && any(startswith.(name, props_to_mod))
                push!(props_to_mod, amp_name)
                push!(props_to_mod, phase_name)
            end
            if drop
                df_orig = df_orig[!, Not(name)]
                if props_to_mod !== nothing
                    props_to_mod = filter(item->item≠name, props_to_mod)
                end
            end
        end
        if props_to_mod !== nothing
            return df_orig, props_to_mod
        else
            return df_orig
        end
    end

    function vec_arrayofarrays!(df::DataFrame)
        for (i,col) in enumerate(eachcol(df))
            if Missings.nonmissingtype(eltype(col)) <: Vector &&
                length(col[1]) == 1
                df[!,i] = vcat(col...)
            end
        end
    end

    
end
