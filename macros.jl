module macros

    using ProgressMeter, ProgressLogging
    export @thropt, @prog
    export thropt, prog

    macro identity(n)
        return n
    end

    macro thropt(control, expr)
        !(control) && return expr
        esc(quote
            Threads.@threads $expr
            end)
    end

    #macro progress(control, expr)
    #    control && return expr
    #    esc(quote
    #        @ $expr
    #        end)
    #end
    #
    #
    macro prog(pkg, expr)
        if isdefined(Main, :PlutoRunner) 
           esc(quote
                @progress $expr
            end)
        else
           esc(quote
                @progressmeter $expr
            end)
       end
    end

end

