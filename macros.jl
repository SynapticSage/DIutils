module macros

    macro identity(n)
        return n
    end

    macro prog(expr)
        if isdefined(Main, :PlutoRunner)
           :(@progres $(esc(expr)))
       else
           :(@identity $(esc(expr)))
       end
    end

end

