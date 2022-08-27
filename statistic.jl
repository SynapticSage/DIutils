module statistic

    export  pfunc

    function pfunc(x::Real)
        if x < 1e-3
            "***"
        elseif x < 1e-2
            "**"
        elseif x < 0.05
            "*"
        elseif x < 0.1
            "†"
        else
            ""  
        end
    end
    pfunc(x::Vector{<:Real}) = pfunc.(x)

end
