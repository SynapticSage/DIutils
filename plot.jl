module plot

    import Utils
    import Plots
    using Dates


    function set_theme_timebased(time::Union{Real,Nothing}=nothing)
        time = time === nothing ? hour(now()) : time
        if Utils.in_range(time, [0,5]) ||
           Utils.in_range(time, [20, 24])
            Plots.theme(:dark)
            theme="dark"
        else
            Plots.theme(:bright)
            theme="bright"
        end
    end
    set_theme_timebased() = set_theme_timebased(nothing)
end
