module plot

    import Utils
    import Plots
    using Dates


    function set_theme_timebased()
        if Utils.in_range(hour(now()), [0,5]) ||
           Utils.in_range(hour(now()), [20, 24])
            Plots.theme(:dark)
            theme="dark"
        else
            Plots.theme(:bright)
            theme="bright"
        end
    end
end
