using Test
using DIutils
using DataFrames

@testset begin "Adaptive Grid and Occpancy"

    d = DataFrame(rand(1000,3), :auto)
    
    epsilon = 1f-6
    widths  = [0.02f0, 0.02f0, 0.02f0]
    grid_kws =
            (;
            widths,
	    radiusinc    = [0.2f0, 0.2f0, 0f0],
            maxrad       = [0.2f0, 0.2f0, 1f0],
            radiidefault = widths./2 .- epsilon,
            steplimit=3,
	    dt=0.1f0,
        ) 

    g = DIutils.binning.get_grid(d, [:x1, :x2, :x3]; grid_kws...)
    o = DIutils.binning.get_occupancy_indexed(d, g)


    function test_ind_consistent_count(o)
        sumocccount = sum(o.count)
        sumdatainds = sum(length.(collect(skipmissing(o.datainds))))
        @test sumocccount == sumdatainds 
    end

    @test begin "Adaptive Indexed Occupancy indices == sum of the count"
        test_ind_consistent_count(o)
    end

end
