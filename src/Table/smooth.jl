
module  smooth

    """
     gauss

    applies a gaussian kernel to a field of a dataframe 
    """
    function gauss(df::DataFrame; fields=["phase"], gaussian=3)
        kernel = Kernel.gaussian((gaussian,))
        for field in fields
            df[!,field] = imfilter(df[!,field], kernel)
        end
        return df
    end

end
