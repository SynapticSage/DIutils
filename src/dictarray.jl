"""
    module DictArray

Paradigm for a dict indexed LIKE an array. Each key is a named tuple. The
values can be gotten like D[i, j, k] for key (i, j, k). Also, sets of
keys can be gotten D[:, j, k] for all i's with key (i, j, k). Importantly too,
the user could say D[k=2, j=3, i=Colon()]. In other words, they can name
dimensions

"""
module dictarray
    struct DictArray
    end
end
