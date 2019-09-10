"""
    sanity_check_volume(vol)

Throws an error if the array `vol` is not 3-dimensional, not cubic
or its first dimension is not even.
"""
function sanity_check_volume(vol)
    if ndims(vol) != 3
        error("[sanity_check_volume] volume dimensions incorrect (must be 3)!")
    elseif size(vol)[1] != size(vol)[2] || size(vol)[1] != size(vol)[3]
        error("[sanity_check_volume] volume has to be cubic!")
    elseif size(vol)[1] % 2 != 0
        error("[sanity_check_volume] volume edge length must be even!")
    end
end
