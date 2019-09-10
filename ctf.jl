"""
    ctf_2d_arr(x, y, size, defocus_u, defocus_v, cs, apix, kv, ac, b_fac)

Calculate the ctf value for the given array parameters at the given
point (`x`, `y`).
"""
function ctf_2d_arr(x,  # x coordinates
                    y,  # y coordinates
                    size,  # unit: pixel
                    defocus_u,  # unit: Angstrom
                    defocus_v,  # unit: Angstrom
                    defocus_a,  # in degree
                    cs,  # unit: mm
                    apix=1.0,  # unit: Angstrom/pixel
                    kv=200,  # unit: kilo voltage
                    ac=0.1,  # arbitrary unit
                    b_fac=0.0)
    # physical length
    length = size * apix

    # convert to correct unit
    azimuthal_a = defocus_a * pi / 180  # convert to rad
    local_v = kv * 1e3  # unit volt
    local_cs = cs * 1e7  # unit A

    # wavelength
    wavelength = 12.2643247 / sqrt(local_v * (1.0 + local_v * 0.978466e-6))
    wavelength_3 = wavelength^3

    # reciprocal k vector
    k_angle = atan(y, x)
    k = sqrt((x / length)^2 + (y / length)^2)
    k_2 = k^2
    k_4 = k_2^2
    a_shift = k_angle - azimuthal_a

    # defocus
    defocus = defocus_u * cos(a_shift) * cos(a_shift) +
              defocus_v * sin(a_shift) * sin(a_shift)

    # phase shift
    ph = pi / 2 * local_cs * wavelength_3 * k_4 -
         pi * wavelength * defocus * k_2

    # ctf
    ctf = -sqrt(1.0 - ac * ac) * sin(ph) + ac * cos(ph)

    # damp ctf based on b_fac
    if b_fac > 0
        ctf *= exp(-b_fac / 4.0 * k_2)
    end

    return ctf
end
