include("./ctf.jl")
include("./rotation.jl")
include("./check.jl")

function sinc_pad!(vol_pad, vol, trilinear)
    sz = size(vol)[1]
    padded_sz = size(vol_pad)[1]
    pad_start = div(padded_sz - sz, 2)
    @inbounds for i in 0:sz-1
        l = i - div(sz, 2)
        for j in 0:sz-1
            k = j - div(sz, 2)
            for t in 0:sz-1
                h = t - div(sz, 2)

                # sinc correct
                r = sqrt(h^2 + k^2 + l^2) / padded_sz
                if r != 0
                    sinc = sin(r * pi) / (r * pi)
                else
                    sinc = 1.0
                end
                if trilinear
                    value = vol[i+1, j+1, t+1] / sinc^2
                else
                    value = vol[i+1, j+1, t+1] / sinc
                end

                # padding
                vol_pad[i+pad_start+1, j+pad_start+1, t+pad_start+1] = value
            end
        end
    end
end

function normalize_volume(vol)
    mean = sum(vol) / length(vol)
    stddev = sqrt(sum((vol .- mean) .^ 2)) / (length(vol) - 1)
    return (vol .- mean) ./ stddev
end

function preprocess_for_projection(vol, pad=2, trilinear=true)
    sanity_check_volume(vol)
    vol = normalize_volume(vol)
    padded_size = size(vol)[1] * pad
    vol_pad = zeros(padded_size, padded_size, padded_size)
    sinc_pad!(vol_pad, vol, trilinear)

    f_3d = rfft(vol_pad, [3, 2, 1])
    return f_3d
end

function projection!(f_2d, f_3d,
                     phi, theta, psi, x_shift, y_shift,
                     pad=2,
                     ctf=false,
                     defocus_u=-1, defocus_v=-1, defocus_a=-1,
                     cs=-1,
                     apix=-1,
                     kv=-1,
                     ac=-1,
                     b_fac=-1,
                     white=false,
                     trilinear=true)
    size_2d = size(f_2d)[1]
    half_size_2d = div(size_2d, 2)
    size_3d = size(f_3d)[1]
    max_r = div(size_3d, 2) - 1
    @inbounds for j in 0:size(f_2d)[2] - 1 # Julia arrays are stored in column major order, so looping over f_2d like this is faster
        h = j

        for i in 0:size_2d - 1
            if i >= half_size_2d
                k = i - size_2d
            else
                k = i
            end
#    for i in 0:size_2d - 1
#        if i >= half_size_2d
#            k = i - size_2d
#        else
#            k = i
#        end
#        for j in 0:size(f_2d)[2] - 1
#            h = j

            rot_mat = zeros(3, 3)
            rotation_matrix_euler!(rot_mat, phi, theta, psi, false)
            x = (rot_mat[1, 1] * h + rot_mat[1, 2] * k) * pad
            y = (rot_mat[2, 1] * h + rot_mat[2, 2] * k) * pad
            z = (rot_mat[3, 1] * h + rot_mat[3, 2] * k) * pad

            if x^2 + y^2 + z^2 < max_r^2
                if x < 0
                    neg = true
                    x = -x
                    y = -y
                    z = -z
                else
                    neg = false
                end
                if trilinear
                    x0 = floor(Int, x)
                    x1 = x0 + 1
                    xd = x - x0
                    y0 = floor(Int, y)
                    y1 = y0 + 1
                    yd = y - y0
                    z0 = floor(Int, z)
                    z1 = z0 + 1
                    zd = z - z0
                    if (y0 < 0) y0 += size_3d end
                    if (y1 < 0) y1 += size_3d end
                    if (z0 < 0) z0 += size_3d end
                    if (z1 < 0) z1 += size_3d end

                    if (z0 + y0 + x0) % 2 == 1
                        flip = -1
                    else
                        flip = 1
                    end

                    d000 = f_3d[z0+1, y0+1, x0+1] * flip
                    d100 = f_3d[z0+1, y0+1, x1+1] * (-flip)
                    d010 = f_3d[z0+1, y1+1, x0+1] * (-flip)
                    d001 = f_3d[z1+1, y0+1, x0+1] * (-flip)
                    d110 = f_3d[z0+1, y1+1, x1+1] * flip
                    d101 = f_3d[z1+1, y0+1, x1+1] * flip
                    d011 = f_3d[z1+1, y1+1, x0+1] * flip
                    d111 = f_3d[z1+1, y1+1, x1+1] * (-flip)
                    d00 = d000 * (1 - xd) + d100 * xd
                    d01 = d001 * (1 - xd) + d101 * xd
                    d10 = d010 * (1 - xd) + d110 * xd
                    d11 = d011 * (1 - xd) + d111 * xd
                    d0 = d00 * (1 - yd) + d10 * yd
                    d1 = d01 * (1 - yd) + d11 * yd
                    value = d0 * (1 - zd) + d1 * zd
                else
                    x0 = Int(round(x))
                    y0 = Int(round(y))
                    z0 = Int(round(z))
                    if (y0 < 0) y0 += size_3d end
                    if (z0 < 0) z0 += size_3d end

                    if (z0 + y0 + x0) % 2 == 1
                        flip = -1
                    else
                        flip = 1
                    end

                    value = f_3d[z0+1, y0+1, x0+1] * flip
                end

                if neg
                    f_2d[i+1, j+1] = conj(value)
                else
                    f_2d[i+1, j+1] = value
                end

                ph = h * x_shift / size_2d + k * y_shift / size_2d
                f_2d[i+1, j+1] *= exp(-2 * pi * ph * 1.0im)

                if ctf
                    ctf_value = ctf_2d_arr(h, k, size_2d,
                                           defocus_u, defocus_v, defocus_a,
                                           cs, apix, kv, ac, b_fac)
                    if white
                        f_2d[i+1, j+1] *= ctf_value
                    else
                        f_2d[i+1, j+1] *= -ctf_value
                    end
                end

                if (h + k) % 2 != 0
                    f_2d[i+1, j+1] *= -1
                end
            end
        end
    end
end
