include("./mrc.jl")
include("./star.jl")
include("./projection.jl")

using FFTW: rfft, irfft

"""
    generate_projections(star, map_mrc, pad, res, apix; trilinear=true, white=true)

Generate a stack of projections (stored in an Mrc object) specified by `star`
of the map `map_mrc`. `star` and `map_mrc` can either be RelionStar and Mrc
objects or paths to .star and .mrc files.
"""
function generate_projections(star, map_mrc, pad, res, apix; trilinear=true, white=true)
    if typeof(star) == String
        star = RelionStar(star)
    end
    if typeof(map_mrc) == String
        map_mrc = Mrc(map_mrc, apix)
    end

    original_size = size(map_mrc.data)[1]
    #println("Original size = ", original_size)
    #println("Finished reading the volume")

    phi, theta, psi, x, y  = read_relion_transformation(star)
    du, dv, da, cs, kv, ac = read_relion_ctf(star)
    cs = cs === nothing ? 0.0 : cs[1]
    kv = kv === nothing ? 200 : kv[1]
    ac = ac === nothing ? 0.0 : ac[1]
    #println("Finished reading the star file")

    f_3d = preprocess_for_projection(map_mrc.data, pad, trilinear)
    #println("Finished fourier tranform of 3D map")

    stack_data = zeros(Float32, length(x), original_size, original_size)
    for i in 1:length(x)
        f_2d = zeros(Complex{Float64}, original_size, div(original_size, 2) + 1)
        projection!(f_2d, f_3d,
                    phi[i], theta[i], psi[i], -x[i], -y[i],
                    pad, true,
                    du[i], dv[i], da[i], cs,
                    apix, kv, ac, -1, white, trilinear)

        #applywhitenoise!(f_2d, 0.1, rms(f_2d)) # apply white noise in frequency domain

        stack_data[i, :, :] = irfft(f_2d, original_size, [2, 1])
        #println("Finished projecting particle $(i)/$(length(x))")
    end

    return Mrc(stack_data, apix)
end

whitenoise_rms(snr, signal_rms) = signal_rms / sqrt(snr)
avg(x) = sum(x) / length(x)
rms(x) = sqrt(avg(map(x -> x^2, x)))

"""
    applywhitenoise!(freq_signal, snr, freq_signal_rms)

Apply white noise to a signal in frequency space `freq_signal` (with rms
`freq_signal_rms`), with a signal-to-noise ratio `snr`.
"""
function applywhitenoise!(freq_signal, snr, freq_signal_rms)
    freq_signal[:] = freq_signal + whitenoise_rms(snr, freq_signal_rms) * randn(size(freq_signal))
end

star_file = "test_data/70S.star"
map_file = "test_data/70S_1_r1.mrc"

pad = 2
res = 1.0
apix = 2.52
trilinear = true
white = true

projs = generate_projections(star_file, map_file, 2, 1.0, 2.52)
write_mrc(projs, "test_output.mrc")
