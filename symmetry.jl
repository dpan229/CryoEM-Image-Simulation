pg_dict = Dict("CI"  => 200,
               "CS"  => 201,
               "CN"  => 202,
               "CNV" => 203,
               "CNH" => 204,
               "SN"  => 205,
               "DN"  => 206,
               "DNV" => 207,
               "DNH" => 208,
               "T"   => 209,
               "TD"  => 210,
               "TH"  => 211,
               "O"   => 212,
               "OH"  => 213,
               "I"   => 214,
               "IH"  => 215,

               "I1"  => 216,
               "I2"  => 217,
               "I3"  => 218,
               "I4"  => 219,
               "I5"  => 220,

               "I1H" => 221,
               "I2H" => 222,
               "I3H" => 223,
               "I4H" => 224,
               "I5H" => 225)

const accuracy = 1e-6

mutable struct SymmetryList#{T<:AbstractFloat}
    L::Array{Float64,2}
    R::Array{Float64,2}
    chain_length::Array{Int64,1}
    true_sym_no ::Int64
    sym_elements::Int64
end

"""
    empty_symmetry_list()

Create an empty SymmetryList object.
"""
empty_symmetry_list() = SymmetryList(Array{Float64}(undef, 0, 4),
                                     Array{Float64}(undef, 0, 4),
                                     Array{Int64}(undef, 0),
                                     0, 0)

"""
    align_with_z(axis)

Compute the 4x4 matrix A ([X X X 0; X X X 0; X X X 0; 0 0 0 1]) that rotates
the given 3-element axis onto the Z axis.
If axis = [x, y, z], then align_with_z(axis) * [x; y; z; 0] == [0; 0; norm(axis); 0].
"""
function align_with_z(axis)
    out = zeros(4, 4)
    out[4, 4] = 1
    norm_axis = norm(axis) > accuracy ? normalize(axis) : zeros(3)
    # length of projection of norm_axis onto YZ plane
    proj_mod = sqrt(norm_axis[2]^2 + norm_axis[3]^2)
    if proj_mod > accuracy
        out[1, 1] = proj_mod
        out[1, 2] = -norm_axis[1] * norm_axis[2] / proj_mod
        out[1, 3] = -norm_axis[1] * norm_axis[3] / proj_mod
        out[2, 2] =  norm_axis[3] / proj_mod
        out[2, 3] = -norm_axis[2] / proj_mod
        out[3, 1:3] = norm_axis
    else
        # proj_mod == 0 => axis is either positive or negative x axis
        out[1, 3] = norm_axis[1] > 0 ? -1 : 1
        out[2, 2] = 1
        out[3, 1] = norm_axis[1] > 0 ? 1 : -1
    end
    return out
end

"""
    rotation_3D_matrix(ang, axis)

Give the 4x4 matrix corresponding to a 3D rotation about the given axis
by the given angle in degrees. The axis `axis` can be a character ('X', 'Y'
or 'Z') or a 3-element Float64 array.
"""
function rotation_3D_matrix(ang, axis::Char)
    result = zeros(4, 4)
    result[4, 4] = 1

    ang = ang * pi / 180
    cosine = cos(ang)
    sine = sin(ang)
    if axis == 'Z'
        result[1, 1] = cosine
        result[1, 2] = sine
        result[2, 1] = -sine
        result[2, 2] = cosine
        result[3, 3] = 1
    elseif axis == 'Y'
        result[1, 1] = cosine
        result[1, 3] = sine
        result[3, 1] = -sine
        result[3, 3] = cosine
        result[2, 2] = 1
    elseif axis == 'X'
        result[2, 2] = cosine
        result[2, 3] = sine
        result[3, 2] = -sine
        result[3, 3] = cosine
        result[1, 1] = 1
    else
        error("[rotation_3D_matrix] error: axis $axis is invalid (should be X, Y or Z)")
    end
    return result
end

function rotation_3D_matrix(ang, axis::Array{Float64, 1})
    result = Array{Float64}(I, 4, 4)
    ang = -(ang * pi / 180)
    sine = sin(ang)
    cosine = cos(ang)
    cosine_1 = 1 - cosine
    x, y, z = axis
    result[1, 1] = cosine + x^2 * cosine_1
    result[1, 2] = x*y * cosine_1 - z * sine
    result[1, 3] = x*z * cosine_1 + y * sine
    result[2, 1] = x*y * cosine_1 + z * sine
    result[2, 2] = cosine + y^2 * cosine_1
    result[2, 3] = y*z * cosine_1 - x * sine
    result[3, 1] = x*z * cosine_1 - y * sine
    result[3, 2] = y*z * cosine_1 + x * sine
    result[3, 3] = cosine + z^2 * cosine_1
    return result
end



"""
    read_sym_or_file(fname_sym)

Read a symmetry file with the given filename and return the resulting SymmetryList
object. Alternatively takes a string representing a symmetry group ("C4", "T", etc.),
returning the resulting SymmetryList object.
"""
function read_sym_or_file(fname_sym)
    true_sym_no = 0
    axis_no = 0
    mirror_plane_no = 0
    inversion_pt_no = 0
    sym_elements = 0
    axis = Array{Float64}(undef, 3)
    file_contents = String[]
    try
        file_contents = collect(eachline(fname_sym))
    catch
        pgGroup, pgOrder = symmetry_group_and_order(fname_sym)
        if pgGroup != -1
            file_contents = fill_symmetry(fname_sym, pgGroup, pgOrder)
        else
            error("[read_sym_or_file] error: could not open file or unrecognized symmetry group $fname_sym")
        end
    end
    for line in file_contents
        (line[1] == ';' || line[1] == '#') && continue
        words = split(line)
        if words[1] == "rot_axis"
            fold = parse(Int64, words[2])
            true_sym_no += fold - 1
            axis_no += 1
        elseif words[1] == "mirror_plane"
            true_sym_no += 1
            mirror_plane_no += 1
        elseif words[1] == "inversion"
            true_sym_no += 1
            inversion_pt_no += 1
        end
    end
    L = Array{Float64}(undef, 4 * true_sym_no, 4)
    R = Array{Float64}(undef, 4 * true_sym_no, 4)
    chain_length = ones(Int64, true_sym_no)

    i = 0
    for line in file_contents
        (line[1] == ';' || line[1] == "#") && continue
        words = split(line)
        if words[1] == "rot_axis"
            fold = parse(Int64, words[2])
            map!(x -> parse(Float64, x), axis, words[3:5])
            this_L = Array{Float64}(I, 4, 4)
            #this_R = Array{Float64}(undef, 4, 4)
            for j in 1:fold-1
                this_R = rotation_3D_matrix(j * 360 / fold, axis)
                small_vals_to_zero!(this_R)
                set_matrices_mtrx!(L, R, i, this_L, permutedims(this_R, (2, 1)))
                i += 1
            end
            sym_elements += 1
        elseif words[1] == "inversion"
            this_L = Array{Float64}(I, 4, 4)
            this_R = Array{Float64}(I, 4, 4)
            this_L[3, 3] = this_R[1, 1] = this_R[2, 2] = this_R[3, 3] = -1.
            set_matrices_mtrx!(L, R, i, this_L, this_R)
            i += 1
            sym_elements += 1
        elseif words[1] == "mirror_plane"
            map!(x -> parse(Float64, x), axis, words[2:4])
            this_L = Array{Float64}(I, 4, 4)
            this_L[3, 3] = -1.
            A = align_with_z(axis)
            A = permutedims(A, (2, 1))
            this_R = A * this_L * inv(A)
            this_L = Array{Float64}(I, 4, 4)
            set_matrices_mtrx!(L, R, i, this_L, this_R)
            i += 1
            sym_elements += 1
        end
    end
    out = SymmetryList(L, R, chain_length, true_sym_no, sym_elements)
    compute_subgroup!(out)
    return out
end

"""
    sym_no(symList)

Give the number of symmetries present in a SymmetryList object.
"""
sym_no(symList::SymmetryList) = size(symList.L)[1] >> 2

# get_matrices and set_matrices! are 0-indexed

"""
    get_matrices(symList, i)

Return the i'th (0-indexed) L and R matrices in a SymmetryList object in a tuple.
"""
function get_matrices(symList::SymmetryList, i)
    return view(symList.L, 4i+1:4i+4, :), view(symList.R, 4i+1:4i+4, :)
end

function set_matrices_mtrx!(L, R, i, new_L, new_R)
    L[4i+1:4i+4, :] = new_L
    R[4i+1:4i+4, :] = new_R
end

"""
    set_matrices!(symList, i, new_L, new_R)

Overwrite the i'th (0-indexed) L and R matrices in a SymmetryList object with
new_L and new_R.
"""
set_matrices!(symList::SymmetryList, i, new_L, new_R) = set_matrices_mtrx!(symList.L, symList.R, i, new_L, new_R)

"""
    add_matrices!(symList, new_L, new_R, new_chain_length)

Append the given symmetry information to the data of a SymmetryList object, adding
a new symmetry.
"""
function add_matrices!(symList::SymmetryList, new_L, new_R, new_chain_length)
    (size(new_L) == (4, 4) && size(new_R) == (4, 4)) || error("[Symmetries] error: transformation matrix not 4x4")
    if symList.true_sym_no == sym_no(symList)
        symList.L = vcat(symList.L, new_L)
        symList.R = vcat(symList.R, new_R)
        push!(symList.chain_length, new_chain_length)
    else
        set_matrices!(symList, symList.true_sym_no, new_L, new_R)
        symList.chain_length[end] = new_chain_length
    end
    symList.true_sym_no += 1
end

# takes paths through `tried` like this:
# XXXvX...
# XXXvX...
# XXXvX...
# <<<<X...
# XXXXX...
# ...
function find_not_tried(tried, true_sym_no)
    i = j = n = 1
    while n <= size(tried)[1]
        if tried[i, j] == 0 && !(i > true_sym_no && j > true_sym_no)
            return (i, j)
        end
        if i != n
            i += 1 # move down
        else
            j -= 1 # move left
            if j == 0 # layer finished; go to next layer
                n += 1
                j = n
                i = 1
            end
        end
    end
    return (-1, -1)
end

"""
    matrix_equal(m1, m2)

Test if two matrices are equal within the global accuracy. m1 and m2 must have the
same dimensions. Return true if every element of m1 is within accuracy of the
corresponding element of m2.
"""
matrix_equal(m1, m2) = maximum(map(abs, m1 - m2)) < accuracy

"""
    is_identity(mtrx)

Test if the given matrix is equal to the idenity matrix within the global accuracy.
"""
is_identity(mtrx) = matrix_equal(mtrx, Matrix{Float64}(I, size(mtrx)))

"""
    small_vals_to_zero!(mtrx)

Set all values that are within the global accuracy of 0 to 0 in the given matrix.
"""
small_vals_to_zero!(mtrx) = map!(x -> abs(x) < accuracy ? 0. : x, mtrx, mtrx)

"""
    compute_subgroup!(symList)

Add new symmetries to a SymmetryList object until it is a complete subgroup.
These symmetries are generated by multiplying its symmetry matrices until no new
matrices can be produced.
"""
function compute_subgroup!(symList::SymmetryList)
    ##Id = Matrix{Float64}(I, 4, 4)
    #L1 = Matrix{Float64}(undef, 4, 4)
    #R1 = Matrix{Float64}(undef, 4, 4)
    #L2 = Matrix{Float64}(undef, 4, 4)
    #L2 = Matrix{Float64}(undef, 4, 4)
    #new_L = Matrix{Float64}(undef, 4, 4)
    #new_R = Matrix{Float64}(undef, 4, 4)
    tried = zeros(Int64, symList.true_sym_no, symList.true_sym_no)
    while ((i, j) = find_not_tried(tried, symList.true_sym_no)) != (-1, -1)
        tried[i, j] = 1
        L1, R1 = get_matrices(symList, i-1)
        L2, R2 = get_matrices(symList, j-1)
        new_L = L1 * L2
        new_R = R1 * R2
        new_chain_length = symList.chain_length[i] + symList.chain_length[j]
        (is_identity(new_L) && is_identity(view(new_R, 1:3, 1:3))) && continue

        found = false
        for n in 0:sym_no(symList)-1
            L1, R1 = get_matrices(symList, n)
            if matrix_equal(new_L, L1) && matrix_equal(new_R, R1)
                found = true
                break
            end
        end
        if !found
            small_vals_to_zero!(new_L)
            small_vals_to_zero!(new_R)
            add_matrices!(symList, new_L, new_R, new_chain_length)
            # expand `tried` in both dimensions by 1
            tried = hcat(vcat(tried, zeros(1, size(tried)[2])), zeros(size(tried)[1] + 1, 1))
        end
    end
end

"""
    symmetry_group_and_order(sym)

Give the group number and order of the symmetry group represented by the given
string. If the string is not a valid symmetry group, return (-1, -1).
"""
# returns (pgGroup, pgOrder) for the symmetry group represented by the input string,
# returns (-1, -1) if unsuccessfull
function symmetry_group_and_order(sym)
    len = length(sym)
    (len < 1 || len > 4) && return (-1, -1)
    sym = uppercase(sym)
    try
        return (pg_dict[sym], -1)
    catch
        if len >= 2 && isdigit(sym[2])
            if len >= 3 && isdigit(sym[3])
                order = parse(Int64, sym[2:3])
                sym = sym[1] * "N" * sym[4:end]
            else
                order = parse(Int64, sym[2])
                sym = sym[1] * "N" * sym[3:end]
            end
            try
                return (pg_dict[sym], order)
            catch
                return (-1, -1)
            end
        end
    end
    return (-1, -1)
end

"""
    fill_symmetry(sym, pgGroup, pgOrder)

Give symmetry information (as an array of strings) for the given symmetry
group number and order.
"""
function fill_symmetry(sym, pgGroup, pgOrder)
    content = Array{String}(undef, 0)
    pgName = findfirst(==(pgGroup), pg_dict)
    if pgName == "CN"
        push!(content, "rot_axis $pgOrder 0 0 1")
    elseif pgName == "CI"
        push!(content, "inversion ")
    elseif pgName == "CS"
        push!(content, "mirror_plane 0 0 1")
    elseif pgName == "CNV"
        push!(content, "rot_axis $pgOrder 0 0 1",
                       "mirror_plane 0 1 0")
    elseif pgName == "CNH"
        push!(content, "rot_axis $pgOrder 0 0 1",
                       "mirror_plane 0 0 1")
    elseif pgName == "SN"
        isodd(pgOrder) && error("[Symmetries] error: order of SN group must be even")
        push!(content, "rot_axis $(pgOrder >> 1) 0 0 1",
                       "inversion ")
    elseif pgName == "DN"
        push!(content, "rot_axis $pgOrder 0 0 1",
                       "rot_axis 2 1 0 0")
    elseif pgName == "DNV"
        push!(content, "rot_axis $pgOrder 0 0 1",
                       "rot_axis 2 1 0 0",
                       "mirror_plane 1 0 0")
    elseif pgName == "DNH"
        push!(content, "rot_axis $pgOrder 0 0 1",
                       "rot_axis 2 1 0 0",
                       "mirror_plane 0 0 1")
    elseif pgName == "T"
        push!(content, "rot_axis 3 0. 0. 1.",
                       "rot_axis 2 0. 0.8164965809 0.5773502692")
    elseif pgName == "TD"
        push!(content, "rot_axis 3 0. 0. 1.",
                       "rot_axis 2 0. 0.8164965809 0.5773502692",
                       "mirror_plane 1.4142136 2.4494897 0.")
    elseif pgName == "O"
        push!(content, "rot_axis 3 .5773502 .5773502 .5773502",
                       "rot_axis 4 0 0 1")
    elseif pgName == "OH"
        push!(content, "rot_axis 3 .5773502 .5773502 .5773502",
                       "rot_axis 4 0 0 1",
                       "mirror_plane 0 1 1")
    elseif pgName == "I" || pgName == "I2"
        push!(content, "rot_axis 2 0 0 1",
                       "rot_axis 5 0.525731114 0. 0.850650807",
                       "rot_axis 3 0. 0.356822076 0.934172364")
    elseif pgName == "I1"
        push!(content, "rot_axis 2 1 0 0",
                       "rot_axis 5 0.85065080702670 0. -0.5257311142635",
                       "rot_axis 3 0.9341723640 0.3568220765 0.")
    elseif pgName == "I3"
        push!(content, "rot_axis 2 -0.5257311143 0. 0.8506508070",
                       "rot_axis 5 0. 0. 1.",
                       "rot_axis 3 -0.4911234778630044 0.3568220764705179 0.7946544753759428")
    elseif pgName == "I4"
        push!(content, "rot_axis 2 0.5257311143 0. 0.8506508070",
                       "rot_axis 5 0.8944271932547096 0. 0.4472135909903704",
                       "rot_axis 3 0.4911234778630044 0.3568220764705179 0.7946544753759428")
    elseif pgName == "I5"
        # ??
    elseif pgName == "IH" || pgName == "I2H"
        push!(content, "rot_axis 2 0 0 1",
                       "rot_axis 5 0.525731114 0. 0.850650807",
                       "rot_axis 3 0. 0.356822076 0.934172364",
                       "mirror_plane 1 0 0")
    elseif pgName == "I1H"
        push!(content, "rot_axis 2 1 0 0",
                       "rot_axis 5 0.85065080702670 0. -0.5257311142635",
                       "rot_axis 3 0.9341723640 0.3568220765 0.",
                       "mirror_plane 0 0 -1")
    elseif pgName == "I3H"
        push!(content, "rot_axis 2 -0.5257311143 0. 0.8506508070",
                       "rot_axis 5 0. 0. 1.",
                       "rot_axis 3 -0.4911234778630044 0.3568220764705179 0.7946544753759428",
                       "mirror_plane 0.850650807 0. 0.525731114")
    elseif pgName == "I4H"
        push!(content, "rot_axis 2 0.5257311143 0. 0.8506508070",
                       "rot_axis 5 0.8944271932547096 0. 0.4472135909903704",
                       "rot_axis 3 0.4911234778630044 0.3568220764705179 0.7946544753759428",
                       "mirror_plane 0.850650807 0. -0.525731114")
    elseif pgName == "I5H"
        # ??
    else
        error("[Symmetries] error: Symmetry $sym (pgGroup $pgGroup) is not known")
    end
end

"""
    write_definition(io, fname_sym)

Write the definition and rotation matrices of the given symmetry file or group
name to `io`, an i/o stream.
"""
function write_definition(io::IO, fname_sym)
    sym = read_sym_or_file(fname_sym)
    write(io, "++++ Using symmetry group $fname_sym with the following $(sym_no(sym)+1) transformation matrices:\n")
    R = Array{Float64}(I, 3, 3)
    write(io, " R(1) = $R\n")
    for i in 0:sym_no(sym)-1
        L, R = get_matrices(sym, i)
        if !is_identity(view(L, 1:3, 1:3))
            write(io, " L($(i+2)) = $(view(L, 1:3, 1:3))\n")
        end
        write(io, " R($(i+2)) = $(view(R, 1:3, 1:3))\n")
    end
end

"""
    non_redundant_ewald_sphere(pgGroup, pgOrder)

Give the area of the non redundant part of the Ewald sphere for the given
symmetry group number and order.
"""
function non_redundant_ewald_sphere(pgGroup, pgOrder)
    pgName = findfirst(==(pgGroup), pg_dict)
    if pgName == "CN" || pgName == "SN"
        return 4*pi/pgOrder
    elseif pgName == "CI" || pgName == "CS"
        return 2*pi
    elseif pgName == "CNV" || pgName == "CNH" || pgName == "DN"
        return 2*pi/pgOrder
    elseif pgName == "DNV" || pgName == "DNH"
        return pi/pgOrder
    elseif pgName == "T"
        return pi/3
    elseif pgName == "TD" || pgName == "TH" || pgName == "O"
        return pi/6
    elseif pgName == "OH"
        return pi/12
    elseif pgName == "I" || pgName == "I1" || pgName == "I2" || pgName == "I3" || pgName == "I4" || pgName == "I5"
        return pi/15
    elseif pgName == "IH" || pgName == "I1H" || pgName == "I2H" || pgName == "I3H" || pgName == "I4H" || pgName == "I5H"
        return pi/30
    else
        error("[non_redundant_ewald_sphere] error: symmetry group $pgGroup not known")
    end
end

periodic_mod(a, b) = (a % b) + (a % b < 0 ? b : 0)

"""
    real_wrap(x, mn, mx)

Wrap `x` to be between `mn` and `mx`
"""
real_wrap(x, mn, mx) = mn + periodic_mod(x - mn, mx - mn)

"""
    copy!(dest, src)

Copy the contents of the array `src` to `dest` in place. The length of `src`
should not be greater than the length of `dest`.
"""
function copy!(dest, src)
    for (i, x) in enumerate(src)
        dest[i] = x
    end
end

"""
    zeros!(dest)

Set the contents of the array `dest` to all zeros in place.
"""
function zeros!(dest)
    for i in 1:length(dest)
        dest[i] = 0
    end
end

"""
    apply_geometry!(out, img, At, inv, wrap)

Apply the 3D transformation represented by the 4x4 matrix `At` to the 3D map
`img` and save the result in `out` in place. `out` will not be resized, so
only the part of the result that fits in `out` will be saved there. Spline
degree is assumed to be 1 (linear interpolation), and the background (outside)
value is assumed to be 0.
"""
function apply_geometry!(out::Array{Float32, 3}, img::Array{Float32, 3},
                         At::Array{Float32, 2}, inv::Bool, wrap::Bool)
    # XSIZE(img) = size(img)[3], YSIZE(img) = size(img)[2], ZSIZE(img) = size(img)[1]
    out_zdim, out_ydim, out_xdim = size(out)
    img_zdim, img_ydim, img_xdim = size(img)
    if is_identity(At) && (out_xdim == 0 || size(out) == size(img))
        copy!(out, img)
        return
    end
    #if img_xdim == 0
        #out = Array{Float64, 3}(undef, 0, 0, 0)
    #    return
    #end

    if (!inv)
        A = inv(At)
    else
        A = At
    end

    #if out_xdim == 0
        #out = Array{Float64}(undef, size(img))
    #end

    zeros!(out)

    cen_z, cen_y, cen_x = (x -> div(x, 2)).(size(out))
    cen_zp, cen_yp, cen_xp = (x -> div(x, 2)).(size(img))
    min_xp, min_yp, min_zp = -cen_xp, -cen_yp, -cen_zp
    max_xp = img_xdim - cen_xp - 1
    max_yp = img_ydim - cen_yp - 1
    max_zp = img_zdim - cen_zp - 1

    for k in 0:out_zdim-1
        for i in 0:out_ydim-1
            x = -cen_x
            y = i - cen_y
            z = k - cen_z

            xp = x*A[1, 1] + y*A[1, 2] + z*A[1, 3] + A[1, 4]
            yp = x*A[2, 1] + y*A[2, 2] + z*A[2, 3] + A[2, 4]
            zp = x*A[3, 1] + y*A[3, 2] + z*A[3, 3] + A[3, 4]

            for j in 0:out_xdim-1
                interp = true
                x_out = xp < min_xp - accuracy || xp > max_xp + accuracy
                y_out = yp < min_yp - accuracy || yp > max_yp + accuracy
                z_out = zp < min_zp - accuracy || zp > max_zp + accuracy
                if wrap
                    x_out && (xp = real_wrap(xp, min_xp - 0.5, max_xp + 0.5))
                    y_out && (yp = real_wrap(yp, min_yp - 0.5, max_yp + 0.5))
                    z_out && (zp = real_wrap(zp, min_zp - 0.5, max_zp + 0.5))
                else
                    (x_out || y_out || z_out) && (interp = false)
                end
                xp = max(xp, min_xp)
                yp = max(yp, min_yp)
                zp = max(zp, min_zp)
                if interp
                    #println("yp=$yp min_yp=$min_yp max_yp=$max_yp")
                    wx = xp + cen_xp
                    m1 = floor(Int, wx)
                    wx -= m1
                    m2 = m1 + 1
                    wy = yp + cen_yp
                    n1 = floor(Int, wy)
                    wy -= n1
                    n2 = n1 + 1
                    wz = zp + cen_zp
                    q1 = floor(Int, wz)
                    wz -= q1
                    q2 = q1 + 1
                    wx_1 = 1 - wx
                    wy_1 = 1 - wy
                    wz_1 = 1 - wz

                    aux1 = wz_1 * wy_1
                    aux2 = aux1 * wx_1
                    tmp = aux2 * img[q1+1, n1+1, m1+1]
                    if wx != 0 && m2 < img_xdim
                        tmp += (aux1 - aux2) * img[q1+1, n1+1, m2+1]
                    end
                    if wy != 0 && n2 < img_ydim
                        aux1 = wz_1 * wy
                        aux2 = aux1 * wx_1
                        tmp += aux2 * img[q1+1, n2+1, m1+1]
                        if wx != 0 && m2 < img_xdim
                            tmp += (aux1 - aux2) * img[q1+1, n2+1, m2+1]
                        end
                    end
                    if wz != 0 && q2 < img_zdim
                        aux1 = wz * wy_1
                        aux2 = aux1 * wx_1
                        tmp += aux2 * img[q2+1, n1+1, m1+1]
                        if wx != 0 && m2 < img_xdim
                            tmp += (aux1 - aux2) * img[q2+1, n1+1, m2+1]
                        end
                        if wy != 0 && n2 < img_ydim
                            aux1 = wz * wy
                            aux2 = aux1 * wx_1
                            tmp += aux2 * img[q2+1, n2+1, m1+1]
                            if wx != 0 && m2 < img_xdim
                                tmp += (aux1 - aux2) * img[q2+1, n2+1, m2+1]
                            end
                        end
                    end
                    out[k+1, i+1, j+1] = tmp
                else
                    out[k+1, i+1, j+1] = 0
                end
                xp += A[1, 1]
                yp += A[2, 1]
                zp += A[3, 1]
            end
        end
    end
end

"""
    symmetrise_map!(img, fname_sym, wrap)

Symmetrise the 3D map `img` according to the symmetry file or group name
`fname_sym`.
"""
function symmetrise_map!(img::Array{Float32, 3}, fname_sym, wrap::Bool)
    sym = read_sym_or_file(fname_sym)
    sum = copy(img)
    aux = Array{Float32, 3}(undef, size(img))
    R_copy = Array{Float32}(undef, 4, 4)
    for i in 0:sym_no(sym)-1
        L, R = get_matrices(sym, i)
        copy!(R_copy, R)
        apply_geometry!(aux, img, R_copy, true, wrap)
        sum += aux
    end
    sum /= sym_no(sym) + 1
    copy!(img, sum)
end
