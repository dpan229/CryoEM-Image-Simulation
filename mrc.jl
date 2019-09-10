using StaticArrays

"""
    Header(nc, nr, ns, mode, ncstart, nrstart, nsstart, nx, ny, nz, [...])

Struct for storing Mrc header data.
"""
mutable struct Header
    nc      ::Int32
    nr      ::Int32
    ns      ::Int32
    mode    ::Int32
    ncstart ::Int32
    nrstart ::Int32
    nsstart ::Int32
    nx      ::Int32
    ny      ::Int32
    nz      ::Int32
    x_length::Float32
    y_length::Float32
    z_length::Float32
    alpha   ::Float32
    beta    ::Float32
    gamma   ::Float32
    mapc    ::Int32
    mapr    ::Int32
    maps    ::Int32
    amin    ::Float32
    amax    ::Float32
    amean   ::Float32
    ispg    ::Int32
    nsymbt  ::Int32
    lskflg  ::Int32
    skwmat  ::SVector{9, Float32} # array length = 9
    skwtrn  ::SVector{3, Float32} # array length = 3
    extra   ::SVector{15, Int32}  # array length = 15
    mapi    ::SVector{4, Int8}    # array length = 4
    machst  ::SVector{4, Int8}    # array length = 4
    rms     ::Float32
    nlabl   ::Int32
    label_n ::SArray{Tuple{10, 80}, Int8} # dimensions = 10, 80
end

"""
    header_byteswap!(header)

Flip the endianness of every field of `header`.
"""
function header_byteswap!(header)
    for field in fieldnames(header)
        setfield!(header, field, bswap(getfield(header, field)))
    end
end

function sys_byteorder()
    return ENDIAN_BOM == 0x04030201 ? "little" : "big"
end

"""
    Mrc(header, data, apix, file)

Struct for storing and manipulating data from a .mrc file.
"""
struct Mrc
    header::Header
    data  ::Union{Array{Int8, 3},
                  Array{Int16, 3},
                  Array{Float32, 3},
                  Array{SVector{2, Int16}, 3},
                  Array{Complex{Float32}, 3}}
    apix  ::Float32
    file  ::String
    function Mrc(obj=nothing, apix=0.0)
        if obj === nothing
            print("[Mrc] No data! Header created!")
            header = _initialize_mrc_header()
            new(header, Array{Float32}(undef, 0, 0, 0), apix, "")
        elseif typeof(obj) == String
            # read from .mrc file
            f = open(obj, "r")
            header = _header_from_file(f)

            if header.machst == SVector{4, Int8}([17, 17, 0, 0]) # need to compare bytes?
                data_endianness = "big"
            else
                data_endianness = "little" # default
            end
            if data_endianness != sys_byteorder()
                header_byteswap!(header)
            end

            if header.mode in [0, 1, 2, 3, 4]
                data_type = [Int8, Int16, Float32, SVector{2, Int16}, Complex{Float32}][header.mode+1] # Array length = 2
            else
                error("[Mrc] Error: unrecognized data mode")
            end
            # update apix
            if header.x_length != 0 && header.y_length != 0 && header.z_length != 0
                out_apix = float(header.x_length / header.nc)
            else
                out_apix = apix
            end
            seek(f, 1024)
            temp_data = Vector{data_type}(undef, div(filesize(obj) - 1024, sizeof(data_type)))
            read!(f, temp_data)
            if header.ns != 0
                # Julia's reshape returns matrix with old entries in column-major order, so we
                # give the new dimensions to reshape in reverse order and permute the dimensions
                temp_data = permutedims(reshape(temp_data, (header.nc, header.nr, header.ns)), [3, 2, 1])
            else
                temp_data = permutedims(reshape(temp_data, (header.nc, header.nr)), [2, 1])
            end
            # swap data if endianness is different
            if data_endianness != sys_byteorder()
                header_byteswap!(header)
            end
            close(f)
            new(header, temp_data, out_apix, obj)
        elseif typeof(obj) <: Array
            # initialize with array
            if ndims(obj) == 2
                obj = reshape(obj, 1, size(obj)...) # turn 2D axb array into 3D 1xaxb array
            else
                if ndims(obj) != 3
                    error("[Mrc] Error: only takes 2 or 3 dimensional data")
                end
            end
            out = new(_initialize_mrc_header(), obj, apix, "")
            update_header!(out)
            return out
        else
            error("[Mrc] Error: construction requires string for filename or multidimensional array for data")
        end
    end
end

function _initialize_mrc_header()
    return Header(vcat(
        zeros(Int32,  10),
        zeros(Float32, 3),
        Float32(90), # alpha
        Float32(90), # beta
        Float32(90), # gamma
        Int32(1),  # mapc
        Int32(2),  # mapr
        Int32(3),  # maps
        zeros(Float32, 3),
        Int32(1),  # ispg
        zeros(Int32,   2),
        [SVector{9}(zeros(Float32, 9))],
        [SVector{3}(zeros(Float32, 3))],
        [SVector{15}(zeros(Int32, 15))],
        [SVector{4, Int8}([77, 65, 80, 32])], # mapi = fromstring(b'MAP ')
        [sys_byteorder == "little" ? SVector{4, Int8}([68, 65, 0, 0]) : SVector{4, Int8}([17, 17, 0, 0])], # machst varies depending on endianness
        zeros(Float32, 1),
        zeros(Int32,   1),
        [SArray{Tuple{10, 80}}(zeros(Int8, 10, 80))])...)
end

function _header_from_file(file)
    args1      = zeros(Int32,  10) # nc, ..., nz
    args2      = zeros(Float32, 6) # x_length, ..., gamma
    args3      = zeros(Int32,   3) # mapc, mapr, maps
    args4      = zeros(Float32, 3) # amin, amax, amean
    args5      = zeros(Int32,   3) # ispg, nsymbt, lskflg
    to_skwmat  = zeros(Float32, 9)
    to_skwtrn  = zeros(Float32, 3)
    to_extra   = zeros(Int32,  15)
    to_mapi    = zeros(Int8,    4)
    to_machst  = zeros(Int8,    4)
    args6      = zeros(Float32, 1) # rms
    args7      = zeros(Int32,   1) # nlabl
    to_label_n = zeros(Int8,   10, 80)
    for c in [args1, args2, args3, args4, args5, to_skwmat, to_skwtrn, to_extra, to_mapi, to_machst, args6, args7, to_label_n]
        read!(file, c)
    end
    return Header(vcat(
            args1, args2, args3, args4, args5,
            [SVector{9, Float32}(to_skwmat)],
            [SVector{3, Float32}(to_skwtrn)],
            [SVector{15, Int32}(to_extra)],
            [SVector{4, Int8}(to_mapi)],
            [SVector{4, Int8}(to_machst)],
            args6, args7,
            [SArray{Tuple{10, 80}, Int8}(to_label_n)])...)
end

"""
    headerset!(mrc, field, val)

Set the given field of the header of `mrc` to val and print the old and new
values if the new value is different.
"""
function headerset!(mrc::Mrc, field::Symbol, val)
    old = getfield(mrc.header, field)
    if old != val
        println("Mrc: ", rpad(field, 16), "NEW: ", rpad(val, 16), "OLD: ", old)
        setfield!(mrc.header, field, typeof(val) == Int64 ? Int32(val) : val)
    end
end

"""
    update_header!(mrc)

Update data in the header of `mrc` to match information about its data.
"""
function update_header!(mrc::Mrc)
    if length(mrc.data) == 0
        error("[Mrc] no data to update header!")
    end
    headerset!(mrc, :nc, size(mrc.data)[3])
    headerset!(mrc, :nr, size(mrc.data)[2])
    headerset!(mrc, :ns, size(mrc.data)[1]) # 1 if data is 2 dimensional. Data always stored in 3D array

    t = eltype(mrc.data)
    headerset!(mrc, :mode, t == Int8    ? 0 :
                           t == Int16   ? 1 :
                           t == Float32 ? 2 :
                           t == SVector{2, Int16} ? 3 :
                           t == Complex{Float32}  ? 4 : error("[Mrc] data type not recognized!"))
    # determine nx, ny, nz
    if size(mrc.data)[1] == size(mrc.data)[2] == size(mrc.data)[3]
        headerset!(mrc, :nx, size(mrc.data)[3])
        headerset!(mrc, :ny, size(mrc.data)[2])
        headerset!(mrc, :nz, size(mrc.data)[1])
    else # image stack or 2D image
        headerset!(mrc, :nx, size(mrc.data)[2])
        headerset!(mrc, :ny, size(mrc.data)[1])
        headerset!(mrc, :nz, 1)
    end
    # determine amin, amax, amean, rms
    headerset!(mrc, :amin,  minimum(mrc.data)) # will this fit the data type of amin?
    headerset!(mrc, :amax,  maximum(mrc.data))
    headerset!(mrc, :amean, sum(mrc.data) / length(mrc.data))
    headerset!(mrc, :rms,   sqrt(sum([((x - mrc.header.amean) / length(mrc.data))^2 for x in mrc.data])))
    # determine x_length, y_length, z_length
    if mrc.apix != 0
        headerset!(mrc, :x_length, mrc.apix * size(mrc.data)[3])
        headerset!(mrc, :y_length, mrc.apix * size(mrc.data)[2])
        headerset!(mrc, :z_length, size(mrc.data)[1] == 1 ? 0 : mrc.apix * size(mrc.data)[1])
    end
end

function write_mrc(mrc, file_name)
    update_header!(mrc)
    file = open(file_name, "w")
    for field in fieldnames(Header)
        write(file, getfield(mrc.header, field))
    end
    write(file, mrc.data)
    close(file)
end
