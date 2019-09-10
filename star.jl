using DataStructures: OrderedDict

"""
    RelionStar(data_dict, block_name_list, key_name_list, star_type, file)

Struct to store and manipulate the data held in a .star file. Can be initialized
with no arguments (creates an empty object), a filename, or a dictionary.
"""
struct RelionStar
    data_dict::OrderedDict{String, OrderedDict{String, Array{Union{Int64, Float64}, 1}}}
    block_name_list::Array{String, 1}
    key_name_list::Array{Array{String, 1}, 1} # Lists may be different lengths
    star_type::String
    file::String
    function RelionStar(args...)
        this_data_dict = OrderedDict()
        this_block_name_list = []
        this_key_name_list = []
        this_star_type = ""
        this_file = ""
        if length(args) != 1
            print("[RelionStar] Empty object created")
            new(this_data_dict, this_block_name_list, this_key_name_list, this_star_type, this_file)
        elseif typeof(args[1]) == String
            f = open(args[1], "r")
            this_file = args[1]

            temp_dict = OrderedDict()
            temp_list = []
            block_name = ""
            is_in_loop = false
            is_parsing = false
            # initialize variables
            key_name_list_index = 0
            strp_line = ""
            for line in readlines(f)
                strp_line = strip(line)
                if startswith(strp_line, "data")
                    block_name = strp_line
                    push!(this_block_name_list, block_name)
                    push!(this_key_name_list, [])
                    key_name_list_index = length(this_key_name_list)
                elseif strp_line == "loop_"
                    is_in_loop = true
                elseif startswith(strp_line, "_")
                    is_parsing = true
                    if is_in_loop
                        temp_dict[strp_line] = []
                        push!(temp_list, strp_line)
                        push!(this_key_name_list[key_name_list_index], strp_line)
                    else
                        data_line = split(line)
                        temp_dict[data_line[1]] = _convert_string_to_number(data_line[2])
                        push!(this_key_name_list[key_name_list_index], data_line[1])
                    end
                elseif strp_line == ""
                    if is_parsing # ignore blank lines when not parsing
                        # skipping over conversion to np array here
                        this_data_dict[block_name] = temp_dict
                        is_parsing = false
                        is_in_loop = false
                        temp_dict = OrderedDict()
                        temp_list = []
                    end
                elseif startswith(strp_line, "#")
                    nothing # ignore comments
                else
                    data_line = split(line)
                    for i in 1:length(data_line)
                        push!(temp_dict[temp_list[i]], _convert_string_to_number(data_line[i]))
                    end
                end
            end
            # end of parsing in case the last line is not empty
            if strp_line != "" && is_parsing && is_in_loop
                # skipping over conversion to np array here
                this_data_dict[block_name] = temp_dict
            end
            close(f)
            this_star_type = _determine_star_type(this_block_name_list)
            new(this_data_dict, this_block_name_list, this_key_name_list, this_star_type, this_file)
        elseif typeof(args[1]) <: Dict
            this_block_name_list = collect(keys(args[1]))
            new(args[1],
                this_block_name_list,
                [collect(keys(v)) for v in values(args[1])],
                _determine_star_type(this_block_name_list),
                this_file)
        else
            error("[RelionStar] Error: argument must be String or Dict")
        end
    end
end


function _determine_star_type(block_name_list)
    if "data_" in block_name_list
        return "data"
    elseif "data_images" in block_name_list
        return "image"
    elseif "data_model_general" in block_name_list
        return "model"
    elseif "data_general" in block_name_list
        return "final_fsc"
    elseif "data_optimiser_general" in self.block_name_list
        return "optimiser"
    elseif "data_sampling_general" in self.block_name_list
        return "sampling"
    else
        return "other"
    end
end

"""
    find(rStar, search_string)

Return a tuple containing block and key names from `rStar` that
contain `search_string`.
"""
function find(rStar, search_string)
    block_string = ""
    pos = findfirst(isequal('|'), search_string)
    if !(pos === nothing)
        block_string = strip(search_string[1:pos - 1])
        search_string = strip(search_string[pos + 1:length(search_string)])
    end
    found = []
    num_hits = 0
    for (i, block_name) in enumerate(rStar.block_name_list)
        if occursin(search_string, block_name)
            num_hits += 1
            push!(found, block_name)
            push!(found, "")
            continue
        end
        for j in rStar.key_name_list[i]
            if occursin(search_string, j)
                num_hits += 1
                push!(found, block_name)
                push!(found, j)
            end
        end
    end
    if num_hits == 0
        println("[RelionStar] No hit found!")
        return (nothing, nothing)
    elseif num_hits > 1
        if !isempty(block_string)
            temp_list = []
            for i in 1:2:length(found) # step goes in the middle
                if block_string in found[i]
                    push!(temp_list, found[i:i+1])
                end
            end
            found = temp_list
        end
        if length(found) > 2
            println("[RelionStar] Multiple hits with the search string")
            return tuple(found...)
        elseif length(found) == 2
            return tuple(found...)
        else
            println("[RelionStar] No hit found!")
            return (nothing, nothing)
        end
    else
        return tuple(found...)
    end
end

"""
    write_star(rStar, file_name)

Write a star file at the specified location with the data_dict of `rStar`,
a RelionStar object.
"""
function write_star(rStar::RelionStar, file_name)
    file = open(file_name, "w")
    for block_name in rStar.block_name_list
        write(file, block_name * "\n\n")
        in_loop = false
        for key_name in rStar.data_dict[block_name].keys
            if typeof(rStar.data_dict[block_name][key_name]) <: Array
                if in_loop
                    write(file, key_name * "\n")
                else
                    write(file, "loop_\n" * key_name * "\n")
                    in_loop = true
                end
            else
                write(file, rpad(key_name, 40) * " " * lpad(rStar.data_dict[block_name][key_name], 12) * "\n")
            end
        end
        if in_loop
            for i in 1:length(rStar.data_dict[block_name][rStar.data_dict[block_name].keys[1]])
                for key_name in rStar.data_dict[block_name].keys
                    write(file, lpad(rStar.data_dict[block_name][key_name][i], 12) * " ")
                end
                write(file, "\n")
            end
        end
        write(file, "\n\n")
    end
    close(file)
end

function _convert_string_to_number(s)
    try
        return parse(Int64, s)
    catch
        try
            return parse(Float64, s)
        catch
            return s
        end
    end
end

function read_relion_parameter(star, param, allowNothing=false)
    block, key = find(star, param)
    if block === nothing || key === nothing
        if allowNothing
            return nothing
        else
            error("[read_relion_parameter] Error: necessary parameter \"$param\" not found")
        end
    end
    return star.data_dict[block][key]
end

function read_relion_transformation(star)
    return tuple(map(x -> read_relion_parameter(star, x),
                     ["data | AngleRot",
                      "data | AngleTilt",
                      "data | AnglePsi",
                      "data | OriginX",
                      "data | OriginY"])...)
end

function read_relion_ctf(star)
    return (read_relion_parameter(star, "data | DefocusU"),
            read_relion_parameter(star, "data | DefocusV"),
            read_relion_parameter(star, "data | DefocusAngle"),
            read_relion_parameter(star, "data | SphericalAberration", true),
            read_relion_parameter(star, "data | Voltage",             true),
            read_relion_parameter(star, "data | AmplitudeContrast",   true))
end

function uniform_random(min, max, entries)
    return [rand() * (max - min) + min for _ in 1:entries]
end

function randomized_star(entries)
    data_images = Dict()
    data_images["_rlnAnglePsi #1"]     = uniform_random(0, 360, entries)
    data_images["_rlnAngleTilt #2"]    = uniform_random(0, 360, entries)
    data_images["_rlnAngleRot #3"]     = uniform_random(0, 360, entries)
    data_images["_rlnOriginX #4"]      = uniform_random(-5, 5, entries)
    data_images["_rlnOriginY #5"]      = uniform_random(-5, 5, entries)
    data_images["_rlnDefocusU #6"]     = uniform_random(5000, 35000, entries)
    data_images["_rlnDefocusV #7"]     = data_images["_rlnDefocusU #6"]
    data_images["_rlnDefocusAngle #8"] = zeros(entries)
    return RelionStar(Dict("data_images" => data_images))
end
