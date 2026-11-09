using JSON
using DataStructures
using Graphs, MetaGraphsNext
includet("../../utils/graph.jl")



function write_network_to_json(mg, is_directed=true)
    formatted_json_str = "{\n"
    
    # about the graph:
    for k in keys(mg[])
        formatted_json_str *= "\t\"" * string(k) * "\": "
        formatted_json_str *= "\"" * string(mg[][k]) * "\""

        formatted_json_str *= ",\n"
    end
    if is_directed
        formatted_json_str *= "\t\"directed\": true,\n"
    else
        formatted_json_str *= "\t\"directed\": false,\n"
    end

    # node:
    formatted_json_str *= "\t\"nodes\": [\n"
    for node in vertices(mg)
        formatted_json_str = formatted_json_str * "\t\t{\"id\": " * string(node) * ","
        for k in keys(mg[node])
            formatted_json_str *= " \"" * string(k) * "\": " * string(mg[node][k]) * ","
        end
        formatted_json_str = formatted_json_str[1:end-1] 
        formatted_json_str *= "},\n"
    end
    formatted_json_str = formatted_json_str[1:end-2] 
    formatted_json_str *= "\n\t],"

    # edges
    formatted_json_str *= "\n\t\"edges\": [\n"
    for edge in edges(mg)
        formatted_json_str = formatted_json_str * "\t\t{\"source\": " * string(src(edge)) * ", \"target\": " * string(dst(edge)) * ", "
        for k in keys(mg[src(edge), dst(edge)])
            formatted_json_str *= " \"" * string(k) * "\": " * string(mg[src(edge), dst(edge)][k]) * ","
        end
        formatted_json_str = formatted_json_str[1:end-1] 
        formatted_json_str *= "},\n"
    end
    formatted_json_str = formatted_json_str[1:end-2] 
    formatted_json_str *= "\n\t]"

    formatted_json_str *= "\n}"

    # Write the formatted JSON string to a file
    filename = mg[][:name] * ".json"
    open(filename, "w") do f
        write(f, formatted_json_str)
    end
end


function read_gml_file(file_path::String)
    open(file_path, "r") do file
        return read(file, String)
    end
end



function get_graph_from_gml(path)
    gml_content = read_gml_file(path)

    g = Graph()

    gml_splitted = split(gml_content, '\n')
    i_line = 1
    name = "unknown-network"

    while i_line < length(gml_splitted)
        line = strip(gml_splitted[i_line])

        if startswith(line, "Network ") 
            key, name = split(line, " ")
            name = strip(name, ['"'])
        end


        if startswith(line, "node") 
            add_vertex!(g)
        end

        if startswith(line, "edge")
            i_line += 1
            line = strip(gml_splitted[i_line])
            key, src = split(line, " ")
            i_line += 1
            line = strip(gml_splitted[i_line])
            key, dst = split(line, " ")
            add_edge!(g, parse(Int, src)+1, parse(Int, dst)+1)
        end

        i_line = i_line+1
    end

    return g, name
end



function read_graph(filename::String)
    println("Ok so")
    json_graph = JSON.parsefile(filename)
    if json_graph["type"] == "virtual"
        g = read_virtual(json_graph)
        return g, 1
    end
    if json_graph["type"] == "substrate"
        g = read_substrate(json_graph)
        return g, 2
    end
end

function read_substrate(json_graph)
    # testing whether it is directed or undirected. By default: undirected.
    if "directed" in keys(json_graph)
        if json_graph["directed"]
            g = MetaGraph(
                DiGraph(),
                Int,
                Dict,
                Dict,
                Dict(:name => json_graph["name"], :type => "substrate", :directed => json_graph["directed"])
            )
        else
            g = MetaGraph(
                Graph(),
                Int,
                Dict,
                Dict,
                Dict(:name => json_graph["name"], :type => "substrate", :directed => json_graph["directed"])
            )
        end
    else
        g = MetaGraph(
            DiGraph(),
            Int,
            Dict,
            Dict,
            Dict(:name => json_graph["name"], :type => "substrate", :directed => true)
        )
    end


    for node in json_graph["nodes"]
        add_vertex!(g, node["id"], Dict{Any, Any}(:cap=> node["cap"], :cost => node["cost"]))
    end

    for edge in json_graph["edges"]
        add_edge!(g, edge["source"], edge["target"], Dict{Any, Any}(:cap => edge["cap"], :cost => edge["cost"]))
    end

    return g
end

function read_virtual(json_graph)


    # testing whether it is directed or undirected. By default: directed.
    if "directed" in keys(json_graph)
        if json_graph["directed"]
            g = MetaGraph(
                DiGraph(),
                Int,
                Dict,
                Dict,
                Dict(:name => json_graph["name"], :type => "virtual", :directed => json_graph["directed"])
            )
        else
            g = MetaGraph(
                Graph(),
                Int,
                Dict,
                Dict,
                Dict(:name => json_graph["name"], :type => "virtual", :directed => json_graph["directed"])
            )
        end
    else
        g = MetaGraph(
            DiGraph(),
            Int,
            Dict,
            Dict,
            Dict(:name => json_graph["name"], :type => "virtual", :directed => true)
        )
    end

    for node in json_graph["nodes"]
        add_vertex!(g, node["id"], Dict{Any, Any}(:dem => node["dem"]))
    end

    for edge in json_graph["edges"]
        add_edge!(g, edge["source"], edge["target"], Dict{Any, Any}(:dem => edge["dem"]))
    end

    return g
end

