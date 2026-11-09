using Graphs, MetaGraphsNext

# Here, I'm doing some stuff extending the graph package. Also, some function related to paths.

#### Paths ####

struct Path
    src
    dst
    edges
    cost
end

function Path()
    return Path(0, 0, [], 10000)
end




function Base.show(
    io::IO, meta_graph::MetaGraph{<:Any,BaseGraph,Label,VertexData,EdgeData}
) where {BaseGraph,Label,VertexData,EdgeData}
    print(
        io,
        "MetaGraph_" * meta_graph[][:name],
    )
    return nothing
end






function get_attribute_node(mg, element, attribute)
    if element in vertices(mg)
        return mg[element][attribute]
    else
        println("Element not found!")
    end
end

function get_attribute_edge(mg, element, attribute)
    if element in edges(mg)
        return mg[src(element), dst(element)][attribute]
    else
        println("Element not found!")
    end
end

function get_attribute_graph(mg, attribute)
    return mg[][attribute]
end

function set_attribute_edge(mg, element, attribute, new_thing)
    if element in edges(mg)
        mg[src(element), dst(element)][attribute] = new_thing
    else
        println("Element not found!")
    end
end

function set_attribute_node(mg, element, attribute, new_thing)
    if element in vertices(mg)
        mg[element][attribute] = new_thing
    else
        println("Element not found!")
    end
end



function get_reverse_edge(g, edge)
    if get_edge(g, dst(edge), src(edge)) ∈ edges(g)
        return get_edge(g, dst(edge), src(edge))
    else
        println("The reverse edge of $edge in $g has not been found...")
    end
end


# Only for undirected graphs !
# !!!!!!!!!!! SELECTOR NEEDS TO BE INT64 TYPE 
function my_induced_subgraph(meta_graph, selector, name)
    induced_graph, osef = induced_subgraph(meta_graph.graph, selector)

    induced_metagraph = MetaGraph(
        Graph(),
        Int,
        Dict,
        Dict,
        Dict(:name => name, :type => meta_graph[][:type], :directed => false)
    )

    for node in vertices(induced_graph)
        add_vertex!(induced_metagraph, node, meta_graph[selector[node]])
    end

    for edge in edges(induced_graph)
        add_edge!(induced_metagraph, src(edge), dst(edge), meta_graph[selector[src(edge)], selector[dst(edge)]])
    end

    return induced_metagraph
end

# Il y a deja outneighbours à utiliser :!!!
function get_out_edges(g, node)
    out_edges = []
    for edge in edges(g)
        if src(edge) == node
            push!(out_edges, edge)
        end  
    end
    return out_edges
end

function get_in_edges(g, node)
    in_edges = []
    for edge in edges(g)
        if dst(edge) == node
            push!(in_edges, edge)
        end  
    end
    return in_edges
end

function get_edge(g, i, j)
    for edge in edges(g)
        if (src(edge) == i) && (dst(edge) == j)
            return edge
        end
    end
    for edge in edges(g)
        if (dst(edge) == i) && (src(edge) == j)
            return edge
        end
    end
    println("Edge not found: " * string(i) * "=>" * string(j))
    print_graph(g)
    
    throw("Edge not found")
end 





# Cette fonction pourrait poser probleme si on met des chemins bizaroides. Pas la plus secure.
function order_path(s_network, used_edges, start, target)
    unordered_path = copy(used_edges)
    ordered_path = []
    path_cost = 0
    current_node = start
    while current_node != target
        has_found_next_edge = false
        for edge in unordered_path
            if src(edge) == current_node
                path_cost += s_network[src(edge), dst(edge)][:cost]
                push!(ordered_path, edge)
                current_node = dst(edge)
                filter!(!=(edge), unordered_path)
                has_found_next_edge = true
                break
            elseif dst(edge) == current_node # FOR UNDIRECTED...
                path_cost += s_network[dst(edge), src(edge)][:cost]
                push!(ordered_path, edge)
                current_node = src(edge)
                filter!(!=(edge), unordered_path)
                has_found_next_edge = true
                break
            end
            
        end
        if has_found_next_edge == false
            error("Failed to find next edge in path.")
            return([])
        end
    end
    return Path(start, target, ordered_path, path_cost)
end

#=

function get_shortest_paths(g, k)
    # Initialize an empty matrix with zeros
    n_vertices = size(vertices(g), 1)
    distmx = zeros(Int, size(vertices(g), 1), n_vertices)

    # Populate the matrix using a nested for loop
    for edge in edges(g)
        if g[][:type] == "substrate"
            distmx[src(edge), dst(edge)] = g[src(edge), dst(edge)][:cost]
        else
            distmx[src(edge), dst(edge)] = 1
        end
    end
    yen_paths = [[yen_k_shortest_paths(g, i, j, distmx, k) 
        for i in range(1, n_vertices)] 
        for j in range(1, n_vertices)
    ]

    shortest_paths_dic = Dict()
    for i in range(1, n_vertices)
        for j in range(1, n_vertices)
            paths_bad = yen_paths[i][j].paths
            paths_good = []
            for path in paths_bad
                push!(paths_good, good_path_of(g, path)) 
            end
            shortest_paths_dic[(i,j)] = paths_good
        end
    end

    return shortest_paths_dic
end


function good_path_of(g, nodes_of_path)
    # edges are supposed to be ordered
    start_node = first(nodes_of_path)
    terminus_node = last(nodes_of_path)
    path = []
    cost = 0

    for i in range(2, size(nodes_of_path, 1))
        push!(path, find_edge(g, nodes_of_path[i-1], nodes_of_path[i]))
        if g[][:type] == "substrate"
            cost += g[src(last(path)), dst(last(path))][:cost]
        else
            cost += 1
        end

        
    end

    return Path(start_node, terminus_node, path, cost)
end

# to delete ?
function path_in_paths(new_path, paths)
    for path in paths
        if length(path.edges) == length(new_path.edges)
            is_same = true
            for i_edge in 1:length(path.edges)
                if (path.edges[i_edge].src != new_path.edges[i_edge].src) || (path.edges[i_edge].dst) != (new_path.edges[i_edge].dst)
                    is_same = false
                end
            end
            if is_same
                #println("Path already in paths : " * string(new_path) * " == " * string(new_path))
                return true
            end
        end
    end
    return false
end


function find_edge(g, start_node, end_node)
    # Find the edge object
    edge = nothing
    for e in edges(g)
        if src(e) == start_node && dst(e) == end_node
            edge = e
            break
        end
    end

    # Check if the edge exists
    if isnothing(edge)
        println("Edge not found!")
    end
    return edge
end
=#

function print_graph(mg)
    println("Metagraph ")
    print("Nodes:")
    for node in vertices(mg)
        print("\nNode " * string(node) * " with ")
        for k in keys(mg[node])
            print(string(k) * ": " * string(mg[node][k]))
        end
    end
    print("\nEdges:")
    for edge in edges(mg)
        print("\nEdge " * string(edge))
        for k in keys(mg[src(edge), dst(edge)])
            print(" with $k : $(mg[src(edge), dst(edge)][k])" )
        end
    end
end



