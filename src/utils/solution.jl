### MAPPING ###

struct Mapping
    virtual_network
    substrate_network
    node_placement
    edge_routing
end





function get_cost_placement(mapping)
    v_network = mapping.virtual_network
    s_network = mapping.substrate_network
    node_placement_cost = 0

    for v_node in vertices(v_network)
        node_placement_cost += s_network[mapping.node_placement[v_node]][:cost] * v_network[v_node][:dem]
    end

    return node_placement_cost
end


function get_usage_placement(mapping)

    s_node_usage = []
    for s_node in vertices(s_network)
        usage = 0
        for v_node in vertices(v_network)
            if node_placement[v_node] == s_node
                usage += 1
            end
        end
        push!(s_node_usage, usage)
    end

    return s_node_usage
end


function get_cost_routing(mapping)
    v_network = mapping.virtual_network

    edge_routing_cost = 0
    for v_edge in edges(v_network)
        path = mapping.edge_routing[v_edge]
        edge_routing_cost += path.cost * v_network[src(v_edge), dst(v_edge)][:dem]
    end
    return edge_routing_cost
end

function get_usage_routing(mapping)

    s_edge_usage = Dict()
    for s_edge in edges(s_network)
        usage = 0
        for v_edge in edges(v_network)
            if s_edge in edge_routing[v_edge].edges
                usage += v_network[src(v_edge), dst(v_edge)][:dem]
            end
        end
        s_edge_usage[s_edge] = usage
    end

    return s_edge_usage 
end




function Base.show(io::IO, mapping::Mapping)
    print(io, "Mapping of $(mapping.virtual_network) on $(mapping.substrate_network) with placement $(mapping.node_placement)")
end



function solution_checker(instance, mapping)

    v_network = instance.v_network
    solution_is_valid = true

    # Capacité noeuds
    allocated_demands = zeros(Float64, nv(instance.s_network))
    nb_virtual_node = zeros(Int32, nv(instance.s_network))
    for v_node in vertices(v_network)
        allocated_demands[mapping.node_placement[v_node]] += v_network[v_node][:dem]
        nb_virtual_node[mapping.node_placement[v_node]] += 1
    end

    for s_node in vertices(instance.s_network)
        dem_totale = allocated_demands[s_node]
        cap_max = instance.s_network[s_node][:cap]
        
        if dem_totale > cap_max
            println("Capacity of node $s_node not respected: $dem_totale used out of $cap_max capacity")
            solution_is_valid = false
        end
        if nb_virtual_node[s_node] > 1
            println("One to one for node $s_node not respected")
            solution_is_valid = false
        end
    end    
    
    # Capacité arretes
    for s_edge in edges(instance.s_network)
        dem_totale = 0
        for v_edge in edges(v_network)
            if s_edge in mapping.edge_routing[v_edge].edges
                dem_totale += v_network[src(v_edge), dst(v_edge)][:dem]
            end
            if get_reverse_edge(instance.s_network_dir, s_edge) ∈ mapping.edge_routing[v_edge].edges
                dem_totale += v_network[src(v_edge), dst(v_edge)][:dem]
            end
        end
        cap_edge = instance.s_network[src(s_edge), dst(s_edge)][:cap]
        if dem_totale > cap_edge
            println("Capacity of edge " * string(s_edge) * " not respected : " * string(dem_totale) * " used for " * string(cap_edge))
            solution_is_valid = false
        end
    end


    if solution_is_valid
        println("Solution is valid")
    else
        throw("Error : invalid solution")
    end
end
