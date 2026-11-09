using Random
using Graphs, MetaGraphsNext
using StatsBase

includet("../utils/import_utils.jl")

# --- Type-stable Custom Weight Type ---
struct DynamicWeightMatrix{M, C, T} <: AbstractMatrix{T}
    base_dist::M
    capacities::C
    demand::T
end

# Outer constructor to infer type easily
function DynamicWeightMatrix(base_dist, capacities, demand)
    return DynamicWeightMatrix{typeof(base_dist), typeof(capacities), typeof(demand)}(base_dist, capacities, demand)
end

# Implement the minimum required interface for Graphs.jl distance matrix
Base.size(d::DynamicWeightMatrix) = size(d.base_dist)
@inline function Base.getindex(d::DynamicWeightMatrix, u::Int, v::Int)
    # If the substrate link cannot support the current virtual demand, it's blocked (Inf)
    if d.capacities[u, v] < d.demand
        return Inf
    else
        return d.base_dist[u, v]
    end
end

function solve_local_search(instance, initial_mapping; nb_local_search = 500, time_max = 10, nb_nodes_removed = 5)




    function shortest_path_routing(v_node_placement)

        capacities_edges_copy = copy(capacities_edges)

        edge_routing = Dict{Edge, Path}()
        overall_routing_costs = 0

        # In which order to route the edges?
        virtual_edges = collect(edges(v_network))
        shuffle!(virtual_edges)

        for v_edge in virtual_edges

            s_src = v_node_placement[src(v_edge)]
            s_dst = v_node_placement[dst(v_edge)]
            demand_curr_edge = edge_demands[src(v_edge), dst(v_edge)]

            # --- THE FAST TRICK: Dynamic Weight View ---
            # We wrap our distances. If an edge lacks capacity for THIS specific demand, 
            # it registers as Inf cost to A*, bypassing the need to modify the graph.
            current_weights = DynamicWeightMatrix(cost_edges, capacities_edges_copy, demand_curr_edge)

            # Run A* using our smart weight view
            edges_of_path = a_star(s_network_dir, s_src, s_dst, current_weights)

            if isempty(edges_of_path)
                return Dict{Edge, Path}(), Inf  # Embedding failed due to capacity constraints
            end

            # Calculate path cost and collect nodes sequentially
            # Pre-allocating nodes_of_path size for speed
            nodes_of_path = Vector{Int}(undef, length(edges_of_path) + 1)
            nodes_of_path[1] = src(edges_of_path[1])
            cost_of_routing_current_edge = 0
            for (idx, edge) in enumerate(edges_of_path)
                u, v = src(edge), dst(edge)
                nodes_of_path[idx + 1] = v
                cost_of_routing_current_edge += cost_edges[u, v] 
                
                # Removing capacities (undirected / bidirected view)
                capacities_edges_copy[u, v] -= demand_curr_edge
                capacities_edges_copy[v, u] -= demand_curr_edge
            end

            overall_routing_costs += cost_of_routing_current_edge * demand_curr_edge

            # Only build the detailed path object if explicitly requested
            # Optimized edge retrieval (avoiding looking it up step-by-step if possible)
            # If your 'get_edge' is slow, just use Graphs.edge(s_graph, src, dst) or reuse edges_of_path
            #path_edges = [get_edge(s_network_dir, src(e), dst(e)) for e in edges_of_path]
            edge_routing[v_edge] = Path(s_src, s_dst, edges_of_path, cost_of_routing_current_edge)
        end

        return edge_routing, overall_routing_costs

    end



    function complete_partial_placement(partial_placement; nb_nodes_to_try=10)

        placement = copy(partial_placement)

        #println("Partial placement : $partial_placement")
        already_placed_v_nodes = filter(v_node -> placement[v_node] != 0, vertices(v_network))


        next_v_nodes = setdiff(
            union([neighbors(v_network, v) for v in already_placed_v_nodes]...),
            already_placed_v_nodes,
        )

        possible_s_nodes = filter(s_node -> s_node ∉ placement, vertices(s_network))

        while !isempty(next_v_nodes)

            # Take a node of the list
            shuffle!(next_v_nodes)
            v_node = popfirst!(next_v_nodes)

            # Get neighbors already placed
            placement_neighbors = [placement[s_neigh] for s_neigh in filter(v_neighbor -> placement[v_neighbor] != 0, neighbors(v_network, v_node))]

            # Choose some random nodes
            number_s_nodes = min(length(possible_s_nodes), nb_nodes_to_try)
            some_s_nodes = sample(possible_s_nodes, number_s_nodes; replace=false)


            # NODE RANKING AND SELECTION
            distances = [ sum(shortest_paths.dists[s_src, s_node] for s_src in placement_neighbors) for s_node in some_s_nodes]
            distances_norm = (distances .- minimum(distances)) ./ (maximum(distances) - minimum(distances) + 1e-9)

            #capacities = [ s_node_scores[s_node] for s_node in some_s_nodes]
            #capacities_norm = (capacities .- minimum(capacities)) ./ (maximum(capacities) - minimum(capacities) + 1e-9)

            #costs = [ node_costs[s_node] for s_node in some_s_nodes]
            #costs_norm = (costs .- minimum(costs)) ./ (maximum(costs) - minimum(costs) + 1e-9)

            final_scores = [ node_capacities[some_s_nodes[idx_node]] >= node_demands[v_node] ? distances_norm[idx_node] : 10e9
                                for idx_node in 1:length(some_s_nodes)]
            selected_idx = argmin(final_scores)
            s_node_selected = some_s_nodes[selected_idx]

            if minimum(final_scores) > 10e8
                return [], Inf
            end

            # Finish the work
            placement[v_node] = s_node_selected
            push!(already_placed_v_nodes, v_node)
            next_v_nodes = next_v_nodes ∪ filter(v_neighbor->v_neighbor ∉ already_placed_v_nodes, neighbors(v_network, v_node) )
            possible_s_nodes = filter(!=(s_node_selected), possible_s_nodes)
        end

        placement_cost = 0
        for v_node in vertices(v_network)
            placement_cost += node_costs[placement[v_node]] * node_demands[v_node]
        end

        return placement, placement_cost

    end




    time_beginning = time()

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir


    
    #---- Usefull things
    node_demands = [get_attribute_node(v_network, v_node, :dem) for v_node in vertices(v_network)]
    edge_demands = zeros(Int, nv(s_network), nv(s_network))
    for v_edge in edges(v_network)
        edge_demands[src(v_edge),dst(v_edge)] = get_attribute_edge(v_network, v_edge, :dem)
    end

    node_capacities = [get_attribute_node(s_network, s_node, :cap) for s_node in vertices(s_network)]
    node_costs = [get_attribute_node(s_network, s_node, :cost) for s_node in vertices(s_network)]
    capacities_edges = zeros(Int, nv(s_network), nv(s_network))
    for s_edge in edges(s_network)
        cap = s_network[src(s_edge), dst(s_edge)][:cap]
        capacities_edges[src(s_edge), dst(s_edge)] = cap
        capacities_edges[dst(s_edge), src(s_edge)] = cap
    end
    cost_edges = zeros(Int, nv(s_network), nv(s_network))
    for s_edge in edges(s_network)
        cost = s_network[src(s_edge), dst(s_edge)][:cost]
        cost_edges[src(s_edge), dst(s_edge)] = cost
        cost_edges[dst(s_edge), src(s_edge)] = cost
    end
    
    # Greedy score based on capacities for nodes
    s_node_scores = [ sum(capacities_edges[s_node,s_neighbor] + node_capacities[s_neighbor] for s_neighbor in neighbors(s_network, s_node) ) 
                    for s_node in vertices(s_network)]


    shortest_paths = floyd_warshall_shortest_paths(s_network_dir, cost_edges)
    


    # Let's go
    #println("Graph r : $(instance.v_network.edges)")
    #println("initial mapping: $(initial_mapping.edge_routing)")
    if initial_mapping == nothing
        return (mapping=nothing, mapping_cost=Inf)
    end
    best_cost = get_cost_placement(initial_mapping) + get_cost_routing(initial_mapping)
    best_placement = initial_mapping.node_placement

    println("Before Local Search, we have solution of $best_cost !")
    
    step=0
    # Do a number of improving iterations
    while step < nb_local_search && time()-time_beginning < time_max
        step+=1
        placement = copy(best_placement)
        
        #= delete some random star in the virtual network
        some_v_node = rand(1:nv(v_network))
        v_nodes_deleted = [some_v_node]
        for v_neighbor in neighbors(v_network, some_v_node) 
            push!(v_nodes_deleted, v_neighbor)
        end
        =#


        
        # 1. Pick a random starting node safely
        some_v_node = rand(1:nv(v_network))
        v_nodes_deleted = [some_v_node]

        # Maintain a tracking set of candidate neighbors we can expand into
        boundary_candidates = Set{Int}(neighbors(v_network, some_v_node))

        while length(v_nodes_deleted) < nb_nodes_removed

            if isempty(boundary_candidates)
                break
            else
                # Standard case: pick a random node sitting on our cluster's boundary
                next_node = rand(boundary_candidates)
            end
            
            push!(v_nodes_deleted, next_node)
            delete!(boundary_candidates, next_node)
            
            for neighbor in neighbors(v_network, next_node)
                if neighbor ∉ v_nodes_deleted
                    push!(boundary_candidates, neighbor)
                end
            end
        end

        for v_node in v_nodes_deleted
            placement[v_node] = 0
        end



        # reconstruct
        placement, placement_cost = complete_partial_placement(placement; nb_nodes_to_try=nv(s_network)) 
        routing, routing_cost = shortest_path_routing(placement)

        current_cost = placement_cost + routing_cost
        if current_cost < best_cost 
            best_cost = current_cost
            best_placement = placement
            println("New best! $current_cost, at iter $step")
        end
    end


    if best_cost>10e9
        return (mapping=nothing,
                    mapping_cost=Inf
        )
    end


    routing, routing_cost= shortest_path_routing(best_placement)
    final_mapping = Mapping(v_network, s_network, best_placement, routing)

    println("Find the solution of $best_cost in $(time()-time_beginning) and $step iters")
    return (mapping = final_mapping,
            mapping_cost = best_cost)
end




function solve_local_search_routing(instance, initial_mapping; nb_local_search = 500, time_max = 10, nb_nodes_removed = 5)




    function shortest_path_routing(v_node_placement, edge_routing)

        capacities_edges_copy = copy(capacities_edges)

        overall_routing_costs = 0
        for v_edge in keys(edge_routing)
            overall_routing_costs += edge_routing[v_edge].cost * edge_demands[src(v_edge), dst(v_edge)]
        end

        # In which order to route the edges?
        virtual_edges = collect(edges(v_network))
        filter!(e -> e ∉ keys(edge_routing), virtual_edges)

        shuffle!(virtual_edges)

        for v_edge in virtual_edges

            s_src = v_node_placement[src(v_edge)]
            s_dst = v_node_placement[dst(v_edge)]
            demand_curr_edge = edge_demands[src(v_edge), dst(v_edge)]

            # --- THE FAST TRICK: Dynamic Weight View ---
            # We wrap our distances. If an edge lacks capacity for THIS specific demand, 
            # it registers as Inf cost to A*, bypassing the need to modify the graph.
            current_weights = DynamicWeightMatrix(cost_edges, capacities_edges_copy, demand_curr_edge)

            # Run A* using our smart weight view
            edges_of_path = a_star(s_network_dir, s_src, s_dst, current_weights)

            if isempty(edges_of_path)
                return Dict{Edge, Path}(), Inf  # Embedding failed due to capacity constraints
            end

            # Calculate path cost and collect nodes sequentially
            # Pre-allocating nodes_of_path size for speed
            nodes_of_path = Vector{Int}(undef, length(edges_of_path) + 1)
            nodes_of_path[1] = src(edges_of_path[1])
            cost_of_routing_current_edge = 0
            for (idx, edge) in enumerate(edges_of_path)
                u, v = src(edge), dst(edge)
                nodes_of_path[idx + 1] = v
                cost_of_routing_current_edge += cost_edges[u, v] 
                
                # Removing capacities (undirected / bidirected view)
                capacities_edges_copy[u, v] -= demand_curr_edge
                capacities_edges_copy[v, u] -= demand_curr_edge
            end

            overall_routing_costs += cost_of_routing_current_edge * demand_curr_edge

            # Only build the detailed path object if explicitly requested
            # Optimized edge retrieval (avoiding looking it up step-by-step if possible)
            # If your 'get_edge' is slow, just use Graphs.edge(s_graph, src, dst) or reuse edges_of_path
            #path_edges = [get_edge(s_network_dir, src(e), dst(e)) for e in edges_of_path]
            edge_routing[v_edge] = Path(s_src, s_dst, edges_of_path, cost_of_routing_current_edge)
        end

        return edge_routing, overall_routing_costs

    end



    function complete_partial_placement(partial_placement; nb_nodes_to_try=10)

        placement = copy(partial_placement)

        #println("Partial placement : $partial_placement")
        already_placed_v_nodes = filter(v_node -> placement[v_node] != 0, vertices(v_network))


        next_v_nodes = setdiff(
            union([neighbors(v_network, v) for v in already_placed_v_nodes]...),
            already_placed_v_nodes,
        )

        possible_s_nodes = filter(s_node -> s_node ∉ placement, vertices(s_network))

        while !isempty(next_v_nodes)

            # Take a node of the list
            shuffle!(next_v_nodes)
            v_node = popfirst!(next_v_nodes)

            # Get neighbors already placed
            placement_neighbors = [placement[s_neigh] for s_neigh in filter(v_neighbor -> placement[v_neighbor] != 0, neighbors(v_network, v_node))]

            # Choose some random nodes
            number_s_nodes = min(length(possible_s_nodes), nb_nodes_to_try)
            some_s_nodes = sample(possible_s_nodes, number_s_nodes; replace=false)


            # NODE RANKING AND SELECTION
            distances = [ sum(shortest_paths.dists[s_src, s_node] for s_src in placement_neighbors) for s_node in some_s_nodes]
            distances_norm = (distances .- minimum(distances)) ./ (maximum(distances) - minimum(distances) + 1e-9)

            #capacities = [ s_node_scores[s_node] for s_node in some_s_nodes]
            #capacities_norm = (capacities .- minimum(capacities)) ./ (maximum(capacities) - minimum(capacities) + 1e-9)

            #costs = [ node_costs[s_node] for s_node in some_s_nodes]
            #costs_norm = (costs .- minimum(costs)) ./ (maximum(costs) - minimum(costs) + 1e-9)

            final_scores = [ node_capacities[some_s_nodes[idx_node]] >= node_demands[v_node] ? distances_norm[idx_node] : 10e9
                                for idx_node in 1:length(some_s_nodes)]
            selected_idx = argmin(final_scores)
            s_node_selected = some_s_nodes[selected_idx]

            if minimum(final_scores) > 10e8
                return [], Inf
            end

            # Finish the work
            placement[v_node] = s_node_selected
            push!(already_placed_v_nodes, v_node)
            next_v_nodes = next_v_nodes ∪ filter(v_neighbor->v_neighbor ∉ already_placed_v_nodes, neighbors(v_network, v_node) )
            possible_s_nodes = filter(!=(s_node_selected), possible_s_nodes)
        end

        placement_cost = 0
        for v_node in vertices(v_network)
            placement_cost += node_costs[placement[v_node]] * node_demands[v_node]
        end

        return placement, placement_cost

    end

    keep_routing = true


    time_beginning = time()

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir


    
    #---- Usefull things
    node_demands = [get_attribute_node(v_network, v_node, :dem) for v_node in vertices(v_network)]
    edge_demands = zeros(Int, nv(s_network), nv(s_network))
    for v_edge in edges(v_network)
        edge_demands[src(v_edge),dst(v_edge)] = get_attribute_edge(v_network, v_edge, :dem)
    end

    node_capacities = [get_attribute_node(s_network, s_node, :cap) for s_node in vertices(s_network)]
    node_costs = [get_attribute_node(s_network, s_node, :cost) for s_node in vertices(s_network)]
    capacities_edges = zeros(Int, nv(s_network), nv(s_network))
    for s_edge in edges(s_network)
        cap = s_network[src(s_edge), dst(s_edge)][:cap]
        capacities_edges[src(s_edge), dst(s_edge)] = cap
        capacities_edges[dst(s_edge), src(s_edge)] = cap
    end
    cost_edges = zeros(Int, nv(s_network), nv(s_network))
    for s_edge in edges(s_network)
        cost = s_network[src(s_edge), dst(s_edge)][:cost]
        cost_edges[src(s_edge), dst(s_edge)] = cost
        cost_edges[dst(s_edge), src(s_edge)] = cost
    end
    
    # Greedy score based on capacities for nodes
    s_node_scores = [ sum(capacities_edges[s_node,s_neighbor] + node_capacities[s_neighbor] for s_neighbor in neighbors(s_network, s_node) ) 
                    for s_node in vertices(s_network)]


    shortest_paths = floyd_warshall_shortest_paths(s_network_dir, cost_edges)
    


    # Let's go
    #println("Graph r : $(instance.v_network.edges)")
    #println("initial mapping: $(initial_mapping.edge_routing)")
    if initial_mapping == nothing
        return (mapping=nothing, mapping_cost=Inf)
    end
    best_cost = get_cost_placement(initial_mapping) + get_cost_routing(initial_mapping)
    best_mapping = initial_mapping

    println("Before Local Search, we have solution of $best_cost !")
    
    step=0
    # Do a number of improving iterations
    while step < nb_local_search && time()-time_beginning < time_max
        step+=1
        placement = copy(best_mapping.node_placement)
        
        #= delete some random star in the virtual network
        some_v_node = rand(1:nv(v_network))
        v_nodes_deleted = [some_v_node]
        for v_neighbor in neighbors(v_network, some_v_node) 
            push!(v_nodes_deleted, v_neighbor)
        end
        =#


        
        # 1. Pick a random starting node safely
        some_v_node = rand(1:nv(v_network))
        v_nodes_deleted = [some_v_node]

        # Maintain a tracking set of candidate neighbors we can expand into
        boundary_candidates = Set{Int}(neighbors(v_network, some_v_node))

        while length(v_nodes_deleted) < nb_nodes_removed

            if isempty(boundary_candidates)
                break
            else
                # Standard case: pick a random node sitting on our cluster's boundary
                next_node = rand(boundary_candidates)
            end
            
            push!(v_nodes_deleted, next_node)
            delete!(boundary_candidates, next_node)
            
            for neighbor in neighbors(v_network, next_node)
                if neighbor ∉ v_nodes_deleted
                    push!(boundary_candidates, neighbor)
                end
            end
        end

        for v_node in v_nodes_deleted
            placement[v_node] = 0
        end



        # reconstruct
        placement, placement_cost = complete_partial_placement(placement; nb_nodes_to_try=nv(s_network)) 

        routing = nothing
        routing_cost = Inf
        if placement_cost < 10e8
            routing = Dict{Edge, Path}()
            if keep_routing
                for v_edge in edges(v_network)
                    if src(v_edge) ∉ v_nodes_deleted && dst(v_edge) ∉ v_nodes_deleted
                        routing[v_edge] = best_mapping.edge_routing[v_edge]
                    end
                end
            end
            routing, routing_cost = shortest_path_routing(placement, routing)
        end

        current_cost = placement_cost + routing_cost
        if current_cost < best_cost 
            best_cost = current_cost
            best_mapping = Mapping(v_network, s_network, placement, routing)
            println("New best! $current_cost, at iter $step")
        end
    end


    if best_cost>10e9
        return (mapping=nothing,
                    mapping_cost=Inf
        )
    end



    println("Find the solution of $best_cost in $(time()-time_beginning) and $step iters")
    return (mapping = best_mapping,
            mapping_cost = best_cost)
end


