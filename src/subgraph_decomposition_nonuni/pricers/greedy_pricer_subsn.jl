using Random
using Graphs, MetaGraphsNext
using StatsBase



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

function solve_greedy_pricer_subsn(v_subgraph, s_subgraph, sub_instance, original_instance, vn_decompo, dual_costs, additional_costs_routing; nb_iterations=500)



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


            # Optimized edge retrieval (avoiding looking it up step-by-step if possible)
            # If your 'get_edge' is slow, just use Graphs.edge(s_graph, src, dst) or reuse edges_of_path
            path_edges = [get_edge(s_network_dir, src(e), dst(e)) for e in edges_of_path]
            edge_routing[v_edge] = Path(s_src, s_dst, path_edges, cost_of_routing_current_edge)
            
        end

        return edge_routing, overall_routing_costs

    end




    function complete_partial_placement(partial_placement)


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

            # Choose some nodes
            number_s_nodes = length(possible_s_nodes)
            some_s_nodes = sample(possible_s_nodes, number_s_nodes; replace=false)


            # Rank them according distance to already placed nodes and capacity
            distances = [ sum(shortest_paths.dists[s_src, s_node] for s_src in placement_neighbors) for s_node in some_s_nodes]
            distances_norm = (distances .- minimum(distances)) ./ (maximum(distances) - minimum(distances) + 1e-9)

            capacities = [ s_node_scores[s_node] for s_node in some_s_nodes]
            capacities_norm = (capacities .- minimum(capacities)) ./ (maximum(capacities) - minimum(capacities) + 1e-9)

            costs = [ node_costs[s_node] + additional_costs_dual[v_node][s_node] for s_node in some_s_nodes]
            costs_norm = (costs .- minimum(costs)) ./ (maximum(costs) - minimum(costs) + 1e-9)

            cut_routing = [ additional_costs_routing[v_node][s_node] for s_node in some_s_nodes]
            cut_norm = (cut_routing .- minimum(cut_routing)) ./ (maximum(cut_routing) - minimum(cut_routing) + 1e-9)




            # Here is the changee
            final_scores = [ node_capacities[some_s_nodes[idx_node]] >= node_demands[v_node] ? 3* distances_norm[idx_node] + 2* cut_norm[idx_node] : 10e9
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
            placement_cost += node_demands[v_node] * node_costs[placement[v_node]] + additional_costs_routing[v_node][placement[v_node]] + additional_costs_dual[v_node][placement[v_node]]
        end

        return placement, placement_cost


    end



    v_network = sub_instance.v_network
    s_network = sub_instance.s_network
    s_network_dir = sub_instance.s_network_dir

    original_s_network = original_instance.s_network
    original_s_network_dir = original_instance.s_network_dir

    bad_result = (sub_mapping=nothing,
        real_cost =Inf,
        reduced_cost=Inf
    )




        
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
    s_node_scores = [ node_capacities[s_node]  * 
                        sum(capacities_edges[s_node,s_neighbor] for s_neighbor in neighbors(s_network, s_node)) 
                    for s_node in vertices(s_network)
    ]




    # shortest path of the substrate network --- WITH DUAL COSTS!
    distmx = zeros(nv(s_network), nv(s_network))
    for s_edge in edges(s_network)
        original_src = s_subgraph.nodes_of_main_graph[src(s_edge)]
        original_dst = s_subgraph.nodes_of_main_graph[dst(s_edge)]
        original_edge = get_edge(original_s_network, original_src, original_dst)
        distmx[src(s_edge), dst(s_edge)] = get_attribute_edge(s_network, s_edge, :cost) - dual_costs.capacity_s_edge[original_edge]
        distmx[dst(s_edge), src(s_edge)] = get_attribute_edge(s_network, s_edge, :cost) - dual_costs.capacity_s_edge[original_edge]
    end
    shortest_paths = floyd_warshall_shortest_paths(s_network_dir, distmx)

    # penalty = dual costs for flow conservation and flow departure constraints
    additional_costs_dual = []
    for v_node in vertices(v_subgraph.graph)

        #current_additional_costs = copy(additional_costs_routing[v_node])
        current_additional_costs = zeros(length(s_network))

        original_v_node = v_subgraph.nodes_of_main_graph[v_node]
        for cut_edge in vn_decompo.v_edges_master

            if original_v_node == src(cut_edge)
                for s_node in vertices(s_network)
                    original_node = s_subgraph.nodes_of_main_graph[s_node]
                    current_additional_costs[s_node] -= dual_costs.flow_conservation[cut_edge][original_node]
                    current_additional_costs[s_node] -= dual_costs.departure[cut_edge][original_node]
                end
            end
            
            if original_v_node == dst(cut_edge)
                for s_node in vertices(s_network)
                    original_node = s_subgraph.nodes_of_main_graph[s_node]
                    current_additional_costs[s_node] += dual_costs.flow_conservation[cut_edge][original_node]
                end
            end

        end
        push!(additional_costs_dual, current_additional_costs)
    end
    
    

    # Getting the most central virtual node
    most_central_v_node = argmin(closeness_centrality(v_network))

    # Loops related things
    best_cost = 10e9
    best_mapping = nothing
    

    # Construct initial mapping
    # We iterate, for the case where it's hard to find a feasible one

    iter = 0
    centrality_nodes = closeness_centrality(s_network)

    while iter < nb_iterations

        # Construct initial mapping
        init_s_nodes_scores = [ node_capacities[s_node] >= node_demands[most_central_v_node] ? centrality_nodes[s_node] + 1. * rand() : 10e9
                                    for s_node in vertices(s_network) ]
        
        s_node_start = argmin(init_s_nodes_scores)
        if minimum(init_s_nodes_scores) > 10e8
            return bad_result
        end
        placement = zeros(Int32, nv(v_network))
        placement[most_central_v_node] = s_node_start
    
        placement, placement_cost = complete_partial_placement(placement) 
        if placement_cost < Inf
            routing, routing_cost = shortest_path_routing(placement)
        else
            routing_cost = Inf
        end        
        total_cost = placement_cost + routing_cost
        
        if total_cost < best_cost 
            best_mapping = Mapping(v_network, s_network, placement, routing)
            best_cost = total_cost
        end

        iter += 1
    end


    
    if best_cost>10e7
        return bad_result
    end


    # Get the true reduced cost
    reduced_cost = best_cost - dual_costs.convexity[v_subgraph]
    for v_node in vertices(v_network)
        reduced_cost -= additional_costs_routing[v_node][best_mapping.node_placement[v_node]]
    end
    # Get original mapping

    real_cost = 0
    real_node_placement = []
    for v_node in vertices(v_network)
        original_s_node = s_subgraph.nodes_of_main_graph[best_mapping.node_placement[v_node]]
        append!(real_node_placement, original_s_node)
        real_cost += original_s_network[original_s_node][:cost] * node_demands[v_node]
    end
    
    real_edge_routing = Dict()
    for v_edge in edges(v_subgraph.graph)
        used_edges = []
        for s_edge in best_mapping.edge_routing[v_edge].edges
            real_s_edge = get_edge(original_s_network_dir, s_subgraph.nodes_of_main_graph[src(s_edge)], s_subgraph.nodes_of_main_graph[dst(s_edge)])
            push!(used_edges, real_s_edge)
            real_cost += original_s_network_dir[src(real_s_edge), dst(real_s_edge)][:cost] * edge_demands[src(v_edge), dst(v_edge)]
        end
        real_edge_routing[v_edge] = order_path(original_s_network_dir, used_edges, real_node_placement[src(v_edge)], real_node_placement[dst(v_edge)]) 
    end


    real_mapping = Mapping(v_subgraph.graph, original_s_network_dir, real_node_placement, real_edge_routing)

    return (sub_mapping=real_mapping,
            real_cost=real_cost,
            reduced_cost=reduced_cost
    )

end

