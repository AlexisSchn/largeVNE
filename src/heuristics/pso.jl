using Random
using Graphs, MetaGraphsNext
using StatsBase
includet("../utils/import_utils.jl")







function solve_PSO_better(instance; initial_particles = [], nb_particle=25, nb_iter=50, time_max=5, print_things=true)


 



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

        already_placed_v_nodes = filter(v_node -> placement[v_node] != 0, vertices(v_network))

        #println("Partial placement : $partial_placement")
        if already_placed_v_nodes == []
            init_s_nodes_scores = [ node_capacities[s_node] >= node_demands[most_central_v_node] ? centrality_nodes[s_node] + 1. * rand() : 10e9
                                    for s_node in vertices(s_network) ]
        
            s_node_start = argmin(init_s_nodes_scores)
            if minimum(init_s_nodes_scores) > 10e8
                return (mappings=[],
                        mapping_cost=Inf
                )
            end
            placement[most_central_v_node] = s_node_start
            already_placed_v_nodes = [most_central_v_node]
        end


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

            final_scores = [ node_capacities[some_s_nodes[idx_node]] >= node_demands[v_node] ? distances_norm[idx_node] : Inf
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



    function minus(pos1, pos2)

        res=[]
        for i in 1:length(pos1)
            if pos1[i] == pos2[i]
                push!(res, 1)
            else
                push!(res, 0)
            end
        end
        return res
    end
    
    
    function plus(vel_inertia, vel_pb, vel_gb)
    
        p_inertia = 0.2
        p_attraction_personal = 0.3
        p_attraction_global = 0.5
    
        new_velocity = []
        for i in 1:length(vel_inertia)
            r = rand()
            if r < p_inertia
                push!(new_velocity, vel_inertia[i])
            elseif r < (p_inertia + p_attraction_personal)
                push!(new_velocity, vel_pb[i])
            else
                push!(new_velocity, vel_gb[i])
            end
        end
    
        return new_velocity

    end

    

    function times(position, velocity)

        new_placement = []
        placement_cost = 0
    
        for i in 1:nv(instance.v_network)
            if velocity[i] == 1
                push!(new_placement, position[i])
                placement_cost += node_costs[position[i]]
            else
                push!(new_placement, 0)
            end
        end

        final_placement, cost = complete_partial_placement(new_placement, nb_nodes_to_try=nv(s_network))
    
    
        return final_placement, cost
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


    

    # Getting the most central virtual node
    most_central_v_node = argmin(closeness_centrality(v_network))

    centrality_nodes = closeness_centrality(s_network)

    # ------ things for the PSO algorithm
    position = []
    velocity = []

    personal_best = []
    personal_best_cost = []

    global_best = nothing
    global_best_cost = Inf

    # initialization
    print_things && print("initialization... ")

    
    nb_particle_init = 0

    for mapping in initial_particles
        nb_particle_init += 1
        cost_mapping_overall = get_cost_placement(mapping) + get_cost_routing(mapping)
        push!(position, mapping.node_placement)
        push!(personal_best, mapping.node_placement)
        push!(personal_best_cost, cost_mapping_overall)
        if cost_mapping_overall < global_best_cost
            global_best = position[nb_particle_init]
            global_best_cost = cost_mapping_overall
            print_things && println("We got a new best solution! value $cost_mapping_overall")
        end
        push!(velocity, ones(nv(v_network)))
    end

    while time()-time_beginning<time_max && nb_particle_init < nb_particle
        partial_placement = zeros(Int, nv(v_network))

        init_s_nodes_scores = [ node_capacities[s_node] >= node_demands[most_central_v_node] ? centrality_nodes[s_node] + 1. * rand() : 10e9
                                    for s_node in vertices(s_network) ]
        
        s_node_start = argmin(init_s_nodes_scores)
        if minimum(init_s_nodes_scores) > 10e8
            return (mappings=[],
                    mapping_cost=Inf
            )
        end
        partial_placement[most_central_v_node] = s_node_start

        placement, placement_cost = complete_partial_placement(partial_placement; nb_nodes_to_try = nv(s_network_dir))
        if placement_cost < 10e7
            routing, routing_cost = shortest_path_routing(placement)
        else
            routing = nothing
            routing_cost = Inf
        end

        overall_cost = placement_cost + routing_cost

        if overall_cost < 10e8 
            nb_particle_init +=1
            push!(position, placement)
            push!(personal_best, placement)
            push!(personal_best_cost, overall_cost)

            if overall_cost < global_best_cost
                global_best = position[nb_particle_init]
                global_best_cost = overall_cost
                print_things && println("We got a new best solution! value $overall_cost")
            end

            push!(velocity, ones(nv(v_network)))
        end
    end



    print_things && println("Initialization done, $nb_particle_init particles, best is $global_best_cost, Starting iterations...")

    # iterations
    iter = 1
    time_total = 0
    while iter < nb_iter && time_total < time_max && nb_particle_init == nb_particle
        
        for particle in 1:nb_particle

            velocity[particle] = plus( velocity[particle], 
                                        minus(personal_best[particle], position[particle]), 
                                        minus(global_best, position[particle]))
            placement, placement_cost = times(position[particle], velocity[particle])


            if placement_cost < 10e7
                routing, routing_cost = shortest_path_routing(placement)
            else
                routing = nothing
                routing_cost = Inf
            end
                
            overall_cost = placement_cost + routing_cost

            if overall_cost < 10e7
                position[particle] = placement
            end

            if overall_cost < personal_best_cost[particle]
                personal_best[particle] = position[particle]
                personal_best_cost[particle] = overall_cost
            end
            if overall_cost < global_best_cost
                    global_best_cost = overall_cost
                    global_best = position[particle]
                    print_things && println("We got a new best solution! value $global_best_cost")
            end


        end
        
        

        iter += 1
        time_total=time()-time_beginning

    end 

    #println("Final best solution: $global_best")
    print_things && println("PSO finished at iteration $iter, finished in $(time()-time_beginning)s, best solution: $global_best_cost")

    if global_best_cost > 10e8
        return (mapping=nothing, mapping_cost=Inf)
    end

    final_mapping = nothing

    return (mapping = final_mapping, mapping_cost = global_best_cost)  
end

