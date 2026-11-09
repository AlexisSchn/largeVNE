using Graphs, MetaGraphsNext
using Statistics
using  NetworkLayout
using Random
using StatsBase

includet("../utils/import_utils.jl")





function graphical_pso(instance; nb_particle = 25, nb_iter = 50)




    function shortest_path_routing(v_node_placement, get_routing=false)

        capacities_edges_copy = copy(capacities_edges)
        s_network_copy_dir_copy_ofgraph = nothing
        is_still_original_s_network = true

        edge_routing = Dict{Edge, Path}()
        overall_routing_costs = 0

        for v_edge in edges(v_network)

            s_src = v_node_placement[src(v_edge)]
            s_dst = v_node_placement[dst(v_edge)]

            nodes_of_path = Int[]
            cost_of_routing_current_edge = 0

            # Get the shortest path
            if is_still_original_s_network
                nodes_of_path = [s_dst]
                v = s_dst
                while v != s_src
                    u = shortest_paths.parents[s_src, v]
                    push!(nodes_of_path, u)
                    v = shortest_paths.parents[s_src, v]
                end
                reverse!(nodes_of_path)
                cost_of_routing_current_edge = shortest_paths.dists[s_src, s_dst]
            else # A smarter thing to do here, would be to use the basic path, and if it's using a removed edge, then compute the astar path. But it does not cost so much time so.
                edges_of_path = a_star(s_network_copy_dir_copy_ofgraph, s_src, s_dst, distmx)
                if edges_of_path == [] # No paths found: the graph is full!
                    return Dict(), 10e9
                end
                push!(nodes_of_path, src(edges_of_path[1]))
                for edge in edges_of_path
                    push!(nodes_of_path, dst(edge))
                    cost_of_routing_current_edge += distmx[src(edge), dst(edge)]
                end
            end
            
            overall_routing_costs += cost_of_routing_current_edge

            # Removing the capacities. If no more capacities, removing the edge.
            for i_node in 1:length(nodes_of_path)-1
                current_src = nodes_of_path[i_node]
                current_dst = nodes_of_path[i_node+1]
                capacities_edges_copy[(current_src), (current_dst)] -= 1
                capacities_edges_copy[(current_dst), (current_src)] -= 1

                if capacities_edges_copy[(current_src), (current_dst)] == 0
                    if is_still_original_s_network
                        s_network_copy_dir_copy_ofgraph = deepcopy(instance.s_network_dir.graph)
                        is_still_original_s_network = false
                    end
                    rem_edge!(s_network_copy_dir_copy_ofgraph, current_src, current_dst)
                    rem_edge!(s_network_copy_dir_copy_ofgraph, current_dst, current_src)
                end
            end

            # Sometime it is necessary to get the actual routing, but not all the time. 
            # Since it takes some time to retrieve the edges, due to poor design choices, I only do it when necessary.
            if get_routing
                edges_of_path = Edge[]
                for i_node in 1:length(nodes_of_path)-1
                    current_src = nodes_of_path[i_node]
                    current_dst = nodes_of_path[i_node+1]
                    push!(edges_of_path, get_edge(s_network_dir, current_src, current_dst))
                end
                path = Path(s_src, s_dst, edges_of_path, cost_of_routing_current_edge)
                edge_routing[v_edge] = path
            end
        end

        return edge_routing, overall_routing_costs
    end



    

    function complete_partial_placement(partial_placement; nb_nodes_to_try=10)

        placement = copy(partial_placement)

        #println("Placement beginning : $partial_placement")
        already_placed_v_nodes = filter(v_node -> placement[v_node] != 0, vertices(v_network))

        # Extreme case where no nodes are kept: we start from a random node again
        if isempty(already_placed_v_nodes)
            placement[most_central_v_node] = rand(capacited_nodes)
            push!(already_placed_v_nodes, most_central_v_node)
        end

        next_v_nodes = setdiff(
            union([neighbors(v_network, v) for v in already_placed_v_nodes]...),
            already_placed_v_nodes,
        )

        possible_s_nodes = filter(s_node -> s_node ∉ placement, capacited_nodes)

        while !isempty(next_v_nodes)

            # Take a node of the list
            shuffle!(next_v_nodes)
            v_node = popfirst!(next_v_nodes)
            #println("Next node to look at: $v_node !")
            # Get neighbors already placed
            placement_neighbors = [placement[s_neigh] for s_neigh in filter(v_neighbor -> placement[v_neighbor] != 0, neighbors(v_network, v_node))]

            # Choose 10 random nodes (or less if the network is not that big)
            number_s_nodes = min(length(possible_s_nodes), nb_nodes_to_try)
            some_s_nodes = sample(possible_s_nodes, number_s_nodes; replace=false)


            # NODE RANKING AND SELECTION
            distances = [ sum(shortest_paths.dists[s_src, s_node] for s_src in placement_neighbors) for s_node in some_s_nodes]
            distances_norm = (distances .- minimum(distances)) ./ (maximum(distances) - minimum(distances) + 1e-9)

            capacities = [ s_node_scores[s_node] for s_node in some_s_nodes]
            capacities_norm = (capacities .- minimum(capacities)) ./ (maximum(capacities) - minimum(capacities) + 1e-9)

            costs = [ node_costs[s_node] for s_node in some_s_nodes]
            costs_norm = (costs .- minimum(costs)) ./ (maximum(costs) - minimum(costs) + 1e-9)

            final_scores = 5. .* distances_norm .+ 0. .* costs_norm .+ 0. .* ( 1 .-capacities_norm) .+ 0. .* rand(length(some_s_nodes))
            selected_idx = argmin(final_scores)
            s_node_selected = some_s_nodes[selected_idx]

            #println("Placed on: $s_node_selected, that has score: $(distances[selected_idx]) !")
            # Finish the work
            placement[v_node] = s_node_selected
            push!(already_placed_v_nodes, v_node)
            next_v_nodes = next_v_nodes ∪ filter(v_neighbor->v_neighbor ∉ already_placed_v_nodes, neighbors(v_network, v_node) )
            possible_s_nodes = filter(!=(s_node_selected), possible_s_nodes)
        end

        #println("Placement end : $placement")


        placement_cost = 0
        for v_node in vertices(v_network)
            placement_cost += node_costs[placement[v_node]] 
        end

        return placement, placement_cost

    end




    println("FINEEEEE, let's do the graphical thing.")
    time_beginning = time()

    
    time_beginning = time()

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir



        
    #---- Usefull things
    node_capacities = [get_attribute_node(s_network, s_node, :cap) for s_node in vertices(s_network)]
    node_costs = [get_attribute_node(s_network, s_node, :cost) for s_node in vertices(s_network)]
    capacities_edges = Dict{Tuple{Int, Int}, Int}()
    for s_edge in edges(s_network)
        cap = s_network[src(s_edge), dst(s_edge)][:cap]
        capacities_edges[(src(s_edge),dst(s_edge))] = cap
        capacities_edges[(dst(s_edge),src(s_edge))] = cap
    end




    #---- Make sure there are enough capacited nodes
    capacited_nodes = [s_node for s_node in vertices(s_network) if node_capacities[s_node] ≥ 1]

    if length(capacited_nodes) < nv(v_network)
        println("What the hell? Not enough capacited nodes...")
        return  (mapping=nothing, 
                mapping_cost=10e9
        )
    end
    
    # Greedy score based on capacities for nodes
    s_node_scores = [ node_capacities[s_node]  * 
                        sum(capacities_edges[(s_node,s_neighbor)] for s_neighbor in neighbors(s_network, s_node) ) 
                    for s_node in vertices(s_network)]


    # shortest path of the substrate network
    distmx = zeros(Int, nv(s_network), nv(s_network))
    for s_edge in edges(s_network_dir)
        distmx[src(s_edge), dst(s_edge)] = get_attribute_edge(s_network_dir, s_edge, :cost)
    end
    shortest_paths = floyd_warshall_shortest_paths(s_network_dir, distmx)



    neighbors_sn_memory = [neighbors(s_network, s_node) for s_node in vertices(s_network)]
    neighbors_vn_memory = [neighbors(v_network, v_node) for v_node in vertices(v_network)]

    # Greedy score based on capacities for nodes: not for now.
    
    s_node_scores = [ node_capacities[s_node]  * 
                        sum(capacities_edges[(s_node,s_neighbor)] for s_neighbor in neighbors_sn_memory[s_node]) 
                    for s_node in vertices(s_network)]
    

    # shortest path of the substrate network
    distmx = zeros(Int, nv(s_network), nv(s_network))
    for s_edge in edges(s_network_dir)
        distmx[src(s_edge), dst(s_edge)] = get_attribute_edge(s_network_dir, s_edge, :cost)
    end
    shortest_paths = floyd_warshall_shortest_paths(s_network_dir, distmx)


    # Getting the most central virtual node
    most_central_v_node = argmin(closeness_centrality(v_network))


    # Algos for layout: stress, spring, 
    coords_s_nodes = stress(instance.s_network; iterations=500)
    coords_s_nodes = normalize_coords(coords_s_nodes)

    pos_particles = []
    vel_particles = []
    best_pos_particles = []
    best_cost_particles = []

    history_pos = []

    best_pos_overall = (0., 0.)
    best_cost_overall = 10e9

    pos_originales = init_particles_grid(nb_particle)

    # initialization
    for particle in 1:nb_particle

        # put a random position
        
        #central_s_node = rand(capacited_nodes)
        #pos = coords_s_nodes[central_s_node]
        
        starting_pos = pos_originales[particle]
        distances = [sqrt((starting_pos[1]-x)^2 + (starting_pos[2]-y)^2) for (x,y) in coords_s_nodes] # NEED TO CHANGE THIS!!!!! TO CAPACITED NODE ONLY
        central_s_node = argmin(distances)
        pos = coords_s_nodes[central_s_node]

        vel = ((rand()*2-1)*0.15, (rand()*2-1)*0.15)
        
        partial_placement = zeros(Int, nv(v_network))
        partial_placement[most_central_v_node] = central_s_node
                
        node_placement, placement_cost = complete_partial_placement(partial_placement, nb_nodes_to_try=length(capacited_nodes))
        edge_routing, routing_cost = shortest_path_routing(node_placement)
        total_cost = placement_cost + routing_cost


        if total_cost < best_cost_overall
            best_pos_overall = pos
            best_cost_overall = total_cost
        end

        push!(pos_particles, pos)
        push!(vel_particles, vel)    
        
        push!(best_pos_particles, pos)
        push!(best_cost_particles, total_cost)

        push!(history_pos, [pos])
    end

    println("Best after init: $best_cost_overall")

    for iter in 1:nb_iter

        for particle in 1:nb_particle

            current_vel = vel_particles[particle]
            current_pos = pos_particles[particle]
            
            new_vel =   0.7 .* current_vel .+ 
                        1.5 * rand() .* (best_pos_particles[particle] .- current_pos) .+ 
                        1.5 * rand() .* (best_pos_overall .- current_pos);
            
            new_vel = clamp.(new_vel, -0.05, 0.05)

            new_pos = new_vel .+ current_pos
            new_pos = clamp.(new_pos, -1, 1.)

            # get closest substrate node to the pos:
            distances = [sqrt((new_pos[1]-x)^2 + (new_pos[2]-y)^2) for (x,y) in coords_s_nodes] # NEED TO CHANGE THIS!!!!! TO CAPACITED NODE ONLY
            central_s_node = argmin(distances)
            partial_placement = zeros(Int, nv(v_network))
            partial_placement[most_central_v_node] = central_s_node
                    
            node_placement, placement_cost = complete_partial_placement(partial_placement, nb_nodes_to_try=ceil(Int,length(capacited_nodes)))

            edge_routing, routing_cost = shortest_path_routing(node_placement)
    
            total_cost = placement_cost + routing_cost


            if total_cost < best_cost_particles[particle]
                best_pos_particles[particle] = coords_s_nodes[central_s_node]
                best_cost_particles[particle] = total_cost
            end
            if total_cost < best_cost_overall
                best_pos_overall = coords_s_nodes[central_s_node]
                best_cost_overall = total_cost
                println("New best found! $best_cost_overall, at iter $iter")
            end
    
            pos_particles[particle] = new_pos
            vel_particles[particle] = new_vel
            
            push!(history_pos[particle], new_pos)

        end

        #visu_particle(s_network.graph, coords_s_nodes, pos_particles, iter)

    end


    #print_history_pos(history_pos, nb_particle, nb_iter)

    println("Final best cost: $best_cost_overall, took me $(time() - time_beginning)")

    return (mapping_cost = best_cost_overall)
end



function normalize_coords(coords)
    xs = [c[1] for c in coords]
    ys = [c[2] for c in coords]
    xmin, xmax = minimum(xs), maximum(xs)
    ymin, ymax = minimum(ys), maximum(ys)

    return [(2*(x - xmin)/(xmax - xmin) - 1,
             2*(y - ymin)/(ymax - ymin) - 1) for (x,y) in coords]
end


function print_history_pos(history_pos, nb_particle, nb_iter)

    #=
    for particle in 1:nb_particle

        println("For particle $particle:")
        for pos in history_pos[particle]
            println("( $(round(pos[1], digits=3)), $(round(pos[2], digits=3)))")
        end

    end
    =#
    println("End position:")
    for particle in 1:nb_particle
        pos = history_pos[particle][nb_iter]
        println("( $(round(pos[1], digits=3)), $(round(pos[2], digits=3)))")
    end

end



function init_particles_grid(num_particles::Int)
    rows = floor(Int, sqrt(num_particles))
    cols = ceil(Int, num_particles / rows)

    println("Number of rows: $rows, number of cols: $cols")
    step_x =2. / cols
    step_y = 2. / rows

    positions = []
    count = 0
    for i in 0:(cols-1)
        for j in 0:(rows-1)
            push!(positions, ( -1 + (i + 0.5) * step_x, -1 +  (j + 0.5) * step_y))
            count += 1
            if count == num_particles
                return positions
            end
        end
    end
    return positions
end

