
using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX




function local_search(instance, vn_decompo, solution)

    s_network = instance.s_network 
    s_network_dir = instance.s_network_dir
    v_network = instance.v_network

    time0 = time()
    iter = 1
    iter_max = 100
    i_subgraph = 1
    stop_maybe = 0
    keep_on = true
    nb_subgraph = length(vn_decompo.subgraphs)

    best_sol_val = get_cost_placement(solution) + get_cost_routing(solution)
    best_sol_placement = solution.node_placement

    println("Starting local search...")

    while keep_on
        println("New iteration...")
        subgraph = vn_decompo.subgraphs[i_subgraph]

        # get node_placement, not taking into account the one subgraph
        placement_restriction = Dict()
        for v_node in vertices(v_network)
            if !(v_node in subgraph.nodes_of_main_graph)
                placement_restriction[v_node] = best_sol_placement[v_node]
            end
        end

        model = Model(CPLEX.Optimizer)
        set_up_problem_restricted(instance, model, placement_restriction)
        set_silent(model)
        set_time_limit_sec(model, 50) # 10 s should be enough? IDK with all the routing
        optimize!(model)

        sol_val = objective_value(model)

        if sol_val < best_sol_val
            println("New best sol found! $sol_val")
            best_sol_val = sol_val
            stop_maybe = 0

            # Get the solution
            best_sol_placement = []
            x_values = value.(model[:x])        
            for v_node in vertices(v_network)
                for s_node in vertices(s_network)
                    if x_values[v_node, s_node] > 0.01
                        push!(best_sol_placement, s_node)
                    end
                end
            end
    
        else
            stop_maybe += 1
        end


        i_subgraph += 1
        if i_subgraph > nb_subgraph
            i_subgraph = 1
        end

        if (iter > iter_max) || (stop_maybe >= nb_subgraph)
            keep_on = false
        end
    end


    println("At the end, the best solution found is $(best_sol_val)")
    time_local_search = time() - time0
    println("Took $time_local_search in local search")
end





function local_search_changin(instance, solution, time_max)

    println("Local search, on a-changin subgraphs!")

    best_sol_val = get_cost_placement(solution) + get_cost_routing(solution)
    best_sol_placement = solution.node_placement

    v_network = instance.v_network
    s_network = instance.s_network

    
    time_beginning = time()
    while (time() - time_beginning) < time_max

        # Get the subgraph we will do
        central_v_node = rand(1:nv(v_network))  # use rand instead of random
        v_nodes_to_change = Set([central_v_node])

        # First neighbors
        for neighbor in neighbors(v_network, central_v_node)
            push!(v_nodes_to_change, neighbor)
        end

        # Add second-level neighbors if still under 10 nodes
        for node in copy(v_nodes_to_change)
            for neighbor in neighbors(v_network, node)
                push!(v_nodes_to_change, neighbor)
                if length(v_nodes_to_change) >= 10
                    break
                end
            end
        end

        # Convert to Vector if needed
        v_nodes_to_change = collect(v_nodes_to_change)

        #println("Selected nodes: $v_nodes_to_change")




        # ===== get node_placement, not taking into account the selected nodes
        placement_restriction = Dict()
        for v_node in vertices(v_network)
            if !(v_node in v_nodes_to_change)
                placement_restriction[v_node] = best_sol_placement[v_node]
            end
        end

        model = Model(CPLEX.Optimizer)
        set_up_problem_restricted(instance, model, placement_restriction)
        set_silent(model)
        set_time_limit_sec(model, 50) # 10 s should be enough? IDK with all the routing
        optimize!(model)

        sol_val = objective_value(model)

        if sol_val < best_sol_val
            println("New best sol found! $sol_val")
            best_sol_val = sol_val

            # Get the solution
            best_sol_placement = []
            x_values = value.(model[:x])        
            for v_node in vertices(v_network)
                for s_node in vertices(s_network)
                    if x_values[v_node, s_node] > 0.01
                        push!(best_sol_placement, s_node)
                    end
                end
            end
        end
        


    end



    return best_sol_val, best_sol_placement

end



# ========== CLASSICAL STUFF
function set_up_problem_restricted(instance, model, placement_restriction)

    v_network = instance.v_network
    s_network_dir = instance.s_network_dir
    s_network = instance.s_network

    ### Variables
    @variable(model, x[vertices(v_network), vertices(instance.s_network)], binary=true);
    @variable(model, y[edges(v_network), edges(s_network_dir)], binary=true);

    

    ### Objective
    placement_cost = @expression(model, sum( instance.s_network[s_node][:cost] * x[v_node, s_node] 
        for v_node in vertices(v_network) for s_node in vertices(instance.s_network) ))
    routing_cost = @expression(model, sum( s_network_dir[src(s_edge), dst(s_edge)][:cost] * y[v_edge, s_edge]
        for v_edge in edges(v_network) for s_edge in edges(s_network_dir) ))
    @objective(model, Min, placement_cost + routing_cost);




    ###=========== Constraints

    ##---- Nodes

    # one substrate node per virtual node
    for v_node in vertices(v_network)
        @constraint(model, sum(x[v_node, s_node] for s_node in vertices(instance.s_network)) == 1)
    end

    # one to one : one virtual node per substrate node
    for s_node in vertices(instance.s_network)
        @constraint(model, sum(x[v_node, s_node] for v_node in vertices(v_network)) <= 1)
    end

    # node capacity : NOT USELESS AHHHHHHHHh
    for s_node in vertices(instance.s_network)
        @constraint(model, sum( x[v_node, s_node] for v_node in vertices(v_network)) <= sum(s_network[s_node][:cap]))
    end

    # respect the placement_restriction pleaseee
    for (v_node, s_node) in placement_restriction
        @constraint(model, x[v_node, s_node] == 1)
    end

    ##--- Edges 
    
    # edge capacity (undirected version !)
    for s_edge in edges(instance.s_network)
        @constraint(model, 
            sum( (y[v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] + y[v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]  )
                for v_edge in edges(v_network)) 
            <= 
            instance.s_network[src(s_edge), dst(s_edge)][:cap] )
    end
    
    # Flow conservation
    for s_node in vertices(instance.s_network)
        for v_edge in edges(v_network)
            @constraint(model, 
                x[src(v_edge), s_node] - x[dst(v_edge), s_node] 
                ==
                sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) - 
                    sum(y[v_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node))
            )
        end
    end

    
    ## Departure constraints    
    for s_node in vertices(instance.s_network)
        for v_edge in edges(v_network)
            @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) 
                >= x[src(v_edge), s_node])
        end
    end
    
    
    
    # Star capacity constraint
    for v_node in vertices(v_network)
        necessary_bw = degree(v_network, v_node) 
        
        for s_node in vertices(s_network)
            s_edges_incident = [get_edge(s_network, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
            available_bw = sum(s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident;init=0.)
            if necessary_bw > available_bw
                @constraint(model, model[:x][v_node, s_node] == 0)
            end 
        end
    end
        
end

