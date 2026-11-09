using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX




function diving_heuristic(instance, vn_decompo, master_problem, time_max)
    
    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    v_subgraphs = vn_decompo.subgraphs

    

    columns_selected = []

    for i_subgraph in 1:(length(v_subgraphs))

        println("\nTime to do it for subgraph $i_subgraph !\n\n")

        model = Model(CPLEX.Optimizer)
        set_up_problem_linear(instance, vn_decompo, model)
        set_time_limit_sec(model, time_max)    


        current_subgraph = v_subgraphs[i_subgraph]

        for i_previous_subgraph in 1:(i_subgraph-1)
            add_column_ip(model, vn_decompo, instance, v_subgraphs[i_previous_subgraph], columns_selected[i_previous_subgraph])
        end

        lambdas = Dict()
        for column in master_problem.columns[current_subgraph]
            lambda = add_column_ip(model, vn_decompo, instance, current_subgraph, column)
            lambdas[column] = lambda
        end

        for i_later_subgraph in (i_subgraph+1):length(v_subgraphs)
            for column in master_problem.columns[v_subgraphs[i_later_subgraph]]
                add_column_linear(model, vn_decompo, instance, v_subgraphs[i_later_subgraph], column)
            end
        end

        optimize!(model)

        status = primal_status(model)
        if status != MOI.FEASIBLE_POINT
            println("Infeasible or unfinished: $status")
            return -999., nothing
        end
        println("Optimal solution : $(objective_value(model))")

        

        for column in keys(lambdas)
            if value(lambdas[column]) > 0.9 # This column is selected
                push!(columns_selected, column)
            end
        end

    end
    

    # Let's do the final mapping. Maybe let re-routing happen here? It does not cost much.

    model = Model(CPLEX.Optimizer)
    set_up_problem_integer(instance, vn_decompo, model)
    set_time_limit_sec(model, time_max)    

    for (i_subgraph, column) in enumerate(columns_selected)
        add_column_ip(model, vn_decompo, instance, v_subgraphs[i_subgraph], column)
    end

    optimize!(model)
    status = primal_status(model)
    if status != MOI.FEASIBLE_POINT
        println("Infeasible or unfinished: $status")
        return -999., nothing
    end
    println("Optimal solution : $(objective_value(model))")

    # ---- GET SOLUTION
    node_placement = zeros(Integer, nv(instance.v_network))
    edge_routing = Dict()

    # Column placement and routing
    for (i_subgraph, subgraph) in enumerate(v_subgraphs)
    
        column = columns_selected[i_subgraph]

        for v_node in vertices(subgraph.graph)
            original_node = subgraph.nodes_of_main_graph[v_node]
            node_placement[original_node] = column.mapping.node_placement[v_node]
        end
        for v_edge in edges(subgraph.graph)
            original_src = subgraph.nodes_of_main_graph[src(v_edge)]
            original_dst = subgraph.nodes_of_main_graph[dst(v_edge)]
            original_edge = get_edge(v_network, original_src, original_dst)
            edge_routing[original_edge] = column.mapping.edge_routing[v_edge]
        end

    end

    # Master problem routing
    y_values = value.(model[:y])
    for v_edge in vn_decompo.v_edges_master
        used_edges = []
        for s_edge in edges(s_network_dir)
            if y_values[v_edge, s_edge] > 0.9
                push!(used_edges, s_edge)
            end
        end
        edge_routing[v_edge] = order_path(s_network_dir, used_edges, node_placement[src(v_edge)], node_placement[dst(v_edge)]) 
    end


    mapping = Mapping(v_network, s_network, node_placement, edge_routing)

    return objective_value(model), mapping
    
end





function set_up_problem_linear(instance, vn_decompo, model)

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    #set_silent(model)

    
    ### Variables    
    @variable(model,  0 <= y[v_edge in vn_decompo.v_edges_master, s_edge in edges(s_network_dir)] <= 1 );
    

    columns = Dict()
    for subgraph in vn_decompo.subgraphs
        columns[subgraph] = []
    end
    
    

    ### Objective
    master_routing_costs = @expression(model, sum( s_network_dir[src(s_edge), dst(s_edge)][:cost] * y[v_edge, s_edge]
        for v_edge in vn_decompo.v_edges_master for s_edge in edges(s_network_dir) ))
    
    @objective(model, Min, master_routing_costs);

    ### Constraints

    # convexity constraints
    @constraint(
        model, 
        mapping_selec[subgraph in vn_decompo.subgraphs],
        0 >= 1
    );


    # capacity of substrate nodes
    @constraint(
        model,
        capacity_s_node[s_node in vertices(s_network)],
        0 <= s_network[s_node][:cap]
    );

    

    # capacity of substrate edges
    # undirected, so both ways !
    @constraint(
        model,
        capacity_s_edge[s_edge in edges(s_network)],
        sum( (y[v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] +  y[v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))] )
            for v_edge in vn_decompo.v_edges_master)
        + 0
        <= s_network[src(s_edge), dst(s_edge)][:cap]
    );


    # flow conservation constraints
    @constraint(
        model,
        flow_conservation[v_edge in vn_decompo.v_edges_master, s_node in vertices(s_network)],
        0 == 
        sum( y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node))
        - sum( y[v_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node))
    );


    # Departure constraints
    @constraint(
        model, 
        departure[v_edge in vn_decompo.v_edges_master, s_node in vertices(s_network)],
        0 
        <=
        sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node))
    )




end



function set_up_problem_integer(instance, vn_decompo, model)

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    #set_silent(model)

    
    ### Variables    
    @variable(model, y[v_edge in vn_decompo.v_edges_master, s_edge in edges(s_network_dir)], binary=true);
    

    columns = Dict()
    for subgraph in vn_decompo.subgraphs
        columns[subgraph] = []
    end
    
    

    ### Objective
    master_routing_costs = @expression(model, sum( s_network_dir[src(s_edge), dst(s_edge)][:cost] * y[v_edge, s_edge]
        for v_edge in vn_decompo.v_edges_master for s_edge in edges(s_network_dir) ))
    
    @objective(model, Min, master_routing_costs);

    ### Constraints

    # convexity constraints
    @constraint(
        model, 
        mapping_selec[subgraph in vn_decompo.subgraphs],
        0 >= 1
    );


    # capacity of substrate nodes
    @constraint(
        model,
        capacity_s_node[s_node in vertices(s_network)],
        0 <= s_network[s_node][:cap]
    );

    

    # capacity of substrate edges
    # undirected, so both ways !
    @constraint(
        model,
        capacity_s_edge[s_edge in edges(s_network)],
        sum( (y[v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] +  y[v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))] )
            for v_edge in vn_decompo.v_edges_master)
        + 0
        <= s_network[src(s_edge), dst(s_edge)][:cap]
    );


    # flow conservation constraints
    @constraint(
        model,
        flow_conservation[v_edge in vn_decompo.v_edges_master, s_node in vertices(s_network)],
        0 == 
        sum( y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node))
        - sum( y[v_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node))
    );


    # Departure constraints
    @constraint(
        model, 
        departure[v_edge in vn_decompo.v_edges_master, s_node in vertices(s_network)],
        0 
        <=
        sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node))
    )




end



function add_column_ip(model, vn_decompo, instance, subgraph, column)

    s_network = instance.s_network
    s_network_dir = instance.s_network_dir


    lambda = @variable(model, binary=true);
    set_objective_coefficient(model, lambda, column.cost)

    # convexity
    set_normalized_coefficient(model[:mapping_selec][subgraph], lambda, 1)

    # capacities
    for s_node in vertices(s_network)

        usage = 0
        for v_node in vertices(subgraph.graph)
            if column.mapping.node_placement[v_node] == s_node
                usage += 1
            end
        end

        set_normalized_coefficient(model[:capacity_s_node][s_node], lambda, usage)
    end

    # Undirected, so it is a bit tricky... (transforming undirected mapping to make it better ?)
    for s_edge in edges(s_network)
        usage = 0
        for v_edge in edges(subgraph.graph)
            s_edge_one = get_edge(s_network_dir, src(s_edge), dst(s_edge))
            if s_edge_one in column.mapping.edge_routing[v_edge].edges
                usage += 1
            end
            s_edge_two = get_edge(s_network_dir, dst(s_edge), src(s_edge))
            if s_edge_two in column.mapping.edge_routing[v_edge].edges
                usage += 1
            end
        end
        set_normalized_coefficient(model[:capacity_s_edge][s_edge], lambda, usage)
    end

    
    # flow conservation 
    for v_edge in vn_decompo.v_edges_master

        if subgraph in keys(vn_decompo.v_nodes_assignment[src(v_edge)])

            v_node_in_subgraph = vn_decompo.v_nodes_assignment[src(v_edge)][subgraph]

            set_normalized_coefficient(
                model[:flow_conservation][v_edge, column.mapping.node_placement[v_node_in_subgraph]], 
                lambda, 
                1)

        end

        if subgraph in keys(vn_decompo.v_nodes_assignment[dst(v_edge)])

            v_node_in_subgraph = vn_decompo.v_nodes_assignment[dst(v_edge)][subgraph]

            set_normalized_coefficient(
                model[:flow_conservation][v_edge, column.mapping.node_placement[v_node_in_subgraph]], 
                lambda, 
                -1 )

        end


    end
    
    # departure
    for v_edge in vn_decompo.v_edges_master

        if subgraph in keys(vn_decompo.v_nodes_assignment[src(v_edge)])

            v_node_in_subgraph = vn_decompo.v_nodes_assignment[src(v_edge)][subgraph]

            set_normalized_coefficient(
                model[:departure][v_edge, column.mapping.node_placement[v_node_in_subgraph]], 
                lambda, 
                1 )

        end
    end    
    
    return lambda


end



function add_column_linear(model, vn_decompo, instance, subgraph, column)

    s_network = instance.s_network
    s_network_dir = instance.s_network_dir


    lambda = @variable(model, lower_bound=0., upper_bound=1.);
    set_objective_coefficient(model, lambda, column.cost)

    # convexity
    set_normalized_coefficient(model[:mapping_selec][subgraph], lambda, 1)

    # capacities
    for s_node in vertices(s_network)

        usage = 0
        for v_node in vertices(subgraph.graph)
            if column.mapping.node_placement[v_node] == s_node
                usage += 1
            end
        end

        set_normalized_coefficient(model[:capacity_s_node][s_node], lambda, usage)
    end

    # Undirected, so it is a bit tricky... (transforming undirected mapping to make it better ?)
    for s_edge in edges(s_network)
        usage = 0
        for v_edge in edges(subgraph.graph)
            s_edge_one = get_edge(s_network_dir, src(s_edge), dst(s_edge))
            if s_edge_one in column.mapping.edge_routing[v_edge].edges
                usage += 1
            end
            s_edge_two = get_edge(s_network_dir, dst(s_edge), src(s_edge))
            if s_edge_two in column.mapping.edge_routing[v_edge].edges
                usage += 1
            end
        end
        set_normalized_coefficient(model[:capacity_s_edge][s_edge], lambda, usage)
    end

    
    # flow conservation 
    for v_edge in vn_decompo.v_edges_master

        if subgraph in keys(vn_decompo.v_nodes_assignment[src(v_edge)])

            v_node_in_subgraph = vn_decompo.v_nodes_assignment[src(v_edge)][subgraph]

            set_normalized_coefficient(
                model[:flow_conservation][v_edge, column.mapping.node_placement[v_node_in_subgraph]], 
                lambda, 
                1)

        end

        if subgraph in keys(vn_decompo.v_nodes_assignment[dst(v_edge)])

            v_node_in_subgraph = vn_decompo.v_nodes_assignment[dst(v_edge)][subgraph]

            set_normalized_coefficient(
                model[:flow_conservation][v_edge, column.mapping.node_placement[v_node_in_subgraph]], 
                lambda, 
                -1 )

        end


    end
    
    # departure
    for v_edge in vn_decompo.v_edges_master

        if subgraph in keys(vn_decompo.v_nodes_assignment[src(v_edge)])

            v_node_in_subgraph = vn_decompo.v_nodes_assignment[src(v_edge)][subgraph]

            set_normalized_coefficient(
                model[:departure][v_edge, column.mapping.node_placement[v_node_in_subgraph]], 
                lambda, 
                1 )

        end
    end    
    
    return lambda


end



