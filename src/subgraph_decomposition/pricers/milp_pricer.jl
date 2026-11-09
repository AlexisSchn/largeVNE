using JuMP, CPLEX


#############============= PRICERRRR

struct SubProblem
    model
    s_network
    v_network
    subgraph
end


function set_up_pricer(instance, subgraph)

    s_network = instance.s_network
    s_network_dir = instance.s_network_dir
    
    #### Model
    model = Model(CPLEX.Optimizer)

    ### Variables
    @variable(model, x[v_node in vertices(subgraph.graph), s_node in vertices(s_network)], binary=true);
    @variable(model, y[v_edge in edges(subgraph.graph), s_edge in edges(s_network_dir)], binary=true);


    ### Constraints

    ## Nodes

    # one substrate node per virtual node
    for v_node in vertices(subgraph.graph)
        @constraint(model, sum(x[v_node, s_node] for s_node in vertices(s_network)) == 1)
    end

    # if one to one : one virtual node per substrate node
    for s_node in vertices(s_network)
        @constraint(model, sum(x[v_node, s_node] for v_node in vertices(subgraph.graph)) <= 1)
    end



    # node capacity
    for s_node in vertices(s_network)
        @constraint(model, 
            sum( x[v_node, s_node] 
                for v_node in vertices(subgraph.graph) ) 
            <= 
            instance.s_network[s_node][:cap] )
    end


    ## Edges 
    
    # edge capacity (undirected version)
    for s_edge in edges(s_network)
        @constraint(model, 
        sum( (y[v_edge, s_edge] + y[v_edge, get_reverse_edge(s_network_dir, s_edge)]  )
                for v_edge in edges(subgraph.graph)) 
            <= 
            s_network[src(s_edge), dst(s_edge)][:cap] )
    end
    
    # Flow conservation
    for s_node in vertices(s_network)
        for v_edge in edges(subgraph.graph)
            @constraint(model, 
                x[src(v_edge), s_node] - x[dst(v_edge), s_node] 
                ==
                sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) - 
                    sum(y[v_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node))
            )
        end
    end


    ## Departure cst : Node + Edge
    for s_node in vertices(s_network)
        for v_node in vertices(subgraph.graph)
            for v_edge in get_out_edges(subgraph.graph, v_node)
                @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) >= x[v_node, s_node])
            end
        end
    end

    # Outgoing edges cap: think about the original virtual graph!
    for v_node in vertices(subgraph.graph)
        for s_node in vertices(s_network)
            necessary_bw = degree(instance.v_network, subgraph.nodes_of_main_graph[v_node])
            s_edges_incident = [get_edge(s_network, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
            available_bw = sum(s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident;init=0.0)
            if necessary_bw > available_bw
                @constraint(model, model[:x][v_node, s_node] == 0)
            end 
        end
    end

    return SubProblem(model, instance.s_network, instance.v_network, subgraph);
end




function update_solve_pricer(instance, vn_decompo, pricer, dual_costs; time_limit = 1000)

    model = pricer.model
    subgraph = pricer.subgraph
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    ### Objective
    placement_cost = @expression(model, 
        sum( ( s_network[s_node][:cost] - dual_costs.capacity_s_node[s_node] ) * model[:x][v_node, s_node] 
            for v_node in vertices(subgraph.graph) for s_node in vertices(s_network) ))

    routing_cost = @expression(model, sum( 
        ( s_network[src(s_edge), dst(s_edge)][:cost] - dual_costs.capacity_s_edge[s_edge] ) 
        *  (model[:y][v_edge, s_edge] + model[:y][v_edge, get_reverse_edge(s_network_dir, s_edge)])
                for v_edge in edges(subgraph.graph) for s_edge in edges(s_network) ))


            
    # flow conservation
    flow_conservation_cost = AffExpr(0.)

    for s_node in vertices(s_network)
        for connecting_edge in vn_decompo.v_edges_master
            if subgraph ∈ keys(vn_decompo.v_nodes_assignment[src(connecting_edge)])
                v_node_subgraph = vn_decompo.v_nodes_assignment[src(connecting_edge)][subgraph]
                add_to_expression!(
                    flow_conservation_cost, 
                    -dual_costs.flow_conservation[connecting_edge][s_node] , 
                    model[:x][v_node_subgraph, s_node])
            end
            if subgraph ∈ keys(vn_decompo.v_nodes_assignment[dst(connecting_edge)])
                v_node_subgraph = vn_decompo.v_nodes_assignment[dst(connecting_edge)][subgraph]
                add_to_expression!(
                    flow_conservation_cost, 
                    +dual_costs.flow_conservation[connecting_edge][s_node], 
                    model[:x][v_node_subgraph, s_node])
            end
        end
    end


    # departure !
    departure_costs = AffExpr(0.)
    for s_node in vertices(s_network)
        for connecting_edge in vn_decompo.v_edges_master
            if subgraph ∈ keys(vn_decompo.v_nodes_assignment[src(connecting_edge)])
                v_node_subgraph = vn_decompo.v_nodes_assignment[src(connecting_edge)][subgraph]
                add_to_expression!(
                    departure_costs, 
                    -dual_costs.departure[connecting_edge][s_node], 
                    model[:x][v_node_subgraph,s_node])
            end
        end
    end


    @objective(model, Min, 
            -dual_costs.convexity[subgraph]
            + placement_cost + routing_cost 
            + flow_conservation_cost 
            + departure_costs);

    set_silent(model)
    set_time_limit_sec(model, time_limit)
    optimize!(model)

    status = primal_status(model)
    if (status != MOI.FEASIBLE_POINT)
        println("Infeasible subproblem... $status")
        return nothing, 10e6, 10e6
    end

    # Get the solution
    x_values = value.(model[:x])
    y_values = value.(model[:y])
    cost_of_column = 0.

    node_placement = []
    for v_node in vertices(subgraph.graph)
        for s_node in vertices(s_network)
            if x_values[v_node, s_node] > 0.99
                append!(node_placement, s_node)
                cost_of_column += pricer.s_network[s_node][:cost]
            end
        end
    end


    edge_routing = Dict()
    for v_edge in edges(subgraph.graph)
        if node_placement[src(v_edge)] == node_placement[dst(v_edge)]
            edge_routing[v_edge] = Path(src(v_edge), dst(v_edge), [], 0)
        end
        used_edges = []
        for s_edge in edges(s_network_dir)
            if y_values[v_edge, s_edge] > 0.99
                push!(used_edges, s_edge)
                cost_of_column += s_network_dir[src(s_edge), dst(s_edge)][:cost]
            end
        end
        edge_routing[v_edge] = order_path(s_network_dir, used_edges, node_placement[src(v_edge)], node_placement[dst(v_edge)]) 
    end
    mapping = Mapping(subgraph.graph, s_network_dir, node_placement, edge_routing)
    #println(mapping)
    #println(cost_of_column)

    dual_value = objective_value(model)
    #println("The price of the column is : $(cost_of_column), and the obj is : $(dual_value)")

    return mapping, cost_of_column, dual_value


end
