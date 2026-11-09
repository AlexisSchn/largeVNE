
using CPLEX, JuMP



struct PricerSubSubstrate
    vn_subgraph
    sn_subgraph
    sub_instance
    original_instance
    model
end

function Base.show(io::IO, pricer::PricerSubSubstrate)
    println(io, "Pricer for subvn $(pricer.vn_subgraph.graph[][:name]) on subsn $(pricer.sub_instance.s_network[][:name])")
end


function solve_pricer_milp_routing(v_subgraph, s_subgraph, original_instance, vn_decompo, additional_costs_routing, dual_costs; time_solver = 60)


    sub_instance = Instance(v_subgraph.graph, s_subgraph.graph)
    model = Model(CPLEX.Optimizer)
    set_up_model_pricer_subsn(model, sub_instance, original_instance, v_subgraph, s_subgraph)
    pricer =  PricerSubSubstrate(v_subgraph, s_subgraph, sub_instance, original_instance, model)

    update_pricer_sn_decompo_penalty(vn_decompo, pricer, dual_costs, additional_costs_routing)
    result = solve_pricers_sn_decompo_penalty(pricer, additional_costs_routing; time_limit = time_solver)
    
    return result
end



#=
function set_up_pricers_sn_partitionning(instance, vn_subgraphs, sn_subgraphs)

    pricers = Dict()
    for v_subgraph in vn_subgraphs
        pricers_v_subgraph = Dict()
        for s_subgraph in sn_subgraphs
            sub_instance = Instance(v_subgraph.graph, s_subgraph.graph)
            model = Model(CPLEX.Optimizer)
            set_up_model_pricer_subsn(model, sub_instance, instance, v_subgraph, s_subgraph)
            pricer =  PricerSubSubstrate(v_subgraph, s_subgraph, sub_instance, instance, model)
            pricers_v_subgraph[s_subgraph] = pricer
        end
        pricers[v_subgraph] = pricers_v_subgraph
    end

    return pricers
end
=#




function set_up_model_pricer_subsn(model, sub_instance, original_instance, vn_subgraph, sn_subgraph)


    s_network = sub_instance.s_network
    s_network_dir = sub_instance.s_network_dir
    v_network = sub_instance.v_network


    original_v_network = original_instance.v_network
    original_s_network = original_instance.s_network


    set_silent(model)
    ### Variables
    @variable(model, x[vertices(v_network), vertices(s_network)], binary=true);
    @variable(model, y[edges(v_network), edges(s_network_dir)], binary=true);
 

    ### Objective
    placement_cost = @expression(model, sum( s_network[s_node][:cost] * x[v_node, s_node] 
        for v_node in vertices(v_network) for s_node in vertices(s_network) ))
    routing_cost = @expression(model, sum( s_network_dir[src(s_edge), dst(s_edge)][:cost] * y[v_edge, s_edge]
        for v_edge in edges(v_network) for s_edge in edges(s_network_dir) ))
    @objective(model, Min, placement_cost + routing_cost);

    ### Constraints

    ## Nodes

    # one substrate node per virtual node
    for v_node in vertices(v_network)
        @constraint(model, sum(x[v_node, s_node] for s_node in vertices(s_network)) == 1)
    end

    # if one to one : one virtual node per substrate node
    for s_node in vertices(s_network)
        @constraint(model, sum(x[v_node, s_node] for v_node in vertices(v_network)) <= 1)
    end



    # node capacity
    for s_node in vertices(s_network)
        @constraint(model, 
            sum( x[v_node, s_node] 
                for v_node in vertices(v_network) ) 
            <= 
            s_network[s_node][:cap] )
    end


    ## Edges 
    
    # edge capacity (undirected version)
    for s_edge in edges(s_network)
        @constraint(model, 
        sum( (y[v_edge, s_edge] + y[v_edge, get_reverse_edge(s_network_dir, s_edge)]  )
                for v_edge in edges(vn_subgraph.graph)) 
            <= 
            s_network[src(s_edge), dst(s_edge)][:cap] )
    end
    
    # Flow conservation
    for s_node in vertices(s_network)
        for v_edge in edges(v_network)
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
        for v_node in vertices(v_network)
            for v_edge in get_out_edges(v_network, v_node)
                @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) >= x[v_node, s_node])
            end
        end
    end
    
    # Removing some x, due to sub_vn and sub_sn connectivity & capacities. THIS IS LOCAL TO THE PRICER!
    for v_node in vertices(v_network)
        necessary_bw = degree(v_network, v_node)
        for s_node in vertices(s_network)
            s_edges_incident = [get_edge(s_network, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
            available_bw = sum(s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident; init=0.0)
            if necessary_bw > available_bw
                @constraint(model, model[:x][v_node, s_node] == 0)
            end 
        end
    end

    # Removing some x, due to the overall v_network and s_network
    for v_node in vertices(v_network)
        necessary_bw = degree(original_v_network, vn_subgraph.nodes_of_main_graph[v_node])
        for s_node in vertices(s_network)
            s_edges_incident = [get_edge(original_s_network, sn_subgraph.nodes_of_main_graph[s_node], neighbor) for neighbor in neighbors(original_s_network, sn_subgraph.nodes_of_main_graph[s_node])]
            available_bw = sum(original_s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident; init=0.0)
            if necessary_bw > available_bw
                @constraint(model, model[:x][v_node, s_node] == 0)
            end 
        end
    end
    

end







function update_pricer_sn_decompo_penalty(vn_decompo, pricer, dual_costs, additional_costs)

    model = pricer.model

    sn_subgraph = pricer.sn_subgraph
    sub_s_network_dir = pricer.sub_instance.s_network_dir
    original_s_network = pricer.original_instance.s_network
    sub_s_network = pricer.sub_instance.s_network
    v_subgraph = pricer.vn_subgraph

    ### Objective
    placement_cost = @expression(model, 
        sum( ( sub_s_network_dir[s_node][:cost] - dual_costs.capacity_s_node[sn_subgraph.nodes_of_main_graph[s_node]] + additional_costs[v_node][s_node])  * model[:x][v_node, s_node] 
            for v_node in vertices(v_subgraph.graph) for s_node in vertices(sub_s_network_dir) ))

    routing_cost = @expression(model, sum( 
        ( sub_s_network[src(s_edge), dst(s_edge)][:cost] - dual_costs.capacity_s_edge[get_edge(original_s_network, sn_subgraph.nodes_of_main_graph[src(s_edge)], sn_subgraph.nodes_of_main_graph[dst(s_edge)])] ) 
        * v_subgraph.graph[src(v_edge), dst(v_edge)][:dem] * (model[:y][v_edge, get_edge(sub_s_network_dir, src(s_edge), dst(s_edge))] + model[:y][v_edge, get_edge(sub_s_network_dir, dst(s_edge), src(s_edge))])
                for v_edge in edges(v_subgraph.graph) for s_edge in edges(sub_s_network) ))


            
    # flow conservation
    flow_conservation_cost = AffExpr(0.)

    for s_node in vertices(sub_s_network)
        original_node = sn_subgraph.nodes_of_main_graph[s_node]
        for connecting_edge in vn_decompo.v_edges_master
            if v_subgraph ∈ keys(vn_decompo.v_nodes_assignment[src(connecting_edge)])
                v_node_subgraph = vn_decompo.v_nodes_assignment[src(connecting_edge)][v_subgraph]
                add_to_expression!(
                    flow_conservation_cost, 
                    -dual_costs.flow_conservation[connecting_edge][original_node] , 
                    model[:x][v_node_subgraph, s_node])
            end
            if v_subgraph ∈ keys(vn_decompo.v_nodes_assignment[dst(connecting_edge)])
                v_node_subgraph = vn_decompo.v_nodes_assignment[dst(connecting_edge)][v_subgraph]
                add_to_expression!(
                    flow_conservation_cost, 
                    +dual_costs.flow_conservation[connecting_edge][original_node], 
                    model[:x][v_node_subgraph, s_node])
            end
        end
    end


    # departure !
    departure_costs = AffExpr(0.)
    for s_node in vertices(sub_s_network)
        for connecting_edge in vn_decompo.v_edges_master
            if v_subgraph ∈ keys(vn_decompo.v_nodes_assignment[src(connecting_edge)])
                v_node_subgraph = vn_decompo.v_nodes_assignment[src(connecting_edge)][v_subgraph]
                add_to_expression!(
                    departure_costs, 
                    -dual_costs.departure[connecting_edge][sn_subgraph.nodes_of_main_graph[s_node]], 
                    model[:x][v_node_subgraph, s_node])
            end
        end
    end


    @objective(model, Min, 
            -dual_costs.convexity[v_subgraph]
            + placement_cost + routing_cost 
            + flow_conservation_cost 
            + departure_costs);


    return
end




function solve_pricers_sn_decompo_penalty(pricer, additional_costs; time_limit = 60)

    original_s_network_dir = pricer.original_instance.s_network_dir
    sn_subgraph = pricer.sn_subgraph
    vn_subgraph = pricer.vn_subgraph
    v_network = vn_subgraph.graph
    s_network_dir = pricer.sub_instance.s_network_dir
    model = pricer.model

    set_time_limit_sec(model, time_limit)
    optimize!(model)
    status = termination_status(model)

    if status == MOI.FEASIBLE_POINT || status == MOI.OPTIMAL
        # Get the solution
        x_values = value.(model[:x])
        y_values = value.(model[:y])
        penalty_routing = 0
        true_cost = 0.

        node_placement = []
        for v_node in vertices(v_network)
            for s_node in vertices(s_network_dir)
                if x_values[v_node, s_node] > 0.99
                    real_s_node = sn_subgraph.nodes_of_main_graph[s_node]
                    append!(node_placement, real_s_node)
                    true_cost += original_s_network_dir[real_s_node][:cost]
                    penalty_routing+= additional_costs[v_node][s_node]
                end
            end
        end


        edge_routing = Dict()
        for v_edge in edges(v_network)
            if node_placement[src(v_edge)] == node_placement[dst(v_edge)]
                edge_routing[v_edge] = Path(src(v_edge), dst(v_edge), [], 0)
            end
            used_edges = []
            for s_edge in edges(s_network_dir)
                if y_values[v_edge, s_edge] > 0.99
                    real_s_edge = get_edge(original_s_network_dir, sn_subgraph.nodes_of_main_graph[src(s_edge)], sn_subgraph.nodes_of_main_graph[dst(s_edge)])
                    push!(used_edges, real_s_edge)
                    true_cost += original_s_network_dir[src(real_s_edge), dst(real_s_edge)][:cost]
                end
            end
            edge_routing[v_edge] = order_path(original_s_network_dir, used_edges, node_placement[src(v_edge)], node_placement[dst(v_edge)]) 
        end
        mapping = Mapping(v_network, original_s_network_dir, node_placement, edge_routing)
        reduced_cost = objective_value(model) - penalty_routing
        return (
            sub_mapping=mapping, 
            real_cost=true_cost, 
            reduced_cost=reduced_cost
        )
    else
        return (
            sub_mapping=nothing, 
            real_cost=10e9, 
            reduced_cost=10e9
        )
    end


end






