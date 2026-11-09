

#using Revise, JuMP, CPLEX, Gurobi
using Revise, JuMP, CPLEX

includet("../utils/import_utils.jl")




# ========== CLASSICAL STUFF
function set_up_problem_ff_plus(instance, model)

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


    ### Variables
    @variable(model, x[vertices(v_network), vertices(s_network_dir)], binary=true);
    @variable(model, y[edges(v_network), edges(s_network_dir)], binary=true);



    ### Objective
    placement_cost = @expression(model, sum( node_demands[v_node] * node_costs[s_node] * x[v_node, s_node] 
                        for v_node in vertices(v_network) for s_node in vertices(instance.s_network) ))
    routing_cost = @expression(model, sum( edge_demands[src(v_edge), dst(v_edge)] * cost_edges[src(s_edge), dst(s_edge)]* y[v_edge, s_edge]
                        for v_edge in edges(v_network) for s_edge in edges(s_network_dir) ))
    @objective(model, Min, placement_cost + routing_cost);




    ###===== Pre treatement on capacities 
    for v_node in vertices(v_network)
        for s_node in vertices(s_network)
            if node_demands[v_node] > node_capacities[s_node]
                @constraint(model, x[v_node, s_node] == 0)
            end
        end
    end

    for v_edge in edges(v_network)
        for s_edge in edges(s_network)
            if edge_demands[src(v_edge), dst(v_edge)] > capacities_edges[src(s_edge), dst(s_edge)]
                @constraint(model, y[v_edge, s_edge] == 0)
                @constraint(model, y[v_edge, get_reverse_edge(s_network_dir, s_edge)] == 0) 
            end
        end
    end
    

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



    ##--- Edges 
    
    # edge capacity (undirected version !)
    for s_edge in edges(instance.s_network)
        @constraint(model, 
            sum( edge_demands[src(v_edge), dst(v_edge)] * (y[v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] + y[v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]  )
                for v_edge in edges(v_network)) 
            <= 
            capacities_edges[src(s_edge), dst(s_edge)] )
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
    
    
    
    #= Star capacity constraint => NOT FOR NOW?
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
    =#  
end



function solve_compact_ffplus(instance; time_solver = 30, stay_silent=true, solution_first = false, random_seed=false)
    
    v_network = instance.v_network
    s_network_dir = instance.s_network_dir

    time_start = time()
    model = Model(CPLEX.Optimizer)
    set_up_problem_ff_plus(instance, model)

    set_time_limit_sec(model, time_solver)
    if stay_silent
        set_silent(model)
    else
        print("Starting solving... ")
    end

    if solution_first
        set_optimizer_attribute(model, "CPXPARAM_Emphasis_MIP", 5)
    end

    if random_seed
        set_optimizer_attribute(model, "CPXPARAM_RandomSeed", rand(1:9999))
    end
    
    optimize!(model)

    status = primal_status(model)
    if status != MOI.FEASIBLE_POINT
        println("Infeasible or unfinished: $status")
        result = Dict()
        result["mapping"] = nothing
        result["sol_value"] = -1
        result["lower_bound"] = objective_bound(model)
        result["gap"]  = -1.
        result["node_count"] = node_count(model)
        result["time_solving"] = time() - time_start    
        return result
    end




    # Get the solution
    x_values = value.(model[:x])
    y_values = value.(model[:y])
    
    node_placement = []
    for v_node in vertices(v_network)
        for s_node in vertices(s_network_dir)
            if x_values[v_node, s_node] > 0.01
                push!(node_placement, s_node)
            end
        end
    end

    edge_routing = Dict()
    for v_edge in edges(v_network)
        edges_of_routing = []
        for s_edge in edges(s_network_dir)
            if y_values[v_edge, s_edge] > 0.01
                push!(edges_of_routing, s_edge)
            end
        end
        path = order_path(s_network_dir, edges_of_routing, node_placement[src(v_edge)], node_placement[dst(v_edge)])
        edge_routing[v_edge] = path
    end

    result = Dict()
    result["mapping"] = Mapping(v_network, s_network_dir, node_placement, edge_routing)
    result["sol_value"] = objective_value(model)
    result["lower_bound"] = objective_bound(model)
    result["gap"]  = relative_gap(model)
    result["node_count"] = node_count(model)
    result["time_solving"] = time() - time_start

    return result
end


function solve_compact_ffplus_linear(instance; time_solver = 1000, stay_silent=true)
    
    v_network = instance.v_network
    s_network_dir = instance.s_network_dir


    model = Model(CPLEX.Optimizer)
    set_up_problem_ff_plus(instance, model)

    set_time_limit_sec(model, time_solver)
    if stay_silent
        set_silent(model)
    else
        print("Starting solving... ")
    end

    relax_integrality(model)

    optimize!(model)

    status = primal_status(model)
    if status != MOI.FEASIBLE_POINT
        println("error! no solution possible...")
        return -999
    end

    # Get the solution
    x_values = value.(model[:x])
    y_values = value.(model[:y])


    println("\nNode placement:")
    for v_node in vertices(v_network)
        println("Of $v_node:")
        for s_node in vertices(s_network_dir)
            if x_values[v_node, s_node] > 0.01
                println("   on $s_node with $(x_values[v_node, s_node])")
            end
        end
    end

    println("\nEdge routing:")
    for v_edge in edges(v_network)
        println("Of $v_edge:")
        for s_edge in edges(s_network_dir)
            if y_values[v_edge, s_edge] > 0.01
                println("   on $s_edge with $(y_values[v_edge, s_edge])")
            end
        end
    end

    return (objective_value(model))
end
