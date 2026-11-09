

#using Revise, JuMP, CPLEX, Gurobi
using Revise, JuMP, CPLEX

includet("../utils/import_utils.jl")




# ========== CLASSICAL STUFF

function set_up_problem_ff_plusplus(instance, model, departure, continuity_degree, continuity_cap)

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
    if departure  
        for s_node in vertices(instance.s_network)
            for v_edge in edges(v_network)
                @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) 
                    >= x[src(v_edge), s_node])
            end
        end
    end

    nb_continuity = 0
    for s_edge_in in edges(s_network_dir)
        for v_edge in edges(v_network)
            s_node = dst(s_edge_in)
            
            if degree(s_network, s_node) < continuity_degree
                if continuity_cap  && (s_network[s_node][:cap]==0)
                    @constraint(model, sum(y[v_edge, s_edge_out] for s_edge_out in get_out_edges(s_network_dir, s_node) ) + x[dst(v_edge), s_node] 
                        >= y[v_edge, s_edge_in] + y[v_edge, get_reverse_edge(s_network_dir, s_edge_in)] )
                    nb_continuity += 1
                elseif !continuity_cap
                    @constraint(model, sum(y[v_edge, s_edge_out] for s_edge_out in get_out_edges(s_network_dir, s_node) ) + x[dst(v_edge), s_node] 
                        >= y[v_edge, s_edge_in] + y[v_edge, get_reverse_edge(s_network_dir, s_edge_in)] )
                    nb_continuity += 1
                end
            end
        end
    end
    
    #=
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
    =#
        
    return nb_continuity
end



function solve_compact_ffplusplus(instance; time_solver = 30, stay_silent=false, departure=true, continuity_degree = 3, continuity_cap=false)
    
    v_network = instance.v_network
    s_network_dir = instance.s_network_dir

    time_start = time()

    model = Model(CPLEX.Optimizer)

    # force single thread 
    #set_optimizer_attribute(model, "CPXPARAM_Threads", 1)  # or "CPX_PARAM_THREADS"

    # vary the internal randomness across runs
    set_optimizer_attribute(model, "CPXPARAM_RandomSeed", rand(1:10^9))  # or "CPX_PARAM_RANDOMSEED"

    nb_continuity = set_up_problem_ff_plusplus(instance, model, departure, continuity_degree, continuity_cap)

    set_time_limit_sec(model, time_solver)
    if stay_silent
        set_silent(model)
    else
        print("Starting solving... ")
    end



    optimize!(model)

    status = primal_status(model)
    if status != MOI.FEASIBLE_POINT
        println("Infeasible or unfinished: $status")
        return ( sol_value= 999999,
            lower_bound = 999999.,
            gap = -1.,
            node_count = 0,
            time_solving = (time() - time_start),
            nb_continuity = nb_continuity
        )
    end

    #=
    if !stay_silent

        x_values = value.(model[:x])
        y_values = value.(model[:y])
    
        println("Node placement:")
        for v_node in vertices(v_network)
            for s_node in vertices(s_network_dir)
                if x_values[v_node, s_node] > 0.01
                    println("$v_node is placed on $s_node")
                end
            end
        end
        println("\nEdge routing:")
        for v_edge in edges(v_network)
            print("Routing of $v_edge : ")
            for s_edge in edges(s_network_dir)
                if y_values[v_edge, s_edge] > 0.01
                    print(" $s_edge")
                end
            end
            print("\n")
        end
    end
    =#

    return ( sol_value= objective_value(model),
        lower_bound = objective_bound(model),
        gap = relative_gap(model),
        node_count = node_count(model),
        time_solving = (time() - time_start),
        nb_continuity = nb_continuity
    )
end



function solve_compact_ffplusplus_linear(instance; time_solver = 30, stay_silent=false, departure=true, continuity_degree = 3, continuity_cap=false)
    
    v_network = instance.v_network
    s_network_dir = instance.s_network_dir

    time_start = time()

    model = Model(CPLEX.Optimizer)


    nb_continuity = set_up_problem_ff_plusplus(instance, model, departure, continuity_degree, continuity_cap)

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
    return (objective_value(model))
end


