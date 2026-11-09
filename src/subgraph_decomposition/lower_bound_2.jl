using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX
using OrderedCollections
using Printf


#general
includet("../utils/import_utils.jl")

# utils colge
includet("utils/master_problem.jl")
includet("utils/graph_decomposition.jl")
includet("../utils/partition-graph.jl")

# pricers
includet("pricers/greedy_pricer_subsn.jl")
includet("pricers/milp_pricer_subsn.jl")

includet("pricers/greedy_pricer.jl")
includet("pricers/milp_pricer.jl")


function lower_bound(instance; nb_virtual_subgraph=0, alpha_colge=0.9)
    
    
    println("Starting...")  

    time_master = 0
    nb_iter=0
    lower_bound = 0
    nb_columns=0
    time_beginning = time()

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir


    # === SOME PARAMETERS === #
    nb_columns_greedy = 1500
    time_max_heuristics = 600
    time_max_overall = 3600


    # ======= SETTING UP THE DECOMPOSITION ======= #
    if nb_virtual_subgraph == 0
        nb_virtual_subgraph = floor(Int, nv(v_network.graph)/10)
    end
        v_node_partitionning = partition_graph(v_network.graph, nb_virtual_subgraph, max_umbalance=1.2)
    vn_decompo = set_up_decompo(instance, v_node_partitionning)
    vn_subgraphs = vn_decompo.subgraphs

    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_network, vn_subgraphs)
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges")





    # === GETTING READY === #

    # master problem things
    master_problem = set_up_master_problem(instance, vn_decompo)
    model = master_problem.model
    print("Master problem set... ")
    empty_dual_costs = get_empty_duals(instance, vn_decompo)


    # partition things
    size_max_v_subgraph = maximum(nv(v_subgraph.graph) for v_subgraph in vn_decompo.subgraphs)
    nb_substrate_nodes_capacited = 0
    for s_node in vertices(s_network)
        if get_attribute_node(s_network, s_node, :cap) >= 1
            nb_substrate_nodes_capacited += 1
        end
    end
    ratio_capacited = nb_substrate_nodes_capacited / nv(s_network)
    nb_substrate_subgraphs = maximum([floor(Int, nv(s_network) / (size_max_v_subgraph*1.5/ratio_capacited)), nb_virtual_subgraph])
    clusters = partition_graph(s_network.graph, nb_substrate_subgraphs; max_umbalance = 1.3)
    sn_subgraphs = []
    for (i_subgraph, cluster) in enumerate(clusters)
        induced_subg = my_induced_subgraph(s_network, cluster, "sub_sn_$i_subgraph")
        push!(sn_subgraphs, Subgraph(induced_subg, cluster))
    end




    # ====== STEP 1: PAVING WITH HEURISTICS (local search then milp)
    println("------- Part 1: Reduced pricers")

    # need to complete the substrate graphs!
    nb_nodes_subgraph = 3. * size_max_v_subgraph
    sn_clusters_2 = complete_clusters(clusters, s_network, nb_nodes_to_have = nb_nodes_subgraph)
    sn_subgraphs_2 = []
    for (i_subgraph, cluster) in enumerate(sn_clusters_2)
        induced_subg = my_induced_subgraph(s_network, cluster, "sub_sn_$i_subgraph")
        push!(sn_subgraphs_2, Subgraph(induced_subg, cluster))
    end
    println("Adding nodes to the previous substrate decomposition done:")
    print_stuff_subgraphs(s_network, sn_subgraphs_2)



    # ======= ITERATION OF MASTER PROBLEM ======= #
    optimize!(model)
    time_master +=  solve_time(model)

    cg_value = objective_value(model)
    time_overall = time() - time_beginning
    @printf("Iter %2d  RMP value: %10.3f  %5d column    time: %5.2fs  \n",
        nb_iter, cg_value, nb_columns, time_overall
    )
    

    used_dual_costs = empty_dual_costs

    keep_on = true
    reason = "I don't know"
    pricer = "greedy"

    current_alpha_colge = 0.

    while keep_on
        nb_iter += 1


        # ---- pricers part
        old_dual_costs = used_dual_costs
        current_dual_costs = get_duals(instance, vn_decompo, master_problem)
        used_dual_costs = average_dual_costs(instance, vn_decompo, old_dual_costs, current_dual_costs, alpha=current_alpha_colge)
        

        sum_pricers_values = 0

        for vn_subgraph in vn_decompo.subgraphs

            time_overall = time() - time_beginning
            time_limit_pricer = time_max_overall - time_overall
            if time_limit_pricer < 0.1
                continue
            end

            result = solve_greedy_pricer(vn_subgraph, instance, vn_decompo, used_dual_costs; nb_iterations=500)

            sub_mapping = result[:sub_mapping]
            true_cost = result[:real_cost]
            reduced_cost = result[:reduced_cost]
            
            if (!isnothing(sub_mapping)) && reduced_cost < -0.0001
                has_found_new_column = true
                column = add_column(master_problem, instance, vn_subgraph, sub_mapping, true_cost)
                nb_columns += 1
            end

            sum_pricers_values += reduced_cost

        end
        
        average_reduced_costs = sum_pricers_values / nb_virtual_subgraph

        # ---- master problem part

        optimize!(model)
        time_master +=  solve_time(model)       
        #println("Time last iter master : $(solve_time(model))")
        cg_value = objective_value(model)

        time_overall = time()-time_beginning


        @printf("Iter %2d  CG bound: %10.3f,  lower bound: %10.3f,    nb columns  %5d,       time: %5.2fs,   average reduced costs %5.2f \n",
                    nb_iter, cg_value, lower_bound, nb_columns, time_overall, average_reduced_costs)
    


        
        # ----- useful things

        if cg_value < 10e4 && average_reduced_costs > - 25 && current_alpha_colge < 0.1
            println("Stabilizing now!")
            current_alpha_colge = alpha_colge
        end

        if average_reduced_costs > -1. && nb_columns > nb_columns_greedy * 0.7
            keep_on = false
            reason="Columns are not that interesting!"
        end


        if nb_columns > nb_columns_greedy
            keep_on = false
            reason="Got enough columns"
        end



        if time_overall > time_max_heuristics
            keep_on = false
            reason="Enough time"
        end



    end
    println("\n\n Step 2 finished, reason: $reason. \n\n\n")
    

    
    
    # ====== STEP 3: full pricers
    println("\n------- Solving method: Exact pricers")

    pricers_full = Dict()
    for subgraph in vn_decompo.subgraphs
        pricers_full[subgraph] = set_up_pricer(instance, subgraph)
    end
    old_cg_value = cg_value
    keep_on = true
    reason = "I don't know"
    while keep_on
        nb_iter += 1

        # ---- Pricers things
        current_dual_costs = get_duals(instance, vn_decompo, master_problem)
        old_dual_costs = used_dual_costs
        used_dual_costs = average_dual_costs(instance, vn_decompo, old_dual_costs, current_dual_costs, alpha=alpha_colge)
        

        has_found_new_column = false
        sum_pricers_values = 0

        for vn_subgraph in vn_decompo.subgraphs
            pricer = pricers_full[vn_subgraph]

            time_overall = time() - time_beginning
            time_limit_pricer = time_max_overall - time_overall
            if time_limit_pricer < 0.1
                continue
            end

            sub_mapping, true_cost, reduced_cost = update_solve_pricer(instance, vn_decompo, pricer, used_dual_costs; time_limit = time_limit_pricer)

            
            if (!isnothing(sub_mapping)) && reduced_cost < -0.0001
                has_found_new_column = true
                column = add_column(master_problem, instance, vn_subgraph, sub_mapping, true_cost)
                nb_columns += 1
            end
            
            if isnothing(sub_mapping)
                println("Pricer with no solution found, stopping the CG...")
                reason="pricer-unfeasible"
                has_found_new_column = false
                break
            end

            sum_pricers_values += reduced_cost

        end

        time_left = time_max_overall - (time() - time_beginning)
        current_lower_bound = (1-alpha_colge) * cg_value + (alpha_colge) * old_cg_value + sum_pricers_values
        if current_lower_bound > lower_bound && time_left > 0.5 # If there is no time left, it means that the pricer was not exact. (well, I could also do this by checking if it was optimal, but i'm too lazy.)
            lower_bound = current_lower_bound
        end


        # ----- Master problem stuff

        optimize!(model)
        time_master +=  solve_time(model)

        old_cg_value = cg_value
        cg_value = objective_value(model)

        time_overall = time()-time_beginning
        average_obj = sum_pricers_values/length(vn_decompo.subgraphs)


        @printf("Iter %2d  CG bound: %10.3f,  lower bound: %10.3f,    nb columns  %5d,       time: %5.2fs,   average reduced costs %5.2f \n",
                    nb_iter, cg_value, lower_bound, nb_columns, time_overall, average_obj)



        if !has_found_new_column
            keep_on = false
            reason="no improving columns"
        end


        if time_left < 0.5
            keep_on = false
            reason="time limit"
        end

        if cg_value < 10e4 && average_obj > - 50
            current_alpha_colge = alpha_colge
        end

    end

    print("\n==================== CG finished ====================\nReason: $reason \n")
    println("Time in MP: $(round(time_master; digits=3)) , time in SP: idk, time overall: $(round(time_overall; digits=3))")
    println("$nb_iter iters, final value: $(round(cg_value; digits=3))")
    println("====================================================\n")



    
    
    # ======= GETTING A SOLUTION ======= #
    
    time_overall = time() - time_beginning

    return (lower_bound = lower_bound,
            rmp_value=cg_value, 
            time_solving = time_overall, 
            nb_iter = nb_iter, 
            nb_columns = nb_columns)
end









function complete_clusters(original_clusters, s_network; nb_nodes_to_have = 30)


    # 2) adding nodes. For now, any adjacent nodes will do.
    i_cluster = 1

    final_clusters = []

    for original_cluster in original_clusters
        cluster = copy(original_cluster)
        all_neighbors = Dict()
        for s_node in cluster
            for neigh in neighbors(s_network, s_node)
                if neigh ∉ cluster
                    if neigh ∉ keys(all_neighbors)
                        all_neighbors[neigh] = 1
                    else
                        all_neighbors[neigh] += 1
                    end
                end
            end
        end

        added = []

        nb_nodes = length(cluster)

        while nb_nodes < nb_nodes_to_have

            # ranking the neighbors
            ranking = sort(collect(keys(all_neighbors)), by = x->all_neighbors[x], rev=true)
            # add the most connected neighbor
            new_s_node = ranking[1]
            nb_nodes+=1

            push!(cluster, new_s_node)
            push!(added, new_s_node)
            delete!(all_neighbors, new_s_node)

            for neigh in neighbors(s_network, new_s_node)
                if neigh ∉ cluster
                    if neigh ∉ keys(all_neighbors)
                        all_neighbors[neigh] = 1
                    else
                        all_neighbors[neigh] += 1
                    end
                end
            end
        end

        sub_s_network = my_induced_subgraph(s_network, cluster, "sub_sn_$i_cluster")
        i_cluster += 1

        push!(final_clusters, cluster)
    end

    return final_clusters

end

