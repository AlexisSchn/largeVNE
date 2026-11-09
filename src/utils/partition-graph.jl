 
includet("import_utils.jl")
includet("kahip_wrapper.jl")
includet("metis_wrapper.jl")


using Statistics




function partition_graph(graph, nb_clusters; max_umbalance=1.25)

    # FIRST, do it with Kahip... 
    clusters = partition_graph_kahip(graph, nb_clusters, inbalance=max_umbalance-1)

    moyenne = mean([length(cluster) for cluster in clusters])
    current_imb = maximum([length(cluster) / moyenne for cluster in clusters] ∪ [moyenne / length(cluster) for cluster in clusters] )

    if current_imb < max_umbalance
        println("Best partition found has imbalance of $current_imb.")
        return clusters
    end

    println("KAHIP isnt balanced enough, using Metis...")
    
    # If very poorly balanced, do it with METIS !
    # Since connectivity is enforced, sometime, it will not the best
    
    
    keep_on=true
    clusters = []
    i_iter = 0
    imbalanced_allowed = max_umbalance
    while keep_on
        partition = partition_metis(graph, nb_clusters, imbalanced_allowed)

        clusters = [Vector{Int64}() for i in 1:nb_clusters]
        for s_node in vertices(graph)
            push!(clusters[partition[s_node]], s_node)
        end

        moyenne = mean([length(cluster) for cluster in clusters])
        current_imb = maximum([length(cluster) / moyenne for cluster in clusters] ∪ [moyenne / length(cluster) for cluster in clusters] )

        #println("Right now, current imb: $current_imb... Should be less than $max_umbalance right?")
        if current_imb < max_umbalance
            println("Best partition found has imbalance of $current_imb.")
            return clusters
        end

        imbalanced_allowed = (1+imbalanced_allowed)/2

        i_iter += 1
        if i_iter > 20
            println("Well, couldn't find a good, balanced mapping... Current one: $current_imb")
            return clusters
        end
    end

end



function partition_graph_kahip(graph, nb_clusters; inbalance = 0.1)

    
    # 1 : Partitionner
    partition = partition_kahip(graph, nb_clusters, inbalance)
    clusters = [Vector{Int64}() for i in 1:nb_clusters]
    for s_node in vertices(graph)
        push!(clusters[partition[s_node]], s_node)
    end

    # 2 : Corriger
    for cluster in clusters
        simple_subgraph, vmap = induced_subgraph(graph, cluster)
        if !is_connected(simple_subgraph)
            #print("Issue with unconnected subgraphs to correct...")
            components = connected_components(simple_subgraph)
            component_sorted = sort(components, by=x->length(x), rev=true)
            for subcluster in component_sorted[2:length(component_sorted)]
                nodes_original = [vmap[node] for node in subcluster]
                subgraph_neighbors = zeros(Int, nb_clusters)
                for node in nodes_original
                    for neighbor in neighbors(graph, node)
                        if neighbor ∉ cluster
                            subgraph_neighbors[partition[neighbor]] += 1
                        end
                    end
                end
                most_connected_subgraph = sortperm(subgraph_neighbors, rev=true)
                cluster_to_put_nodes_in = most_connected_subgraph[1]
                append!(clusters[cluster_to_put_nodes_in], nodes_original)
                for node in nodes_original
                    partition[node] = cluster_to_put_nodes_in
                end
                filter!(e->e∉nodes_original, cluster)
            end
        end
    end
    
    
    #= equilibrating virtual subgraphs, to get all of them pretty close...
    target_size = nv(instance.v_network) / nb_clusters
    max_size = maximum(length(cluster) for cluster in clusters)
    min_size = minimum(length(cluster) for cluster in clusters)
    iter = 1
    while max_size > 1.2 * target_size || min_size < 0.8 * target_size

        if max_size > 1.2 * target_size # let's remove a node from the biggest subgraph
            oversized_clusters = findall(c -> length(c) > 1.2 * target_size, clusters)
            cluster_to_reduce = clusters[oversized_clusters[1]]

            neighboring_clusters = Dict()
            for node in cluster_to_reduce
                for neighbor in neighbors(graph, node)
                    if neighbor ∉ cluster_to_reduce
                        if partition[neighbor] ∈ keys(neighboring_clusters)
                            push!(neighboring_clusters[partition[neighbor]], node)
                        else
                            neighboring_clusters[partition[neighbor]] = [node]
                        end
                    end
                end
            end
            #println("hum sooo... $neighboring_clusters")
            smallest_cluster_around = argmin(i -> length(clusters[i]), keys(neighboring_clusters))
            #println("Well, let's add it to $smallest_cluster_around !")
            node_to_change_cluster = neighboring_clusters[smallest_cluster_around][1]
            deleteat!(cluster_to_reduce, findfirst(==(node_to_change_cluster), cluster_to_reduce)) 
            push!(clusters[smallest_cluster_around], node_to_change_cluster)
            partition[node_to_change_cluster] = smallest_cluster_around
            # put a node in the neighbor cluster with fewer nodes



        else # let's add a node to the smallest subgraph


            # do stuff..


        end




        # Update max size
        max_size = maximum(length(cluster) for cluster in clusters)
        min_size = minimum(length(cluster) for cluster in clusters)

        if iter == 1
            println("Have to do extra balancing of the subgraphs, due to KaHIP...")
        end
        iter += 1
        if iter > 30
            println("Coudn't balance the subgraphs well... Sorry.")
            break
        end
    end
    =#

    return clusters

end




function partition_vn_metis(instance, nb_clusters)

    graph = instance.v_network.graph
    
    # Partitionning. Since connectivity is enforced, sometime, it will not the best
    println("$nb_clusters clusters to do... Partitionning with METIS!")
    best_clusters = nothing
    best_imb = 10000
    imb = [1.01, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3]
    #imb = [1.05]
    for imbalance in imb
        partition = partition_metis(graph, nb_clusters, imbalance)

        clusters = [Vector{Int64}() for i in 1:nb_clusters]
        for s_node in vertices(graph)
            push!(clusters[partition[s_node]], s_node)
        end

        moyenne = mean([length(cluster) for cluster in clusters])
        current_imb = maximum([length(cluster) / moyenne for cluster in clusters])
        if current_imb < 1.10
            best_clusters = clusters
            best_imb = current_imb
            break
        end
        if current_imb < best_imb
            best_imb = current_imb
            best_clusters = clusters
        end
    end
    
    println("Best partition found has imbalance of $best_imb.")

    return best_clusters
end





