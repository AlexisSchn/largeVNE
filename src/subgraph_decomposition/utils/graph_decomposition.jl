using Graphs, MetaGraphsNext

struct Subgraph
    graph
    nodes_of_main_graph
end



function print_stuff_subgraphs(original_graph, subgraphs)

    println("For $original_graph, there is $(length(subgraphs)) subgraphs:")
    for subgraph in subgraphs
        println("       $(subgraph.graph[][:name]) with $(nv(subgraph.graph)) nodes and $(ne(subgraph.graph)) edges")
    end
    
end


struct NetworkDecomposition
    subgraphs
    v_nodes_assignment
    v_edges_master
end

function set_up_decompo(instance, node_partitionning)

    vn = instance.v_network

        
    node_assignment = Dict()
    for v_node in vertices(vn)
        node_assignment[v_node] = Dict()
    end

    # getting the subgraphs and the node assignment
    # i couldnt make the base induced_graph function work so I did adapt it
    subgraphs = []
    for (i_subgraph, v_nodes) in enumerate(node_partitionning)
        subgraph = Subgraph(my_induced_subgraph(vn, v_nodes, "subgraph_$i_subgraph"), v_nodes)
        
        for (i_node, v_node) in enumerate(v_nodes)
            node_assignment[v_node][subgraph] = i_node
        end
        push!(subgraphs, subgraph)
        #println("Look at my nice graph for the nodes $v_nodes")
        #print_graph(subgraph.graph)
    end


    # finding out the master virtual edges
    v_edge_master = [] 
    for v_edge in edges(vn)
        in_master = true
        for subgraph_src in keys(node_assignment[src(v_edge)])
            for subgraph_dst in keys(node_assignment[dst(v_edge)])
                if subgraph_src == subgraph_dst
                    in_master = false
                end
            end
        end
        if in_master
            push!(v_edge_master, v_edge)
        end
    end

    

    vn_decompo = NetworkDecomposition(subgraphs, node_assignment, v_edge_master)



    return vn_decompo
end

