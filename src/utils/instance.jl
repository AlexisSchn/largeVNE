using Graphs, MetaGraphsNext

### Instance ##
struct Instance
    v_network
    s_network
    s_network_dir
end


function Instance(v_network, s_network)

    s_network_dir = generate_dir_sn(s_network)

    return Instance(v_network, s_network, s_network_dir)
end


function generate_dir_sn(s_network)

    s_network_dir = MetaGraph(
        DiGraph(),
        Int,
        Dict,
        Dict,
        Dict(:name => s_network[][:name] * "_dir", :type => "Substrate", :directed => true)
    )    

    for s_node in vertices(s_network)
        add_vertex!(s_network_dir, s_node, s_network[s_node])
    end

    for s_edge in edges(s_network)
        add_edge!(s_network_dir, src(s_edge), dst(s_edge), copy(s_network[src(s_edge), dst(s_edge)]))
        s_network_dir[src(s_edge), dst(s_edge)][:cost] = s_network[src(s_edge), dst(s_edge)][:cost]
        add_edge!(s_network_dir, dst(s_edge), src(s_edge), copy(s_network[src(s_edge), dst(s_edge)]))
        s_network_dir[dst(s_edge), src(s_edge)][:cost] = s_network[src(s_edge), dst(s_edge)][:cost]
    end

    return s_network_dir

end
