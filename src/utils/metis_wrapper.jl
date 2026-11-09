using Base.Libc

# Path to the KaHIP shared library
using Libdl
using METIS_jll
Libdl.dlopen(METIS_jll.libmetis_path)
#using CEnum

#const idx_t = Int32
#const real_t = Cfloat

using Graphs





function graph_to_csr(g::Graph)
    n = nv(g)  # Number of vertices

    # Build xadj and adjncy arrays
    xadj = Int[]
    adjncy = Int[]  # This will hold all adjacency lists concatenated

    push!(xadj, 0)
    for i in 1:n
        append!(adjncy, neighbors(g, i))
        push!(xadj, length(adjncy))
    end

    return xadj, adjncy
end

function partition_metis(g, nb_clusters, imbalance)
    xadj, adjncy = graph_to_csr(g)

    # adapt it for c, because they start indexing at 0 lol
    for i in 1:length(adjncy)
        adjncy[i] = adjncy[i] - 1
    end

    edgecut, part = simplified_kway(xadj, adjncy, nb_clusters, imbalance)

    # adapt it for julia, because we start indexing at 1 lol
    for i in 1:length(part)
        part[i] = part[i] + 1
    end

    #write_some_stuff(part, nb_clusters)

    return part
end


function simplified_kway(xadj, adjncy, nparts, imbalance)
    n = Int32(length(xadj)-1)  # Number of vertices
    
    nparts = Int32(nparts)
    imbalance = Float64(imbalance)

    xadj = Vector{Int32}(xadj)
    adjncy = Vector{Int32}(adjncy)

    imbalance_vec = Vector{Float32}()
    push!(imbalance_vec, imbalance)

    status, objval, part = metis_kway(n, xadj, adjncy, nparts, nothing, nothing, nothing, nothing, imbalance_vec, nothing)
    #println("Partition: $part")

    return(objval, part)

end


function metis_kway(n::Int32,
    xadj::Vector{Int32},
    adjncy::Vector{Int32},
    nparts::Int32,
    vwgt::Union{Vector{Int32},Nothing}=nothing, # weight of vertices, useless
    vsize::Union{Vector{Int32},Nothing}=nothing, # size of the vertices for computing the total communication volume..?
    adjwgt::Union{Vector{Int32},Nothing}=nothing, #The weights of the edges 
    tpwgts::Union{Vector{Float32},Nothing}=nothing, # desired weight for each partition and constraint
    ubvec::Union{Vector{Float32},Nothing}=nothing, # load imbalance tolerance for each constraint. Gotta modify this one!
    options::Union{Vector{Int32},Nothing}=nothing) # Some options. Gotta use METIS_OPTION_CONTIG...

    ncon = Ref(Int32(1)) # From metis : The number of balancing constraints. It should be at least 1
    objval = Ref(Int32(0)) # not sure
    part = Vector{Int32}(undef, n) # The partition



    
    METIS_NOPTIONS = 40
    options = Vector{Int32}(undef, METIS_NOPTIONS)

    ccall(
        (:METIS_SetDefaultOptions, METIS_jll.libmetis),
        Cvoid,
        (Ptr{Int32},),
        options
    )

    
    
    METIS_OPTION_CONTIG = 11
    options[METIS_OPTION_CONTIG+1] = 1
    
    METIS_OPTION_SEED = 12
    options[METIS_OPTION_SEED+1] = 9323


    METIS_OPTION_UFACTOR = 16
    options[METIS_OPTION_UFACTOR+1] = 1 # this has no effect in my tests, remove it ?

    
    # Convert Nothing to C_NULL for optional arguments
    function to_cptr(x)
        x === nothing ? C_NULL : x
    end


    status = ccall(
        (:METIS_PartGraphKway, METIS_jll.libmetis),
        Int32,
        (Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
        Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
        Ptr{Float32}, Ptr{Float32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}),
        Ref(n), ncon, xadj, adjncy,
        to_cptr(vwgt), to_cptr(vsize), to_cptr(adjwgt),
        Ref(nparts), to_cptr(tpwgts), to_cptr(ubvec),
        to_cptr(options), objval, part
    )

    return status, objval[], part
end



function write_some_stuff(part, nb_part)

    trucs = [0 for i in 1:nb_part]

    for node in 1:length(part)
        trucs[part[node]] += 1
    end

    println("The partition is:")
    for t in 1:nb_part
        println("For subraph $t, $(trucs[t]) nodes")
    end
end