# Large Virtual Network Embedding instances and solvers

This folder contains the instances and code used for the paper
> A. Benhamiche, P. Fouilhoux, L. Létocart, N. Perrot, and A. Schneider.
> **A Decomposition Approach for Large Virtual Network Embedding.**
> *Networks*, 2026. https://doi.org/10.1002/net.70062

The aim of this paper is to solve very large instances of VNE (up to 100 virtual nodes, 200 substrate nodes), using heuristics and mathematical programming methods. In particular, we introduce the Virtual Partition Formulation, a new integer formulation of VNE that can provide much better bound than the classical (compact) Flow Formulation.

The instances proposed are constructed as indicated in the paper, based on real networks from SNDlib and Internet Topology Zoo.


## Citation

If you use this code or these instances, please cite:

```bibtex
@article{benhamiche2026decomposition,
  author  = {Benhamiche, Amal and Fouilhoux, Pierre and L{\'e}tocart, Lucas and Perrot, Nancy and Schneider, Alexis},
  title   = {A Decomposition Approach for Large Virtual Network Embedding},
  journal = {Networks},
  year    = {2026},
  doi     = {10.1002/net.70062},
  url     = {https://doi.org/10.1002/net.70062}
}
```

## Dependencies

The code is written in Julia, using Graphs.jl, MetaGraphsNext.jl, JuMP.jl, CPLEX.jl and other packages.
