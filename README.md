<!-- LTeX: language=en-US -->
# Discontinuous Galerkin Methods for Atmospheric Simulations with p4est in Julia

This repository contains the code to reproduce the numerical results in my master's thesis.
The presented code uses [Trixi.jl](https://github.com/trixi-framework/Trixi.jl),
a software for adaptive numerical simulations of hyperbolic partial differential equations written in Julia.
To reproduce the numerical experiments using Trixi.jl, you need to install
[Julia](https://julialang.org/).
The numerical experiments were carried out using Julia v1.7.0.

Elixirs can be run with `trixi_include`, e.g.,
```julia
using Trixi
trixi_include("examples/structured_2d_dgsem/elixir_euler_free_stream.jl")
```
and variables in the elixir can be overwritten by passing keyword arguments, e.g.,
```julia
using Trixi
trixi_include("examples/structured_2d_dgsem/elixir_euler_free_stream.jl", polydeg=4, cfl=1.4)
```

Simulation results can be converted to the VTK format with
[Trixi2Vtk.jl](https://github.com/trixi-framework/Trixi2Vtk.jl)
and visualized with [ParaView](https://www.paraview.org/).

Convergence tests can be conducted with `convergence_test`, e.g.,
```julia
using Trixi
convergence_test("examples/structured_1d_dgsem/elixir_euler_source_terms.jl", 4)
```
and variables in the elixir can be overwritten by passing keyword arguments, e.g.,
```julia
using Trixi
convergence_test("examples/structured_1d_dgsem/elixir_euler_source_terms.jl", 4, polydeg=3, surface_flux=flux_hll)
```


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
