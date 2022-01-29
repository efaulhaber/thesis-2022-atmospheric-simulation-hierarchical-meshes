<!-- LTeX: language=en-US -->
# Discontinuous Galerkin Methods for Atmospheric Simulations on Hierarchical Meshes in Julia

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)

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

## List of all tables and corresponding commands

We provide a list of commands to reproduce every table in the thesis:

Table 5.1:
```julia
using Trixi
convergence_test("examples/structured_1d_dgsem/elixir_advection_basic.jl", 7)
```

Table 5.2:
```julia
using Trixi
convergence_test("examples/structured_1d_dgsem/elixir_euler_source_terms.jl", 4, polydeg=3)
```

Table 5.3:
```julia
using Trixi
convergence_test("examples/structured_1d_dgsem/elixir_euler_source_terms.jl", 4)
```

Table 5.4:
```julia
using Trixi
convergence_test("examples/structured_1d_dgsem/elixir_euler_source_terms.jl", 4, polydeg=3, surface_flux=flux_hll)
```

Table 5.5:
```julia
using Trixi
convergence_test("examples/structured_1d_dgsem/elixir_euler_source_terms.jl", 4, surface_flux=flux_hll)
```

Table 5.6:
```julia
using Trixi
convergence_test("examples/structured_2d_dgsem/elixir_advection_basic.jl", 5)
```

Table 5.7:
```julia
using Trixi
convergence_test("examples/structured_2d_dgsem/elixir_advection_waving_flag.jl", 5)
```

Table 5.8:
```julia
using Trixi
convergence_test("examples/structured_2d_dgsem/elixir_advection_basic_warped.jl", 5)
```

Table 5.9:
```julia
using Trixi
convergence_test("examples/structured_2d_dgsem/elixir_euler_source_terms.jl", 4)
```

Table 5.10:
```julia
using Trixi
convergence_test("examples/structured_2d_dgsem/elixir_euler_source_terms_rotated.jl", 4)
```

Table 5.11:
```julia
using Trixi
convergence_test("examples/structured_2d_dgsem/elixir_euler_source_terms_parallelogram.jl", 4)
```

Table 5.12:
```julia
using Trixi
convergence_test("examples/structured_2d_dgsem/elixir_euler_source_terms_waving_flag.jl", 4)
```

Table 5.13:
```julia
using Trixi
trixi_include("examples/structured_2d_dgsem/elixir_euler_free_stream.jl")
trixi_include("examples/structured_2d_dgsem/elixir_euler_free_stream.jl", polydeg=4, cfl=1.4)
```

Table 5.14:
```julia
using Trixi
convergence_test("examples/structured_3d_dgsem/elixir_advection_basic.jl", 4)
```

Table 5.15:
```julia
using Trixi
convergence_test("examples/structured_3d_dgsem/elixir_advection_nonperiodic_curved.jl", 4)
```

Table 5.16:
```julia
using Trixi
convergence_test("examples/structured_3d_dgsem/elixir_euler_source_terms.jl", 3)
```

Table 5.17:
```julia
using Trixi
convergence_test("examples/structured_3d_dgsem/elixir_euler_source_terms_nonperiodic_curved.jl", 3)
```

Table 5.18:
```julia
using Trixi
trixi_include("examples/structured_3d_dgsem/elixir_euler_free_stream.jl")
```

Table 5.19:
```julia
using Trixi
convergence_test("examples/structured_2d_dgsem/elixir_euler_source_terms_ring_coupled.jl", 4)
```

Table 5.20:
```julia
using Trixi
convergence_test("examples/structured_3d_dgsem/elixir_euler_source_terms_cubed_sphere_coupled.jl", 3)
```

Table 5.21:
```julia
using Trixi
convergence_test("examples/p4est_2d_dgsem/elixir_advection_nonconforming_flag.jl", 4)
```

Table 5.22:
```julia
using Trixi
convergence_test("examples/p4est_2d_dgsem/elixir_advection_unstructured_flag.jl", 4)
```

Table 5.23:
```julia
using Trixi
convergence_test("examples/p4est_2d_dgsem/elixir_advection_amr_unstructured_flag.jl", 3, amr=false)
trixi_include("examples/p4est_2d_dgsem/elixir_advection_amr_unstructured_flag.jl")
```

Table 5.24:
```julia
using Trixi
convergence_test("examples/p4est_2d_dgsem/elixir_advection_amr_solution_independent.jl", 3)
```

Table 5.25:
```julia
using Trixi
convergence_test("examples/p4est_2d_dgsem/elixir_euler_source_terms_nonconforming_unstructured_flag.jl", 4)
```

Table 5.26:
```julia
using Trixi
trixi_include("examples/p4est_2d_dgsem/elixir_euler_free_stream.jl")
```

Table 5.27:
```julia
using Trixi
convergence_test("examples/p4est_3d_dgsem/elixir_advection_nonconforming.jl", 4)
```

Table 5.28:
```julia
using Trixi
convergence_test("examples/p4est_3d_dgsem/elixir_advection_unstructured_curved.jl", 4)
```

Table 5.29:
```julia
using Trixi
convergence_test("examples/p4est_3d_dgsem/elixir_advection_amr_unstructured_curved.jl", 3, amr=false)
trixi_include("examples/p4est_3d_dgsem/elixir_advection_amr_unstructured_curved.jl")
```

Table 5.30:
```julia
using Trixi
convergence_test("examples/p4est_3d_dgsem/elixir_euler_source_terms_nonconforming_unstructured_curved.jl", 3)
```

Table 5.31:
```julia
using Trixi
trixi_include("examples/p4est_3d_dgsem/elixir_euler_free_stream.jl")
```

Table 5.32:
```julia
using Trixi
trixi_include("examples/p4est_3d_dgsem/elixir_euler_free_stream_extruded.jl")
```

Table 6.1:
```julia
using Trixi
convergence_test("examples/p4est_3d_dgsem/elixir_euler_circular_wind_nonconforming.jl", 5)
```

Table 6.2:
```julia
using Trixi
convergence_test("examples/p4est_3d_dgsem/elixir_euler_box_gravity.jl", 3, conforming=true)
convergence_test("examples/p4est_3d_dgsem/elixir_euler_box_gravity.jl", 3)
```

## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
