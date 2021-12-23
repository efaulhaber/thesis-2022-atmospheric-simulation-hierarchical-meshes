
using OrdinaryDiffEq
using Trixi

###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleEulerEquations2D(1.4)

initial_condition = initial_condition_convergence_test

source_terms = source_terms_convergence_test

# BCs must be passed as Dict
boundary_condition = BoundaryConditionDirichlet(initial_condition)
boundary_conditions = Dict(
  :all => boundary_condition
)

solver = DGSEM(polydeg=4, surface_flux=flux_hll)

# Deformed rectangle that looks like a waving flag,
# lower and upper faces are sinus curves, left and right are vertical lines.
f1(s) = SVector(-1.0, s - 1.0)
f2(s) = SVector( 1.0, s + 1.0)
f3(s) = SVector(s, -1.0 + sin(0.5 * pi * s))
f4(s) = SVector(s,  1.0 + sin(0.5 * pi * s))
faces = (f1, f2, f3, f4)

Trixi.validate_faces(faces)
mapping_flag = Trixi.transfinite_mapping(faces)

# Get the uncurved mesh from a file
# Unstructured mesh with 24 cells of the square domain [-1, 1]^n
mesh_file = joinpath(@__DIR__, "square_unstructured_2.inp")

mesh = P4estMesh{2}(mesh_file, polydeg=3,
                    mapping=mapping_flag)

# Refine bottom left quadrant of each tree to level 3
function refine_fn(p4est, which_tree, quadrant)
  if quadrant.x == 0 && quadrant.y == 0 && quadrant.level < 3
    # return true (refine)
    return Cint(1)
  else
    # return false (don't refine)
    return Cint(0)
  end
end

# Refine recursively until each bottom left quadrant of a tree has level 2
# The mesh will be rebalanced before the simulation starts
refine_fn_c = @cfunction(refine_fn, Cint, (Ptr{Trixi.p4est_t}, Ptr{Trixi.p4est_topidx_t}, Ptr{Trixi.p4est_quadrant_t}))
Trixi.refine_p4est!(mesh.p4est, true, refine_fn_c, C_NULL)
Trixi.balance!(mesh)

# Refine everything again for convergence tests
refine_fn2(p4est, which_tree, quadrant) = Cint(1)
refine_fn_c2 = @cfunction(refine_fn2, Cint, (Ptr{Trixi.p4est_t}, Ptr{Trixi.p4est_topidx_t}, Ptr{Trixi.p4est_quadrant_t}))

# Abuse the variable `initial_refinement_level` for this because this will be increased
# during convergence tests
initial_refinement_level = 0
for i in 1:initial_refinement_level
  Trixi.refine_p4est!(mesh.p4est, false, refine_fn_c2, C_NULL)
end

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms=source_terms,
                                    boundary_conditions=boundary_conditions)


###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 2.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval=analysis_interval)

alive_callback = AliveCallback(analysis_interval=analysis_interval)

save_restart = SaveRestartCallback(interval=1000,
                                   save_final_restart=true)

save_solution = SaveSolutionCallback(interval=1000,
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     solution_variables=cons2prim)

stepsize_callback = StepsizeCallback(cfl=0.8)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        save_restart, save_solution,
                        stepsize_callback)
###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
            dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep=false, callback=callbacks);
summary_callback() # print the timer summary
