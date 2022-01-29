
using OrdinaryDiffEq
using Trixi
using LinearAlgebra

###############################################################################
# semidiscretization of the compressible Euler equations
gamma = 1.4
equations = CompressibleEulerEquations3D(gamma)


# Initial condition for an idealized baroclinic instability test
# https://doi.org/10.1002/qj.2241, Section 3.2 and Appendix A
function initial_condition_pressure_gradient(x, t, equations::CompressibleEulerEquations3D)
  r = x[3]

  radius_earth = 6.371229e6
  # Make sure that r is not smaller than radius_earth
  z = max(r - radius_earth, 0.0)

  # Unperturbed basic state
  rho, u, p = basic_state_baroclinic_instability_longitudinal_velocity(0.0, 0.0, z)

  # No velocity
  v1 = 0.0
  v2 = 0.0
  v3 = 0.0

  return prim2cons(SVector(rho, v1, v2, v3, p), equations)
end


# Unperturbed balanced steady-state.
# Returns primitive variables with only the velocity in longitudinal direction (rho, u, p).
# The other velocity components are zero.
function basic_state_baroclinic_instability_longitudinal_velocity(lon, lat, z)
  # Parameters from Table 1 in the paper
  # Corresponding names in the paper are commented
  radius_earth                   = 6.371229e6  # a
  half_width_parameter           = 2           # b
  gravitational_acceleration     = 9.80616     # g
  k                              = 3           # k
  surface_pressure               = 1e5         # p₀
  gas_constant                   = 287         # R
  surface_equatorial_temperature = 310.0       # T₀ᴱ
  surface_polar_temperature      = 240.0       # T₀ᴾ
  lapse_rate                     = 0.005       # Γ
  angular_velocity               = 7.29212e-5  # Ω

  # Distance to the center of the Earth
  r = z + radius_earth

  # In the paper: T₀
  temperature0 = 0.5 * (surface_equatorial_temperature + surface_polar_temperature)
  # In the paper: A, B, C, H
  const_a = 1 / lapse_rate
  const_b = (temperature0 - surface_polar_temperature) /
            (temperature0 * surface_polar_temperature)
  const_c = 0.5 * (k + 2) * (surface_equatorial_temperature - surface_polar_temperature) /
                            (surface_equatorial_temperature * surface_polar_temperature)
  const_h = gas_constant * temperature0 / gravitational_acceleration

  # In the paper: (r - a) / bH
  scaled_z = z / (half_width_parameter * const_h)

  # Temporary variables
  temp1 = exp(lapse_rate/temperature0 * z)
  temp2 = exp(-scaled_z^2)

  # In the paper: ̃τ₁, ̃τ₂
  tau1 = const_a * lapse_rate / temperature0 * temp1 + const_b * (1 - 2 * scaled_z^2) * temp2
  tau2 = const_c * (1 - 2 * scaled_z^2) * temp2

  # In the paper: ∫τ₁(r') dr', ∫τ₂(r') dr'
  inttau1 = const_a * (temp1 - 1) + const_b * z * temp2
  inttau2 = const_c * z * temp2

  # Temporary variables
  temp3 = r/radius_earth * cos(lat)
  temp4 = temp3^k - k/(k + 2) * temp3^(k+2)

  # In the paper: T
  temperature = 1 / ((r/radius_earth)^2 * (tau1 - tau2 * temp4))

  # In the paper: U, u (zonal wind, first component of spherical velocity)
  big_u = gravitational_acceleration/radius_earth * k * temperature * inttau2 * (temp3^(k-1) - temp3^(k+1))
  temp5 = radius_earth * cos(lat)
  u = -angular_velocity * temp5 + sqrt(angular_velocity^2 * temp5^2 + temp5 * big_u)

  # Hydrostatic pressure
  p = surface_pressure * exp(-gravitational_acceleration/gas_constant * (inttau1 - inttau2 * temp4))

  # Density (via ideal gas law)
  rho = p / (gas_constant * temperature)

  return rho, u, p
end


@inline function source_terms_gravity(u, x, t, equations::CompressibleEulerEquations3D)
  radius_earth               = 6.371229e6  # a
  gravitational_acceleration = 9.80616     # g
  angular_velocity           = 7.29212e-5  # Ω

  r = x[3]
  # Make sure that r is not smaller than radius_earth
  z = max(r - radius_earth, 0.0)
  r = z + radius_earth

  du1 = zero(eltype(u))

  # Gravity term
  temp = -gravitational_acceleration * radius_earth^2 / r^3
  du2 = zero(eltype(u))
  du3 = zero(eltype(u))
  du4 = temp * u[1] * x[3]
  du5 = temp * u[4] * x[3]

  return SVector(du1, du2, du3, du4, du5)
end


function indicator_test(u::AbstractArray{<:Any,5},
                        mesh, equations, dg::DGSEM, cache;
                        kwargs...)
  alpha = zeros(Int, nelements(dg, cache))
  radius_earth = 6.371229e6

  for element in eachelement(dg, cache)
    # Calculate coordinates at Gauss-Lobatto nodes
    for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
      x = Trixi.get_node_coords(cache.elements.node_coordinates, equations, dg, i, j, k, element)
      lambda = x[1] / radius_earth

      if pi/16 + 1e-3 < lambda < 3*pi/16 - 1e-3
        alpha[element] = 1
      end
    end
  end

  return alpha
end

function Trixi.get_element_variables!(element_variables, indicator::typeof(indicator_test), ::AMRCallback)
  return nothing
end

initial_condition = initial_condition_pressure_gradient

boundary_conditions = Dict(
  :z_neg  => boundary_condition_slip_wall,
  :z_pos  => boundary_condition_slip_wall
)

# This is a good estimate for the speed of sound in this example.
# Other values between 300 and 400 should work as well.
surface_flux = FluxLMARS(340)
volume_flux  = flux_kennedy_gruber
solver = DGSEM(polydeg=3, surface_flux=surface_flux, volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

radius_earth = 6.371229e6
circumference_earth = 2 * pi * radius_earth

trees_per_dimension = (4, 1, 4)
mesh = P4estMesh(trees_per_dimension, polydeg=1,
                 coordinates_min=(0.0, -1/64 * circumference_earth, radius_earth),
                 coordinates_max=(circumference_earth / 8, 1/64 * circumference_earth, radius_earth + 30000.0),
                 periodicity=(true, true, false), initial_refinement_level=0)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms=source_terms_gravity,
                                    boundary_conditions=boundary_conditions)


###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 5.0e4)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval=analysis_interval)

alive_callback = AliveCallback(analysis_interval=analysis_interval)

save_solution = SaveSolutionCallback(interval=1000,
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     solution_variables=cons2prim)

# Abuse the variable `initial_refinement_level` for this because this will be increased
# during convergence tests
initial_refinement_level = 0
amr_controller = ControllerThreeLevel(semi, indicator_test,
                                      base_level=initial_refinement_level,
                                      max_level=initial_refinement_level + 1, max_threshold=0.6)
amr_callback = AMRCallback(semi, amr_controller,
                           interval=0,
                           adapt_initial_condition=true,
                           adapt_initial_condition_only_refine=true)

# This variable can be changed by passing the keyword argument conforming=true
# to trixi_include or convergence_test
conforming=false
if conforming
  callbacks = CallbackSet(summary_callback,
                          analysis_callback,
                          save_solution,
                          alive_callback)
else
  callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        amr_callback,
                        save_solution,
                        alive_callback)
end


###############################################################################
# run the simulation

# Use a Runge-Kutta method with automatic (error based) time step size control
# Enable threading of the RK method for better performance on multiple threads
sol = solve(ode, RDPK3SpFSAL49(thread=OrdinaryDiffEq.True()), abstol=1.0e-6, reltol=1.0e-6,
            save_everystep=false, callback=callbacks);

summary_callback() # print the timer summary
