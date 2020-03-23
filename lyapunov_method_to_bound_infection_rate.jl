# Install required dependies
using Pkg
Pkg.activate(".")
Pkg.instantiate()


include("covid-19-model.jl")

##############################################
# Lyapunov method for upper bounding the number 
# of infected people
##############################################

function lyapunnov_upper_bound(;r, r₀,
        deg_V, deg_putinar, 
        solver) 
    
    model = SOSModel(solver)

    # Define the state variable and the vector field
    @polyvar S E Ir Iu t
    vars = [S, E, Ir, Iu]
    tvars = [S, E, Ir, Iu, t]
    f = transmission_vf(S, E, Ir, Iu)

    # Omega
    space = [ [r₀ - v for v=vars]..., 
             [r₀ + v for v=vars]...]
    
    # X0
    safe = [ [r - v for v=vars]..., 
             [r + v for v=vars]...]        

    
    ϕ = Ir + Iu
    

    # define the Putinar multipliers

    mons = monomials(tvars, 0:deg_putinar)
    p_Vdot = @variable model [1:size(space,1)] Poly(mons)
    p_ϕ_V = @variable model [1:size(space,1)] Poly(mons)
    p_V_max = @variable model [1:size(space,1)] Poly(mons)
    
    # V
    @variable model V Poly(monomials(tvars, 0:deg_V))
    @variable model δ
    
    ∇V = differentiate(V, vars)
    V_t = differentiate(V, t)
    Vdot = V_t + ∇V' * f
    

    @constraint model -Vdot- space' * p_Vdot >= 0
    @constraint model -(ϕ - V) - space' * p_ϕ_V >= 0
    @constraint model δ - subs(V, t=>0)  - safe' * p_V_max >= 0

    @constraint model p_Vdot .>= 0
    @constraint model p_ϕ_V .>= 0
    @constraint model p_V_max .>= 0

    @objective model Min δ 
    
    # Solve the model
    optimize!(model)

    # Returns true if the model is successfully solved
    st = termination_status(model)   
    @show (st)
    
    return value(δ)

end

##############################################
# Find the optimal δ corresponding to different 
# values of deg(V) and deg(putinar).
##############################################

if ~isinteractive()
  println("Computing the table of δs...")

  # deg(B) and deg putinar to try
  deg_Vs = 2:6
  deg_putinars = 2:2:6
  # 2D array to save the result
  δs = zeros(size(deg_Vs,1), size(deg_putinars,1))

  for (i, deg_V) in enumerate(deg_Vs)
      for (j, deg_putinar) in enumerate(deg_putinars)

          δ = lyapunnov_upper_bound(;deg_V=deg_V,
                                  deg_putinar=deg_putinar,
                                  r₀=r₀, r=r, solver=solver)

          δs[i, j] = δ

          @show deg_V, deg_putinar, δ

      end
  end

  @show round.(δs, digits=3)

end
