# Install required dependies
using Pkg
Pkg.activate(".")
Pkg.instantiate()


include("covid-19-model.jl")

##############################################
# Barrier method for upper bounding the number 
# of infected people
##############################################

function is_valid_upper_bound(;δ, r, r₀,
        deg_B, deg_putinar, 
        ϵ, solver) 
    # Returns true is δ can be certified to be an upper bound on `I_r
    # + I_u` using a barrier function of degree `deg_B` and putinar
    # certificates of degree `deg_putinar`.

    
    model = SOSModel(solver)

    # Define the state variable and the vector field
    @polyvar S E Ir Iu
    vars = [S, E, Ir, Iu]
    f = transmission_vf(S, E, Ir, Iu)


    # The defining inequalities of the whole space, the safe space,
    # and the unsafe space

    space = [ [r₀ - v for v=vars]..., 
             [r₀ + v for v=vars]...]
    
    safe = [ [r - v for v=vars]..., 
             [r + v for v=vars]...]        

    
    unsafe = [Ir + Iu - δ]
    

    # define the Putinar multipliers
    mons = monomials(vars, 0:deg_putinar)
    p_safe = @variable model [1:size(safe,1)] Poly(monomials(vars, 0:deg_putinar))
    p_unsafe = @variable model [1:size(unsafe,1)] Poly(monomials(vars, 0:deg_putinar))
    p_space = @variable model [1:size(space,1)] Poly(monomials(vars, 0:deg_putinar))

    # define the barrier function B
    @variable model B Poly(monomials(vars, 0:deg_B))
    ∇B = differentiate(B, vars)
    Bdot = ∇B' * f
    

    # Constraints on the barrier function to ensure that the unsafe
    # set is never reached
    @constraint model -Bdot - p_space'*space >= 0
    # -B > ϵ on safe
    @constraint model -B-ϵ - p_safe'*safe >= 0
    # B >= 0 on unsafe
    @constraint model B - p_unsafe'*unsafe >= 0
    
    @constraint model p_space .>= 0
    @constraint model p_safe .>= 0
    @constraint model p_unsafe .>= 0


    # Solve the model
    optimize!(model)

    # Returns true if the model is successfully solved
    st = termination_status(model)   
    (st == MathOptInterface.ALMOST_OPTIMAL) | 
    (st == MathOptInterface.OPTIMAL) #|
    (st == MathOptInterface.SLOW_PROGRESS) 
end

   
##############################################
# Bisection for finding the optimal δ
##############################################

function find_δ_by_bisection(δmin=1e-2, δmax=10.; 
                             tol, r₀, r, ϵ, solver,
                             deg_B, deg_putinar)

    # Find the optimal δ by performing a bisection.

    best_δ = 1e6 # proxy for infinity
    
    while δmax - δmin > tol
        δ = (δmax+δmin)/2
        st_optimal = is_valid_upper_bound(δ=δ, r₀=r₀, r=r, ϵ=ϵ,
                deg_B=deg_B, deg_putinar=deg_putinar,  
                solver=solver)

        if st_optimal
            best_δ = δ
            δmax = δ
        else
            δmin = δ
        end

    end
    
    best_δ
end


##############################################
# Find the optimal δ corresponding to different 
# values of deg(B) and deg(putinar).
##############################################

if ~isinteractive()
  println("Computing the table of δs...")

  # deg(B) and deg putinar to try
  deg_Bs = 2:6
  deg_putinars = 2:2:6
  # 2D array to save the result
  δs = zeros(size(deg_Bs,1), size(deg_putinars,1))

  for (i, deg_B) in enumerate(deg_Bs)
      for (j, deg_putinar) in enumerate(deg_putinars)

          δ = find_δ_by_bisection(0., 1.; deg_B=deg_B,
                                  deg_putinar=deg_putinar, tol=tol,
                                  r₀=r₀, r=r, ϵ=ϵ, solver=solver)

          δs[i, j] = δ

          @show deg_B, deg_putinar, δ

      end
  end

  @show round.(δs, digits=3)

end
