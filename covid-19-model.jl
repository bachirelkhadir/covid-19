##############################################
# Import libraries
##############################################

using SumOfSquares
using DynamicPolynomials 
using DifferentialEquations
using Plots
using JuMP
using MathOptInterface
using CSDP
using MosekTools


##############################################
# Parameters
##############################################

β, μ, α = 1.0, 0.8, .1
D, Z  = 4, 4 # Days
N = 1.0
default_parameters=[β, μ, α, D, Z]

r₀=1.
r=.1
ϵ = .1
tol = 1e-3
#solver = optimizer_with_attributes(CSDP.Optimizer, "printlevel" => 0)
solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
@show solver

##############################################
# Vector Field
##############################################

function transmission_vf(S, E, Ir, Iu; 
    paramters=default_parameters)
    β, μ, α, D, Z = paramters
    [
        β * S/N * Ir - μ * β * S/N * Iu,
        β * S/N * Ir + μ * β * S/N * Iu - E/Z,
        α * E/Z - Ir/D,
        (1-α) * E/Z - Iu/D
    ]
end

