module ChebyshevFiltering

using KrylovKit
using SparseArrays
using Polynomials
using LinearAlgebra
using NumericalIntegration
using Dates

export eigen_cheb

include(joinpath(@__DIR__, "hamiltonian_rescaling/renormalization_hamiltonian.jl"))
include(joinpath(@__DIR__, "KPM_density/KPM_density.jl"))
include(joinpath(@__DIR__, "optimal_search_vectors_polynomial_order/optimal_values_parameters.jl"))
include(joinpath(@__DIR__, "orthogonalization_routines/orthogonalization.jl"))
include(joinpath(@__DIR__, "expansion_damping_coeffs/expansion_damping_coeffs.jl"))
include(joinpath(@__DIR__, "polynomial_filtering/polynomial_filtering.jl"))
include(joinpath(@__DIR__, "wrappers/eigen_cheb.jl"))

end # module ChebyshevFiltering
