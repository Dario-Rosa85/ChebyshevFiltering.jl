using ChebyshevFiltering
using Test
using SparseArrays
using LinearAlgebra
include("/home/dario/Dropbox/delocalization_and_growth/revision_chaos_paper_code/QMB_module.jl")
using .ImpurityUtilities

const number_of_majorana = 20
majorana_matrices = majorana_matrices_chiral(number_of_majorana)
hamiltonian_matrix, edges_graph, weights_graph = hamiltonian_chiral(number_of_majorana, majorana_matrices, 4, 1.0)

@testset "renormalization_hamiltonian!" begin
    ChebyshevFiltering.renormalization_hamiltonian!(hamiltonian_matrix)
    spectrum = eigvals(Matrix(hamiltonian_matrix))
    e_min = spectrum[1]
    e_max = spectrum[end]
    @test isapprox(e_min, -1.0)
    @test isapprox(e_max, 1.0)
end

@testset "KPM_density estimation" begin
    ChebyshevFiltering.renormalization_hamiltonian!(hamiltonian_matrix)
    lambda_min = -0.1
    lambda_max = 0.1
    real_number = length([value  for value in eigvals(Matrix(hamiltonian_matrix)) if lambda_min < value < lambda_max])
    approximated_density = KPM_density(hamiltonian_matrix)
    x_values = lambda_min:((lambda_max - lambda_min) / 200):lambda_max
    y_values = [approximated_density(x) for x in x_values]

    using NumericalIntegration
    expected_value = NumericalIntegration.integrate(x_values, y_values)

    @test (real_number - 5 <= expected_value <= real_number + 5)
    
end