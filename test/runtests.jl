using ChebyshevFiltering
using Test
using SparseArrays
using LinearAlgebra

hamiltonian_matrix = sparse(randn(Float64, (300, 300)))
hamiltonian_matrix = hamiltonian_matrix + hamiltonian_matrix'

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
    max_degree_KPM = 150
    stochastic_dimension_KPM = 30
    real_number = length([value  for value in eigvals(Matrix(hamiltonian_matrix)) if lambda_min < value < lambda_max])
    approximated_density = KPM_density(hamiltonian_matrix, max_degree_KPM, stochastic_dimension_KPM)
    x_values = lambda_min:((lambda_max - lambda_min) / 200):lambda_max
    y_values = [approximated_density(x) for x in x_values]

    using NumericalIntegration
    expected_value = NumericalIntegration.integrate(x_values, y_values)

    @test (real_number - 2 <= expected_value <= real_number + 2)
    
end

@testset "optimal_parameters_estimate" begin
    ChebyshevFiltering.renormalization_hamiltonian!(hamiltonian_matrix)
    max_degree_KPM = 150
    stochastic_dimension_KPM = 30
    lambda_min = -0.6
    lambda_max = -0.4
    search_target_ratio = 3 
    e_max = 1.
    e_min = -1. 
    N_0 = 6.23
    real_number = length([value  for value in eigvals(Matrix(hamiltonian_matrix)) if lambda_min < value < lambda_max])
    search_vector_numbers, polynomial_degree_optim = ChebyshevFiltering.optimal_search_degree(hamiltonian_matrix, lambda_min, lambda_max, search_target_ratio, max_degree_KPM, stochastic_dimension_KPM, e_max, e_min, N_0)

    @test abs(search_vector_numbers - search_target_ratio * real_number) <= 10
end