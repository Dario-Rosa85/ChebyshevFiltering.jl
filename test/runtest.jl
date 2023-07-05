using ChebyshevFiltering
using Test
using SparseArrays
using LinearAlgebra

hamiltonian = sparse(rand(ComplexF64, (100, 100)))
hamiltonian .= hamiltonian + hamiltonian' 

@testset "renormalization_hamiltonian!" begin
    renormalization_hamiltonian!(hamiltonian)
    spectrum = eigvals(Matrix(hamiltonian))
    e_min = spectrum[1]
    e_max = spectrum[end]
    @test isapprox(e_min, -1.0)
    @test isapprox(e_max, 1.0)
end