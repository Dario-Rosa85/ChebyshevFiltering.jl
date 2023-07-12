# ChebyshevFiltering.jl

A Julia code for high-performance Chebyshev polynomial diagonalization. In the current version, the package is built with the purpose to find all the eigenvalues and the corresponding eigenvectors -- of a given sparse and Hermitian (both real and complex) matrix -- lying in a chosen interval $\left[\lambda_{\mathrm{min}}, \, \lambda_{\mathrm{max}} \right]$.

It uses the mathematical algorithm developed in [Journal of Computational Physics 325, 226 (2016)](https://www.sciencedirect.com/science/article/abs/pii/S0021999116303837?via%3Dihub). As such, it renormalizes authomatically the matrix such that its spectrum lies in the interval $\left[-1, \, 1 \right]$.

## Example of usage

In the current version, the package provides the function:

    eigen_cheb(hamiltonian_matrix, lambda_min::Float64=-0.1, lambda_max::Float64=0.1, log_path="none", log_file_name="none", search_target_ratio::Float64=3.0, max_degree_KPM::Int64=150, stochastic_dimension_KPM::Int64=30, N_0::Float64=6.23, epsilon_convergence::Float64=10^(-7))

that requires as input the matrix ```hamiltonian_matrix``` of interest and gives back a tuple containing a vector with all the eigenvalues found in the interval between ```lambda_min``` and ```lambda_max``` and another vector containing all the corresponding eigenvectors. Optionally, the user can provide the path to the folder where the log file is located and the name of the log file itself. The extra variables ```search_target_ratios```, ```max_degree_KPM```, ```stochastic_dimension_KPM```, ``N_0``` and ```epsilon_convergence``` are all defined in [Journal of Computational Physics 325, 226 (2016)](https://www.sciencedirect.com/science/article/abs/pii/S0021999116303837?via%3Dihub). The default values should work optimally for most eigenvalue problems of interest in many-body physics.
