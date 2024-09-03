# ChebyshevFiltering.jl

A Julia code for high-performance Chebyshev polynomial diagonalization. In the current version, the package is built with the purpose to find all the eigenvalues and the corresponding eigenvectors, lying in a chosen interval $\left[\lambda_{\mathrm{min}}, \ \lambda_{\mathrm{max}} \right]$, - for a given *sparse* and *Hermitian* (real or complex) matrix.

The code makes use of multithreading to speed up the computation and so it works at its best on single-node, multicore, configurations. From its very definition, the algorithm works in such a way to keep the sparsity of the problem at every step, and it works at its best when the interval chosen, $\left[\lambda_{\mathrm{min}}, \ \lambda_{\mathrm{max}} \right]$, contains just few hundreds of eigenvalues. In particular, it *does not* work well in finding large chunks of eigenvalues and eigenvectors.

Mathematically, it uses the algorithm developed in [Journal of Computational Physics 325, 226 (2016)](https://www.sciencedirect.com/science/article/abs/pii/S0021999116303837?via%3Dihub). As such, it renormalizes authomatically the matrix such that its spectrum lies in the interval $\left[-1, \ 1 \right]$. The current Julia implementation has been used in ["From Dyson Models to Many-Body Quantum Chaos", arXiv:2302.00917 [quant-ph]](https://arxiv.org/abs/2302.00917) to discuss the quantum chaotic properties of a family of impurity SYK models.

## How to install it

```ChebyshevFiltering.jl``` can be installed with the Julia package manager. From the REPL, just type ```]``` to enter in PKG mode and run:

    pkg> add ChebyshevFiltering

The current version has been tested on Julia v1.9.2.

## Example of usage

In the current version, the package provides the function:

    eigen_cheb(hamiltonian_matrix; lambda_min::Float64=-0.1, lambda_max::Float64=0.1, n_chunks_thread::Int64=2, n_first_iterations::Int64=3, epsilon_convergence::Float64=10^(-7), log_path="none", log_file_name="none", search_target_ratio::Float64=3.0, max_degree_KPM::Int64=150, stochastic_dimension_KPM::Int64=30, N_0::Float64=6.23)

that requires as input the sparse and hermitian matrix ```hamiltonian_matrix``` of interest and gives back a tuple made by a vector with all the eigenvalues found in the interval between ```lambda_min``` and ```lambda_max``` and by another vector containing all the corresponding eigenvectors. Optionally, the user can provide the path to the folder where the log file is located and the name of the log file itself. The extra variables ```search_target_ratios```, ```max_degree_KPM```, ```stochastic_dimension_KPM```, ```N_0``` and ```epsilon_convergence``` are all defined in [Journal of Computational Physics 325, 226 (2016)](https://www.sciencedirect.com/science/article/abs/pii/S0021999116303837?via%3Dihub). The default values should work optimally for most eigenvalue problems of interest in many-body physics.
The variable ```n_chunks_thread``` set the number of chunks that a single thread should work on during parallelization, and the the variable ```n_first_iteration``` set the number of filtering steps the algorithm performs before checking for convergence of the Ritz vectors.

## Questions and Contributions

Contributions, suggestions and, especially, feature requests are very much appreciated. Please open an issue for any feedback.

## How to cite

In case you find this code useful for your research projects, please add a citation to the following papers in your reference list:

* [A. Andreanov, M. Carrega, J. Murugan, J. Olle, D. Rosa, R. Shir, "*From Dyson Models to Many-Body Quantum Chaos*", arXiv:2302.00917 [quant-ph].](https://arxiv.org/abs/2302.00917)

* [A. Pieper et al., "*High-performance implementation of Chebyshev filter diagonalization for interior eigenvalue computations*", Journal of Computational Physics 325, 226 (2016).](https://www.sciencedirect.com/science/article/abs/pii/S0021999116303837?via%3Dihub)
