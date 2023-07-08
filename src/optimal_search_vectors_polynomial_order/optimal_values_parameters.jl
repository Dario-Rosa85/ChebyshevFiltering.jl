function delta_prime_estimate(approximated_density, search_vector_numbers, delta_prime_list, lambda_min, lambda_max, e_min, e_max)
    @inbounds for i in delta_prime_list 
        x_list_integral = range(max(lambda_min - i, e_min + 10^-5), min(lambda_max + i, e_max - 10^-5), length=500)
        y_list_integral = map(x -> approximated_density(x), x_list_integral)
        expected_vector = NumericalIntegration.integrate(x_list_integral, y_list_integral)
        if expected_vector > search_vector_numbers
            return i
        end
    end
    return "Error"
end

function optimal_search_degree(hamiltonian_matrix, lambda_min, lambda_max, search_target_ratio, max_degree_KPM, stochastic_dimension_KPM, e_max, e_min, N_0)
    approximated_density = KPM_density(hamiltonian_matrix, max_degree_KPM, stochastic_dimension_KPM)
    x_list_integral = range(lambda_min, lambda_max, length=500)
    y_list_integral = map(x -> approximated_density(x), x_list_integral)
    search_vector_numbers = search_target_ratio * Int64(floor(NumericalIntegration.integrate(x_list_integral, y_list_integral)))
    delta_prime_list = range(0,  search_target_ratio * (lambda_max - lambda_min) / 2, 50)
    delta_prime = delta_prime_estimate(approximated_density, search_vector_numbers, delta_prime_list, lambda_min, lambda_max, e_min, e_max)
    polynomial_degree_optim = Int64(ceil(N_0 * (e_max - e_min) / delta_prime))
    return search_vector_numbers, polynomial_degree_optim
end