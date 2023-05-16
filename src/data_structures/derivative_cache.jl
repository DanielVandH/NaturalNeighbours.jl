struct DerivativeCache{I, F}
    iterated_neighbourhood::Set{I}
    linear_matrix::ElasticMatrix{F, Vector{F}}
    quadratic_matrix::ElasticMatrix{F, Vector{F}}
    rhs_vector::Vector{F}
    linear_sol::Vector{F}
    quadratic_sol::Vector{F}
end
get_iterated_neighbourhood(cache::DerivativeCache) = cache.iterated_neighbourhood
get_linear_matrix(cache::DerivativeCache) = cache.linear_matrix
get_quadratic_matrix(cache::DerivativeCache) = cache.quadratic_matrix
get_rhs_vector(cache::DerivativeCache) = cache.rhs_vector
get_linear_sol(cache::DerivativeCache) = cache.linear_sol
get_quadratic_sol(cache::DerivativeCache) = cache.quadratic_sol
function DerivativeCache(tri::Triangulation{P,Ts,I,E,Es,BN,BNM,B,BIR,BPL}) where {P,Ts,I,E,Es,BN,BNM,B,BIR,BPL}
    iterated_neighbourhood = Set{I}()
    F = number_type(tri)
    linear_matrix = ElasticMatrix{F, Vector{F}}(undef, 2, 0)
    quadratic_matrix = ElasticMatrix{F, Vector{F}}(undef, 9, 0)
    rhs_vector = zeros(F, 0)
    linear_sol = zeros(F, 2)
    quadratic_sol = zeros(F, 9)
    sizehint_m = 2^4
    sizehint!(iterated_neighbourhood, sizehint_m)
    sizehint!(linear_matrix, (2, sizehint_m))
    sizehint!(quadratic_matrix, (9, sizehint_m))
    sizehint!(rhs_vector, 2^4)
    return DerivativeCache(iterated_neighbourhood, linear_matrix, quadratic_matrix, rhs_vector, linear_sol, quadratic_sol)
end