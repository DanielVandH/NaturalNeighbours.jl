function eval_gradient(
    ::Iterative,
    tri::Triangulation,
    z,
    p,
    d_cache=DerivativeCache(tri),
    n_cache=NaturalNeighboursCache(tri);
    use_sibson_weight=false, kwargs...)
    return eval_gradient(Direct(), tri, z, p, d_cache, n_cache; use_sibson_weight, kwargs...) # Nothing to iterate in the order 1 case
end

function eval_gradient_and_hessian(
    ::Iterative,
    tri::Triangulation,
    z,
    p,
    d_cache=DerivativeCache(tri),
    n_cache=NaturalNeighboursCache(tri);
    use_sibson_weight=false, alpha=1/2, kwargs...)
    return _eval_gradient_and_hessian_iterative(tri, z, p, d_cache, n_cache; use_sibson_weight, alpha, kwargs...)
end

function _eval_gradient_and_hessian_iterative(
    tri,
    z,
    p,
    d_cache=DerivativeCache(tri),
    n_cache=NaturalNeighboursCache(tri);
    use_sibson_weight=false, alpha=1/2, kwargs...)
    g₁¹, g₂¹ = eval_gradient(Direct(), tri, z, p, d_cache, n_cache; use_sibson_weight, kwargs...)

end