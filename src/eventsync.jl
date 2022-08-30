"""
    escor(X, tau, method)

Compute the event synchronization covariance matrix of the matrix `X`.

# Arguments
- `method::Function`: A method for thresholding potential spikes (e.g., `mean`, `median`)
- `X::AbstractMatrix{<:Union{Real, Missing}}`: A matrix of the form nvar by nobs (e.g., nodes by time).
- `tau::Int`: Window size (e.g., 6 hours)

# Returns
- `Matrix{Real}`: a square matrix representing the event synchronization correlation

# TODO
- Update the expected dimensionality of X to nobs by nvar to better fit with the rest of Simulation.jl
- Support a dim argument to support both X orientations (similar to how `cor` works)

# References:
- https://www.invenia.ca/wp-content/uploads/Correlations-and-Clustering-in-Wholesale-Electricity-Markets.pdf
- https://www.researchgate.net/publication/11025829_Event_Synchronization_A_simple_and_fast_method_to_measure_synchronicity_and_time_delay_patterns
"""
function escor(
    method::Union{typeof(median), typeof(mean)},
    X::AbstractMatrix{<:Union{Real, Missing}},
    tau::Int,
)
    # Remove columns that are entirely missing
    cols = .!missing_columns(X)
    X = @view X[:, cols]

    # Determine spike time series.
    Z = detectspikes(method, X)

    # Compute correlation inside a function barrier to handle type-instability associated
    # with `X`.
    return _function_barrier(Z, tau)
end

function _function_barrier(Z, tau)
    F = cumsum(Z; dims=2)
    H = Matrix{Float32}(undef, size(F))
    L = size(Z, 2)
    for col in 1:size(Z, 2), row in 1:size(Z, 1)
        if col - tau - 1 <= 0
            H[row, col] = sign(F[row, min(L, col + tau)])
        else
            H[row, col] = sign(F[row, min(L, col + tau)] - F[row, max(1, col - tau - 1)])
        end
    end
    J = H * H'
    D = pinv(Diagonal(sqrt.(Float64.(diag(J)))))
    return D * J * D
end

"""
    sqrtescor(method, X, tau)

Computes the square root of the event synchronization correlation matrix of `X`.
This has historically been useful if we want to use the event synchronization matrix with portfolio optimization.

See [`escor`](@ref) for more details.
"""
function sqrtescor(
    method::Union{typeof(median), typeof(mean)},
    X::AbstractMatrix{<:Union{Real, Missing}},
    tau::Int,
)
    M = escor(method, X, tau)
    M = (M + M') / 2
    return real(sqrt(M))
end

function detectspikes(
    rf::Union{typeof(median), typeof(mean)},
    X::AbstractMatrix{<:Union{Real, Missing}},
)
    Y = collect(X)'
    rfs_ub = map(col -> rf(Iterators.filter(x -> !ismissing(x) && x >= 0, col)), eachcol(Y))
    rfs_lb = map(col -> rf(Iterators.filter(x -> !ismissing(x) && x <= 0, col)), eachcol(Y))
    return coalesce.(X .> rfs_ub, false) .- coalesce.(X .< rfs_lb, false)
end

"""
    missing_columns(M) -> Vector{Bool}

Return `true` for columns which contain only `missing`. Used as a more performant
implementation of `vec(all(ismissing, M, 1))`.
"""
function missing_columns(M::AbstractMatrix)
    # Avoid using a BitVector with the current algorithm as we typically need to write to
    # every element.
    all_missing_cols = fill(true, size(M, 2))::Vector{Bool}
    rowinds, colinds = axes(M)
    for j in colinds
        for i in rowinds
            @inbounds if !ismissing(M[i,j])
                all_missing_cols[j] = false
                break
            end
        end
    end
    return all_missing_cols
end
