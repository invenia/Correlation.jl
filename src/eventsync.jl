"""
    kron!(X, v)

Compute ``X + (v \\otimes v')``, i.e. the sum of the matrix `X` and the Kronecker product
of `v` with its transpose, and store the result in `X`.

!!! warn
    For performance reasons, no error condition checking is performed in this function; it
    assumes that `X` is square with dimensions `length(v) × length(v)`.
"""
function kron!(X::AbstractMatrix, v::AbstractVector)
    # NOTE: The lack of error checking is okay for our use case, since `X` is square by
    # construction with dimensions conformant with `v`
    @inbounds for i = eachindex(v), j = eachindex(v)
        X[i,j] += v[i] * v[j]
    end
    X
end

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

    N = size(X, 1)
    L = size(X, 2)
    Z = detectspikes(method, X)
    J = zeros(N, N)

    for k=1:L
        Zt = view(Z, :, max(1, k-tau):min(L,k+tau))
        t = vec(sign.(sum(Zt, dims=2)))
        kron!(J, t)
    end

    D = pinv(Diagonal(sqrt.(diag(J))))
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
    return real(sqrtm(M))
end

function detectspikes(
    rf::Union{typeof(median), typeof(mean)},
    X::AbstractMatrix{<:Union{Real, Missing}},
)
    Z1 = threshold_p(rf, X)
    Z2 = threshold_n(rf, X)

    return Z1 .- Z2
end

# This is easier to read without Missings.
# threshold_p(rf::Function, X::AbstractMatrix) = mapslices(x -> x.>rf(x[x.>=0]), P, 2)
# The additional mess is to replicate MATLAB behaviour: [1,2, NaN] > 0 == [true, true, false]
# This is the Missings.replace(x.>=0, false) part.
# The collect statements are because Missings is lazy
# The final Missings.replace(x, -Inf)).>rf(...) is because x.>rf(...) appears to fail
function threshold_p(rf::Union{typeof(median), typeof(mean)}, X::AbstractMatrix)
    rfs = mapslices(row -> rf([x for x in row if !ismissing(x) && x >= 0]), X, dims=2)
    collect(Missings.replace(X .> rfs, false))
end

function threshold_n(rf::Union{typeof(median), typeof(mean)}, X::AbstractMatrix)
    rfs = mapslices(row -> rf([x for x in row if !ismissing(x) && x <= 0]), X, dims=2)
    collect(Missings.replace(X .< rfs, false))
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
