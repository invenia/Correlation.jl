using Base: @deprecate

@deprecate escor(
    X::AbstractMatrix{<:Union{Real, Missing}},
    tau::Int,
    method::Union{typeof(median), typeof(mean)}
) escor(method, X, tau)

@deprecate sqrtescor(
    X::AbstractMatrix{<:Union{Real, Missing}},
    tau::Int,
    method::Union{typeof(median), typeof(mean)}
) sqrtescor(method, X, tau)

@deprecate detectspikes(
    X::AbstractMatrix{<:Union{Real, Missing}},
    rf::Union{typeof(median), typeof(mean)}
) detectspikes(rf, X)

@deprecate threshold_p(
    X::AbstractMatrix,
    rf::Union{typeof(median), typeof(mean)}
) threshold_p(rf, X)

@deprecate threshold_n(
    X::AbstractMatrix,
    rf::Union{typeof(median), typeof(mean)}
) threshold_n(rf, X)
