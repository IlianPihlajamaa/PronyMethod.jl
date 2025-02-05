module Prony

using PolynomialRoots, LinearAlgebra

export prony
export PronyMethod, ApproximatePronyMethod

# Fitting object returned by prony: 
struct DampedCosineFit{T} 
    p::Int64
    x_array::Array{T,1}
    y_array::Array{Complex{T},1}

    # y = A*B^{ (x-x0)(N-1)/(xmax-xmin) }
    bases::Array{Complex{T}, 1}
    amplitudes::Array{Complex{T}, 1}
    
    # y = A*exp(-λ(t-t0))
    exponents::Array{Complex{T}, 1}

    function DampedCosineFit(p::Int64, x_array::Array{T,1}, y_array::Array{Complex{T},1}, bases::Array{Complex{T}, 1}, amplitudes::Array{Complex{T}, 1}) where T
        exponents::Array{Complex{T},1} = log.(bases) * (length(x_array)-1)/(x_array[end]-x_array[1])
        new{T}(p, x_array, y_array, bases, amplitudes, exponents)
    end
end

struct PronyMethod end
struct ApproximatePronyMethod 
    decay::Bool
end

ApproximatePronyMethod() = ApproximatePronyMethod(false)

function (fit::DampedCosineFit{T})(x::T) where T
    xmin = fit.x_array[1]
    xmax = fit.x_array[end]
    N = length(fit.x_array)
    indx = (x-xmin)/(xmax-xmin)*(N-1)
    result = zero(T)
    for i = 1:round(Int64,fit.p)
        result += (fit.amplitudes[i] * fit.bases[i]^indx)
    end
    return result
end

function construct_prediction_matrix(y, p)  # H_ij = y_{i+j+1}
    N = length(y)
    prediction_matrix = zeros(eltype(y), N-p, p)
    for row = 1:p
        for column = 1:N-p
            prediction_matrix[column, row] = y[p + (column-1) - (row-1)]
        end
    end
    return prediction_matrix
end

function construct_observation_vector(y, p)
    N = length(y)
    observation_vector = zeros(eltype(y), N-p)
    for (i,ii) in  enumerate(p+1:N)
        observation_vector[i] = -y[ii]
    end
    return observation_vector
end

function is_equidistant(x) # x input must be equidistant or code returns error
    N = length(x)
    dx = x[2] - x[1]
    for i = 2:N-1
        if !(x[i+1] - x[i] ≈ dx)
            return false
        end
    end
    return true
end

function transposed_vandermonde(z, length_vander)
    N = length(z)
    vandermonde = zeros(eltype(z), length_vander, N)
    for row in axes(vandermonde, 2)
        for column in axes(vandermonde, 1)
            vandermonde[column, row] = z[row]^(column-1)
        end
    end
    return vandermonde
end

function find_coefficients(y, p)
    N = length(y)
    prediction_matrix = construct_prediction_matrix(y, p)
    observation_vector = construct_observation_vector(y, p)
    if N == 2p
        return prediction_matrix\observation_vector
    elseif N > 2p
        return least_squares(prediction_matrix, observation_vector)
    end
end

function check_args(x, y, p)
    N = length(x)
    if length(y) != N
        error("length of x must equal that of y.")
    end
    if 2p > N
        error("2p must be smaller than the number of data points.")
    end
    if !issorted(x)
        error("x must be a sorted array.")
    end
    if !is_equidistant(x)
        error("x must be an equidistant grid.")
    end
end

function find_amplitudes(bases, y, ::PronyMethod)
    p = length(bases)
    vandermonde_matrix = transposed_vandermonde(bases, p)
    return vandermonde_matrix\y[1:size(vandermonde_matrix,1)]
end

function find_bases(y, p, ::PronyMethod)
    coeffictients = find_coefficients(y, p)
    reverse!(coeffictients)
    push!(coeffictients, one(eltype(coeffictients)))
    polynomialroots = roots(coeffictients)
    return polynomialroots
end

# Prony method with length(x) exponentials
function prony(x::AbstractVector, y::AbstractVector, ::PronyMethod)
    if isodd(length(x))
        pop!(x)
        pop!(y)
    end
    p = round(Int, length(x)/2)
    check_args(x, y, p)
    bases = find_bases(y, p, PronyMethod())
    amplitudes = find_amplitudes(bases, y, PronyMethod())
    return DampedCosineFit(p, x, Complex.(y), bases, amplitudes)
end

# Prony method with M exponentials
function prony(x::AbstractVector, y::AbstractVector, M::Int64, m::ApproximatePronyMethod)
    if isodd(length(x))
        pop!(x)
        pop!(y)
    end
    N = round(Int, (length(x))/2)
    if M > N
        error("Number of requested functions is too large. Must be smaller than N, where 2N+1 is the number of data-points")
    end

    bases = find_bases(y, N, PronyMethod())
    amplitudes = find_amplitudes(bases, y, PronyMethod())

    bases[isnan.(amplitudes)] .= 0.0
    amplitudes[isnan.(amplitudes)] .= 0.0
    if m.decay == true
        amplitudes[abs.(bases) .> 1] .= 0.0
    end

    sortedindx = sortperm(abs.(amplitudes), rev=true)[1:M]
    bases = bases[sortedindx]
    amplitudes = amplitudes[sortedindx]

    return DampedCosineFit(length(bases), x, Complex.(y), bases, amplitudes)
end


function prony(x::AbstractVector, y::AbstractVector, N::Int64; decay=false)
    return prony(x, y, N, ApproximatePronyMethod(decay))
end

function prony(x::AbstractRange, y::AbstractVector, N::Int64; decay=false)
    return prony(collect(x), y, N, ApproximatePronyMethod(decay))
end

function prony(x::AbstractVector, y::AbstractVector)
    return prony(x, y, PronyMethod())
end

function prony(x::AbstractRange, y::AbstractVector)
    return prony(collect(x), y, PronyMethod())
end

end # module
