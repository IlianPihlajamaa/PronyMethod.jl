module Prony

using PolynomialRoots, LinearAlgebra

export prony
export PronyMethod, PronyMethodLS, PronyMethodMPM

struct DampedCosineFit{T}
    p::Int64
    x_array::Array{T,1}
    y_array::Array{T,1}
    bases::Array{Complex{T}, 1}
    amplitudes::Array{Complex{T}, 1}
end

struct PronyMethod end
struct PronyMethodLS end
struct PronyMethodMPM end

function (fit::DampedCosineFit{T})(x::T) where T
    xmin = fit.x_array[1]
    xmax = fit.x_array[end]
    N = length(fit.x_array)
    indx = (x-xmin)/(xmax-xmin)*(N-1)
    result = zero(Complex{T})
    for i = 1:round(Int64,fit.p)
        result += fit.amplitudes[i] * fit.bases[i]^indx
    end
    return result
end

function construct_prediction_matrix(y, p)
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

function is_equidistant(x)
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
    p = length(z)
    vandermonde = zeros(eltype(z), length_vander, p)
    for row = 1:p
        for column = 1:length_vander
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
    if isodd(N)
        error("The number of data points must be even.")
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

# function least_squares(A::AbstractMatrix, b::AbstractVector)
#     @assert size(A, 1) == length(b)
#     N = size(A, 2)
#     C = [A b]
#     _, _, V = svd(C)
#     Vsub = V[1:N, 1+N:end]
#     Vlast = V[end, end]
#     if Vlast == 0 
#         error("Least squares not converged")
#     end
#     return -Vsub[:,1]/Vlast
# end

# function hankel(y, p)
#     N = length(y)
#     hank = zeros(eltype(y), N-p, p+1)
#     for row in axes(hank,2)
#         for col in axes(hank, 1)
#             hank[col, row] = y[(col-1)+(row-1)+1]
#         end
#     end
#     return hank
# end

# function find_amplitudes(bases, y, ::PronyMethodLS)
#     N = length(y)
#     vandermonde_matrix = transposed_vandermonde(bases, N)
#     return least_squares(vandermonde_matrix, y[1:size(vandermonde_matrix,1)])
# end

function find_amplitudes(bases, y, ::PronyMethod)
    p = length(bases)
    vandermonde_matrix = transposed_vandermonde(bases, p)
    return vandermonde_matrix\y[1:size(vandermonde_matrix,1)]
end

# function find_amplitudes(bases, y, ::PronyMethodMPM)
#     N = length(y)
#     vandermonde_matrix = transposed_vandermonde(bases, N)
#     return vandermonde_matrix\y
# end

# function find_bases(y, p, ::PronyMethodMPM)
#     Y = hankel(y, p)
#     Y1 = Y[:, 1:end-1]
#     Y2 = Y[:, 2:end]
#     bases = eigvals(pinv(Y1)*Y2)
#     return bases
# end

function find_bases(y, p, ::Any)
    coeffictients = find_coefficients(y, p)
    reverse!(coeffictients)
    push!(coeffictients, one(eltype(coeffictients)))
    polynomialroots = roots(coeffictients)
    return polynomialroots
end

# function prony(x::AbstractVector, y::AbstractVector, p::Integer, ::PronyMethodMPM)
#     check_args(x, y, p)
#     bases = find_bases(y, p, PronyMethodMPM())
#     amplitudes = find_amplitudes(bases, y, PronyMethodMPM())
#     return DampedCosineFit(p, x, y, bases, amplitudes)
# end

# function prony(x::AbstractVector, y::AbstractVector, p::Integer, ::PronyMethodLS)
#     check_args(x, y, p)
#     bases = find_bases(y, p, PronyMethodLS())
#     amplitudes = find_amplitudes(bases, y, PronyMethodLS())
#     return DampedCosineFit(p, x, y, bases, amplitudes)
# end

function prony(x::AbstractVector, y::AbstractVector, ::PronyMethod)
    p = round(Int, length(x)/2)
    check_args(x, y, p)
    bases = find_bases(y, p, PronyMethod())
    amplitudes = find_amplitudes(bases, y, PronyMethod())
    return DampedCosineFit(p, x, y, bases, amplitudes)
end

function prony(x::AbstractVector, y::AbstractVector, p::Integer)
    return prony(x, y, p, PronyMethodLS())
end

function prony(x::AbstractVector, y::AbstractVector)
    return prony(x, y, PronyMethod())
end

end # module
