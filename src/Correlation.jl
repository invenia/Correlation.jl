module Correlation

using LinearAlgebra
using Missings
using Statistics: mean, median

if isdefined(LinearAlgebra, :sqrtm)
    sqrtm = LinearAlgebra.sqrtm
else
    sqrtm = LinearAlgebra.sqrt
end

include("eventsync.jl")

end # module
