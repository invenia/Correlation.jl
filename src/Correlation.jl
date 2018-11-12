module Correlation

using Compat: Compat, axes, sum
using Compat.LinearAlgebra
using Compat.Statistics: mean, median
using Missings

if isdefined(Compat.LinearAlgebra, :sqrtm)
    sqrtm = Compat.LinearAlgebra.sqrtm
else
    sqrtm = Compat.LinearAlgebra.sqrt
end

include("eventsync.jl")
include("deprecated.jl")

end # module
