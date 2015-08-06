module Align

export distance, distance2

import Base: getindex, setindex!, size

using Bio.Seq

include("model.jl")
include("algorithm.jl")
include("distance.jl")

end
