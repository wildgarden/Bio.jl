module Align

export
    distance, distance!, distance2, distance2!,
    AlignmentMatrix,
    UnitCost

import Base: getindex, setindex!, size, resize!, empty!

using Bio.Seq

include("model.jl")
include("algorithm.jl")
include("distance.jl")

end
