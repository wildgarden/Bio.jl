module Align

export
    distance, distance!,
    AlignmentMatrix,
    UnitCost,
    # algorithms
    NaiveDP,
    ShortDetourDP

import Base: getindex, setindex!, size, resize!, empty!

using Bio.Seq

include("model.jl")
include("algorithm.jl")
include("distance.jl")

end
