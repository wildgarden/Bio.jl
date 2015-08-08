module Align

export
    distance, distance!,
    align,
    AlignmentMatrix,
    AlignmentVector,
    GAP,
    UnitCost,
    CostModel,
    # algorithms
    NaiveDP,
    ShortDetourDP

import Base:
    getindex,
    setindex!,
    size,
    resize!,
    empty!,
    reverse!,
    push!,
    length,
    endof,
    show

using Bio.Seq

include("model.jl")
include("pairwise.jl")

end
