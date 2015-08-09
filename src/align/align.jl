module Align

export
    align,
    GAP,
    distance, distance!,
    UnitCost,
    CostModel,
    score, score!,
    UnitScore,
    # algorithms
    NeedlemanWunsch,
    ShortDetour,
    SmithWaterman,
    SemiGlobal

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
include("pairwise/pairwise.jl")

end
