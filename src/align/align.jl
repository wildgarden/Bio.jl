module Align

export
    align,
    GAP,
    distance, distance!,
    UnitScore,
    BLOSUM62,
    UnitCost,
    CostModel,
    score, score!,
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
