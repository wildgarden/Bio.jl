module Align

export
    align,
    AlignmentMatrix,
    AlignmentVector,
    GAP,
    distance, distance!,
    UnitCost,
    CostModel,
    score, score!,
    score_local, score_local!,
    UnitScore,
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
