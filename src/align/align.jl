module Align

export
    align,
    GAP,
    distance, distance!,
    UnitScore,
    ScoreModel,
    AffineScoreModel,
    BLOSUM62,
    UnitCost,
    CostModel,
    score, score!,
    # algorithms
    NeedlemanWunsch,
    ShortDetour,
    SmithWaterman,
    Gotoh,
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

include("substition_matrix.jl")
include("model.jl")
include("pairwise/pairwise.jl")

end
