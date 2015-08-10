module Align

export
    align,
    GAP,
    distance, distance!,
    UnitScore,
    ScoreModel,
    AffineScoreModel,
    UnitCost,
    CostModel,
    # substitution matrices
    BLOSUM45,
    BLOSUM50,
    BLOSUM62,
    BLOSUM80,
    BLOSUM90,
    PAM30,
    PAM70,
    PAM250,
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
