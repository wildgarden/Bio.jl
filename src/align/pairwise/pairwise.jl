abstract PairwiseAlignmentAlgorithm

function score(a, b)
    return score!(NeedlemanWunsch(UnitScore), a, b)
end

function score!(alg::PairwiseAlignmentAlgorithm, a, b)
    return score!(alg, a, 1, length(a), b, 1, length(b))
end

function distance(a, b)
    return distance!(NeedlemanWunsch(UnitCost), a, b)
end

function distance!(alg::PairwiseAlignmentAlgorithm, a, b)
    return distance!(alg, a, 1, length(a), b, 1, length(b))
end

include("workspace.jl")
include("needleman_wunsch.jl")
include("short_detour.jl")
include("smith_waterman.jl")
include("semi_global.jl")
include("gotoh.jl")
