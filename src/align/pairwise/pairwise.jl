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
     return distance!(alg, a, 1, length(a), b, 1, length(b), false)
end


"returns an optimal alignment. Other optimal alignments may exist. See funtion alignAll()"
function align(a, b)
    return align(NeedlemanWunsch(UnitCost), a, b)
end


"returns an optimal alignment. Other optimal alignments may exist. See funtion alignAll()"
function align(alg::PairwiseAlignmentAlgorithm, a, b)
    distance!(alg, a, 1, length(a), b, 1, length(b), false)

    return backtrackSingle(alg, a, 1, b, 1)
end


"Returns all optimal alignments of two sequences using the UniCost model"
function alignAll(a, b)
    return alignAll(NeedlemanWunsch(UnitCost), a, b)
end


"returns all optimal alignments"
function alignAll(alg::PairwiseAlignmentAlgorithm, a, b)
    distance!(alg, a, 1, length(a), b, 1, length(b), false)
    (row, col) = size(alg.trace_matrix)
    result = []
    backtrackAll!(alg, result, row, col, a, b)

    return result
end


include("workspace.jl")
include("needleman_wunsch.jl")
include("short_detour.jl")
include("smith_waterman.jl")
include("semi_global.jl")
include("gotoh.jl")
