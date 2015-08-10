abstract PairwiseAlignmentAlgorithm

function score!{A<:PairwiseAlignmentAlgorithm}(alg::A, a, b)
    return score!(alg, a, 1, length(a), b, 1, length(b))
end

function distance!{A<:PairwiseAlignmentAlgorithm}(alg::A, a, b)
    return distance!(alg, a, 1, length(a), b, 1, length(b))
end

include("workspace.jl")
include("needleman_wunsch.jl")
include("short_detour.jl")
include("smith_waterman.jl")
include("semi_global.jl")
