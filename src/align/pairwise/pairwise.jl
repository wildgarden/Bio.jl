abstract PairwiseAlignmentAlgorithm

function score{A<:PairwiseAlignmentAlgorithm}(a, b, score::AbstractScoreModel=UnitScore, ::Type{A}=NeedlemanWunsch)
    return score!(A(), a, b, score)
end

function score!{A<:PairwiseAlignmentAlgorithm}(alg::A, a, b, score::AbstractScoreModel=UnitScore, ::Type{A}=NeedlemanWunsch)
    return score!(alg, a, 1, length(a), b, 1, length(b), score)
end

function distance{A<:PairwiseAlignmentAlgorithm}(a, b, cost::AbstractCostModel=UnitCost, ::Type{A}=NeedlemanWunsch)
    return distance!(A(), a, b, cost)
end

function distance!{A<:PairwiseAlignmentAlgorithm}(alg::A, a, b, cost::AbstractCostModel)
    return distance!(alg, a, 1, length(a), b, 1, length(b), cost)
end

include("workspace.jl")
include("needleman_wunsch.jl")
include("short_detour.jl")
include("smith_waterman.jl")
