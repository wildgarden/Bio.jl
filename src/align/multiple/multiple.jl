abstract MultipleAlignmentAlgorithm

#--------------------------------------------------------------------------------------------------

include("center_star.jl")


"returns all optimal alignments"
function align(sequences::Array)
    return align(CenterStar(UnitCost), sequences)
end
