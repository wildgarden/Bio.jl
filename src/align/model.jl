abstract AlignmentModel

# cost model is about a minimizing problem
abstract AbstractCostModel <: AlignmentModel

immutable GAP; end
typealias Character Union{Char,Nucleotide,AminoAcid}

# a gap and a mismatch costs 1; otherwise 0
immutable UnitCostModel <: AbstractCostModel; end
const UnitCost = UnitCostModel()

getindex(::UnitCostModel,  ::Character,  ::Type{GAP}) = 1
getindex(::UnitCostModel,  ::Type{GAP},  ::Character) = 1
getindex(::UnitCostModel, x::Character, y::Character) = ifelse(x === y, 0, 1)

minimum_indel_cost(::UnitCostModel) = 1

type CostModel{T<:Character} <: AbstractCostModel
    cost::Matrix{Int}
    char2gap::Int
    gap2char::Int
end

getindex{T}(m::CostModel{T}, ::T, ::Type{GAP}) = m.char2gap
getindex{T}(m::CostModel{T}, ::Type{GAP}, ::T) = m.gap2char
getindex{T}(m::CostModel{T}, x::T, y::T) = m.cost[x, y]


# score model is a maximizing problem
abstract AbstractScoreModel <: AlignmentModel

# PAM, BLOSUM, etc.
