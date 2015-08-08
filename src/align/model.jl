abstract AlignmentModel

# cost model is about a minimizing problem
abstract AbstractCostModel <: AlignmentModel

# gap character: -
immutable GAP; end
typealias Character Union{Char,Nucleotide,AminoAcid}

# a gap and a mismatch costs 1; otherwise 0
immutable UnitCostModel <: AbstractCostModel; end
const UnitCost = UnitCostModel()

getindex(::UnitCostModel,  ::Character,  ::Type{GAP}) = 1
getindex(::UnitCostModel,  ::Type{GAP},  ::Character) = 1
getindex(::UnitCostModel, x::Character, y::Character) = ifelse(x === y, 0, 1)

minimum_indel_cost(::UnitCostModel) = 1

# arbitrary cost model
type CostModel{T<:Character} <: AbstractCostModel
    cost::Matrix{Int}
    char2gap::Int
    gap2char::Int
    function CostModel(alphabetsize)
        new(zeros(Int, alphabetsize, alphabetsize), 0, 0)
    end
end

getindex{T}(m::CostModel{T},  ::T,          ::Type{GAP}) = m.char2gap
getindex{T}(m::CostModel{T},  ::Type{GAP},  ::T        ) = m.gap2char
getindex{T}(m::CostModel{T}, x::T,         y::T        ) = m.cost[convert(UInt8,x)+1,convert(UInt8,y)+1]

macro check_cost(ex)
    quote
        local cost = $(esc(ex))
        cost ≥ 0 || error("cost should be non-negative (given ", cost, ")")
        cost
    end
end

setindex!{T}(m::CostModel{T}, c::Int,  ::Type{GAP}              ) = m.char2gap = m.gap2char = @check_cost c
setindex!{T}(m::CostModel{T}, c::Int,  ::T,          ::Type{GAP}) = m.char2gap = @check_cost c
setindex!{T}(m::CostModel{T}, c::Int,  ::Type{GAP},  ::T        ) = m.gap2char = @check_cost c
setindex!{T}(m::CostModel{T}, c::Int, x::T,         y::T        ) = m.cost[convert(UInt8,x)+1,convert(UInt8,y)+1] = @check_cost c


# score model is a maximizing problem
abstract AbstractScoreModel <: AlignmentModel

immutable UnitScoreModel <: AbstractScoreModel; end
const UnitScore = UnitScoreModel()

getindex(::UnitScoreModel,  ::Character,  ::Type{GAP}) = -1
getindex(::UnitScoreModel,  ::Type{GAP},  ::Character) = -1
getindex(::UnitScoreModel, x::Character, y::Character) = ifelse(x === y, 1, 0)

type ScoreModel{T<:Character} <: AbstractScoreModel
    score::Matrix{Int}
    char2gap::Int
    gap2char::Int
    function ScoreModel(alphabetsize)
        new(zeros(Int, alphabetsize, alphabetsize), 0, 0)
    end
end

getindex{T}(m::ScoreModel{T},  ::T,           ::Type{GAP}) = m.char2gap
getindex{T}(m::ScoreModel{T},  ::Type{GAP},  ::T         ) = m.gap2char
getindex{T}(m::ScoreModel{T}, x::T,         y::T         ) = m.score[convert(UInt8,x)+1,convert(UInt8,y)+1]

setindex!{T}(m::ScoreModel{T}, c::Int,  ::Type{GAP}              ) = m.char2gap = m.gap2char = c
setindex!{T}(m::ScoreModel{T}, c::Int,  ::T,          ::Type{GAP}) = m.char2gap = c
setindex!{T}(m::ScoreModel{T}, c::Int,  ::Type{GAP},  ::T        ) = m.gap2char = c
setindex!{T}(m::ScoreModel{T}, c::Int, x::T,         y::T        ) = m.score[convert(UInt8,x)+1,convert(UInt8,y)+1] = c

# PAM, BLOSUM, etc.
