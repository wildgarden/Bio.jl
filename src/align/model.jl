abstract AlignmentModel

# gap character: -
immutable GAP; end
typealias Character Union{Char,Nucleotide,AminoAcid}

# score model is a maximizing problem
abstract AbstractScoreModel{T<:Real} <: AlignmentModel

immutable UnitScoreModel{T} <: AbstractScoreModel{T}; end
const UnitScore = UnitScoreModel{Int}()

getindex{T}(::UnitScoreModel{T},  ::Character,  ::Type{GAP}) = T(-1)
getindex{T}(::UnitScoreModel{T},  ::Type{GAP},  ::Character) = T(-1)
getindex{T}(::UnitScoreModel{T}, x::Character, y::Character) = ifelse(x === y, T(1), T(-1))

type ScoreModel{C<:Character,T} <: AbstractScoreModel
    score::Matrix{T}
    char2gap::T
    gap2char::T
    function ScoreModel(alphabetsize)
        new(zeros(T, alphabetsize, alphabetsize), 0, 0)
    end
end

getindex{C}(m::ScoreModel{C},  ::C,           ::Type{GAP}) = m.char2gap
getindex{C}(m::ScoreModel{C},  ::Type{GAP},  ::C         ) = m.gap2char
getindex{C}(m::ScoreModel{C}, x::C,         y::C         ) = m.score[convert(UInt8,x)+1,convert(UInt8,y)+1]

setindex!{C}(m::ScoreModel{C}, s::Real,  ::Type{GAP}              ) = m.char2gap = m.gap2char = s
setindex!{C}(m::ScoreModel{C}, s::Real,  ::C,          ::Type{GAP}) = m.char2gap = s
setindex!{C}(m::ScoreModel{C}, s::Real,  ::Type{GAP},  ::C        ) = m.gap2char = s
setindex!{C}(m::ScoreModel{C}, s::Real, x::C,         y::C        ) = m.score[convert(UInt8,x)+1,convert(UInt8,y)+1] = s

# TODO: PAM, BLOSUM, etc.

# cost model is about a minimizing problem
abstract AbstractCostModel{T<:Real} <: AlignmentModel

# a gap and a mismatch costs 1; otherwise 0
immutable UnitCostModel{T} <: AbstractCostModel{T}; end
const UnitCost = UnitCostModel{Int}()

getindex{T}(::UnitCostModel{T},  ::Character,  ::Type{GAP}) = T(1)
getindex{T}(::UnitCostModel{T},  ::Type{GAP},  ::Character) = T(1)
getindex{T}(::UnitCostModel{T}, x::Character, y::Character) = ifelse(x === y, T(0), T(1))

minimum_indel_cost{T}(::UnitCostModel{T}) = T(1)

# arbitrary cost model
type CostModel{C<:Character,T} <: AbstractCostModel{T}
    cost::Matrix{T}
    char2gap::T
    gap2char::T
    function CostModel(alphabetsize)
        new(zeros(T, alphabetsize, alphabetsize), 0, 0)
    end
end

getindex{C}(m::CostModel{C},  ::C,          ::Type{GAP}) = m.char2gap
getindex{C}(m::CostModel{C},  ::Type{GAP},  ::C        ) = m.gap2char
getindex{C}(m::CostModel{C}, x::C,         y::C        ) = m.cost[convert(UInt8,x)+1,convert(UInt8,y)+1]

setindex!{C}(m::CostModel{C}, c::Real,  ::Type{GAP}              ) = m.char2gap = m.gap2char = c
setindex!{C}(m::CostModel{C}, c::Real,  ::C,          ::Type{GAP}) = m.char2gap = c
setindex!{C}(m::CostModel{C}, c::Real,  ::Type{GAP},  ::C        ) = m.gap2char = c
setindex!{C}(m::CostModel{C}, c::Real, x::C,         y::C        ) = m.cost[convert(UInt8,x)+1,convert(UInt8,y)+1] = c
