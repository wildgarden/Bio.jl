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

immutable SubstitutionMatrix{C<:Character,T}
    score::Matrix{T}
end

call{C,T}(::Type{SubstitutionMatrix{C}}, score::Matrix{T}) = SubstitutionMatrix{C,T}(score)

getindex{C}(m::SubstitutionMatrix{C}, x::C, y::C) = m.score[convert(UInt8,x)+1,convert(UInt8,y)+1]

function show{T}(io::IO, sm::SubstitutionMatrix{AminoAcid,T})
    aarange = 0x00:0x14
    println(io, "SubstitutionMatrix{AminoAcid,", T, "}")
    print(io, "  ")
    for j in aarange
        b = convert(AminoAcid, j)
        print(io, "  ", b)
    end
    println()
    for i in aarange
        a = convert(AminoAcid, i)
        print(io, " ", a)
        for j in aarange
            b = convert(AminoAcid, j)
            @printf io "%3d" sm[a,b]
        end
        println(io)
    end
end

function parse_blosum_blast(io::IO)
    matrix = Matrix{Int}(21, 21)
    columns = Dict{Int,AminoAcid}()
    i = 0
    for line in eachline(io)
        line[1] == '#' && continue
        values = split(strip(line), r"\s+")
        if i == 0
            # read columns
            for (j, val) in enumerate(values)
                aa = tryparse(AminoAcid, val)
                if !isnull(aa)
                    columns[j] = get(aa)
                end
            end
        else
            for (j, val) in enumerate(values)
                if haskey(columns, i) && haskey(columns, j)
                    i′ = convert(UInt8, columns[i])
                    j′ = convert(UInt8, columns[j])
                    matrix[i′+1,j′+1] = parse(Int, val)
                end
            end
        end
        i += 1
    end
    return SubstitutionMatrix{AminoAcid}(matrix)
end

const BLOSUM62 = open(parse_blosum_blast, joinpath(dirname(@__FILE__), "blosum62.bla"))

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
