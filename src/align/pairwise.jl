# matrix used for dynamic programming
# O(m*n) space
type AlignmentMatrix{T<:Real} <: AbstractMatrix{T}
    nrows::Int
    ncols::Int
    matrix::Matrix{T}
    function AlignmentMatrix{T}(::Type{T}, nrows::Integer, ncols::Integer)
        resize!(new(), nrows, ncols)
    end
end

call{T}(::Type{AlignmentMatrix{T}}, nrows::Integer, ncols::Integer) = AlignmentMatrix{T}(T, nrows, ncols)

# index is 0-based
getindex(m::AlignmentMatrix, i::Integer, j::Integer) = m.matrix[i+1,j+1]
setindex!(m::AlignmentMatrix, x, i::Integer, j::Integer) = m.matrix[i+1,j+1] = x

size(m::AlignmentMatrix) = (m.nrows, m.ncols)

function size(m::AlignmentMatrix, d::Integer)
    if d == 1
        return m.nrows
    elseif d == 2
        return m.ncols
    end
    return 1
end

function resize!{T}(mtx::AlignmentMatrix{T}, nrows::Integer, ncols::Integer)
    mtx.matrix = Array{T}(nrows + 1, ncols + 1)
    mtx.nrows = nrows
    mtx.ncols = ncols
    mtx
end

function fitsize!(mtx::AlignmentMatrix, m, n)
    if mtx.nrows < m || mtx.ncols < n
        resize!(mtx, m, n)
    end
    mtx.nrows = m
    mtx.ncols = n
    return mtx
end

function empty!{T}(mtx::AlignmentMatrix{T})
    resize!(mtx, 0, 0)
end

# O(min(m,n)) space
type AlignmentVector{T<:Real} <: AbstractVector{T}
    nrows::Int
    vector::Vector{T}
    function AlignmentVector{T}(::Type{T}, nrows::Integer)
        resize!(new(), nrows)
    end
end

call{T}(::Type{AlignmentVector{T}}, nrows::Integer) = AlignmentVector{T}(T, nrows)

getindex(v::AlignmentVector, i::Integer) = v.vector[i+1]
setindex!(v::AlignmentVector, x, i::Integer) = v.vector[i+1] = x

length(v::AlignmentVector) = v.nrows

function resize!{T}(v::AlignmentVector{T}, m)
    if isdefined(v, :vector)
        resize!(v.vector, m + 1)
    else
        v.vector = Vector{T}(m + 1)
    end
    v.nrows = m
    v
end

function fitsize!(v::AlignmentVector, m)
    if v.nrows ≥ m
        v.nrows = m
    else
        resize!(v, m)
    end
    return v
end

function empty!(v::AlignmentMatrix)
    resize!(v, 0)
end


# global algorithms

# Conventions of variables:
#   a: sequence A
#   b: sequence B
#   m: aligned length of sequence A
#   n: aligned length of sequence B
#   p: start position of sequence A
#   q: start position of sequence B
# This means a[p:p+m-1] and b[q:q+n-1] are aligned to each other.

abstract PairwiseAlignmentAlgorithm

# Needleman-Wunsch
immutable NaiveDP <: PairwiseAlignmentAlgorithm; end

function fill_matrix!(mtx::AlignmentMatrix, a, p::Int, m::Int, b, q::Int, n::Int, cost::AbstractCostModel, ::Type{NaiveDP})
    fitsize!(mtx, m, n)
    mtx[0,0] = 0
    for i in 1:m
        mtx[i,0] = mtx[i-1,0] + cost[a[i+p-1],GAP]
    end
    for j in 1:n
        mtx[0,j] = mtx[0,j-1] + cost[GAP,b[j+q-1]]
    end
    for j in 1:n
        for i in 1:m
            mtx[i,j] = min(
                mtx[i-1,j-1] + cost[a[i+p-1],b[j+q-1]],
                mtx[i-1,j  ] + cost[a[i+p-1],GAP     ],
                mtx[i,  j-1] + cost[GAP,     b[j+q-1]]
            )
        end
    end
    return mtx
end

function fill_vector!(vec::AlignmentVector, a, p::Int, m::Int, b, q::Int, n::Int, cost::AbstractCostModel, ::Type{NaiveDP})
    @assert m ≤ n
    fitsize!(vec, m)
    vec[0] = 0
    for i in 1:m
        vec[i] = vec[i-1] + cost[a[i+p-1],GAP]
    end
    for j in 1:n
        diag = vec[0]
        vec[0] += cost[GAP,b[j+q-1]]
        for i in 1:m
            tmp = vec[i]
            vec[i] = min(
                diag     + cost[a[i+p-1],b[j+q-1]],
                vec[i-1] + cost[a[i+p-1],GAP     ],
                vec[i  ] + cost[GAP,     b[j+q-1]]
            )
            diag = tmp
        end
    end
    return vec
end

immutable ShortDetourDP <: PairwiseAlignmentAlgorithm; end

immutable AbberationError <: Exception; end

# fill cells within the "diagonal zone"
function fill_matrix!{T}(mtx::AlignmentMatrix{T}, a, p::Int, m::Int, b, q::Int, n::Int, t::T, cost::AbstractCostModel, ::Type{ShortDetourDP})
    fitsize!(mtx, m, n)
    # TODO: remove this restriction
    @assert m ≤ n
    d = n - m
    Δ = minimum_indel_cost(cost)
    if t < d * Δ
        throw(AbberationError())
    end
    mtx[0,0] = zero(T)
    # the diagonal zone is [-x..n-m+x], where x = ceil(t/(2Δ) - (n - m)/2)
    x = ceil(Int, t / (2Δ) - d / 2)
    for i in 1:x
        mtx[i,0] = mtx[i-1,0] + cost[a[i+p-1],GAP]
    end
    for j in 1:d+x
        mtx[0,j] = mtx[0,j-1] + cost[GAP,b[j+q-1]]
    end
    for j in 1:n
        l = max(1, j - (d + x))
        u = min(j + x, m)
        min_cost = typemax(T)
        for i in l:u
            c = mtx[i-1,j-1] + cost[a[i+p-1],b[j+q-1]]
            if i != l
                c = min(c, mtx[i-1,j] + cost[a[i+p-1],GAP])
            end
            if i != u
                c = min(c, mtx[i,j-1] + cost[GAP,b[j+q-1]])
            end
            min_cost = min(min_cost, c)
            mtx[i,j] = c
        end
        if min_cost > t
            throw(AbberationError())
        end
    end
    if mtx[m,n] > t
        throw(AbberationError())
    end
    return mtx
end

type Alignment
    len::Int
    data::Vector{UInt8}
end

Alignment() = Alignment(0, UInt8[])

const Mat = 0x00
const Ins = 0x01
const Del = 0x02

const code_char = ['M', 'I', 'D']

function show(io::IO, aln::Alignment)
    println(io, "Alignment")
    print(' ')
    for i in 1:endof(aln)
        op = aln[i]
        if op & 0x80 != 0
            @assert op == (Mat | 0x80)
            print(io, "M̂")
        else
            print(io, code_char[op+1])
        end
    end
end

length(aln::Alignment) = aln.len
endof(aln::Alignment) = aln.len

getindex(aln::Alignment, i::Integer) = aln.data[i]

function push!(aln::Alignment, op)
    if aln.len + 1 > length(aln.data)
        resize!(aln.data, aln.len + 1)
    end
    aln.data[aln.len+=1] = op
    aln
end

function reverse!(aln::Alignment)
    reverse!(aln.data, 1, aln.len)
    aln
end

function reset!(aln::Alignment)
    aln.len = 0
    aln
end

function traceback!(aln::Alignment, mtx::AlignmentMatrix, a, b, cost::AbstractCostModel)
    i, j = size(mtx)
    reset!(aln)
    while i > 0 || j > 0
        if i ≥ 1 && j ≥ 1 && mtx[i,j] == mtx[i-1,j-1] + cost[a[i],b[j]]
            push!(aln, Mat | (a[i] == b[j] ? 0x00 : 0x80))
            i -= 1
            j -= 1
        elseif i ≥ 1 && mtx[i,j] == mtx[i-1,j] + cost[a[i],GAP]
            push!(aln, Del)
            i -= 1
        elseif j ≥ 1 && mtx[i,j] == mtx[i,j-1] + cost[GAP,b[j]]
            push!(aln, Ins)
            j -= 1
        else
            @assert false
        end
    end
    reverse!(aln)
    return aln
end

function traceback(mtx::AlignmentMatrix, a, b, cost::AbstractCostModel)
    aln = Alignment()
    traceback!(aln, mtx, a, b, cost)
end

function align(a, b, cost::AbstractCostModel=UnitCost)
    m = length(a)
    n = length(b)
    mtx = AlignmentMatrix{Int}(m, n)
    fill_matrix!(mtx, a, 1, m, b, 1, n, cost, NaiveDP)
    traceback(mtx, a, b, cost)
end


# `a` and `b` are something like a sequence
function distance{A<:PairwiseAlignmentAlgorithm}(a, b, cost::AbstractCostModel=UnitCost, alg::Type{A}=NaiveDP)
    return distance(a, b, cost, alg)
end

function distance{A<:PairwiseAlignmentAlgorithm}(a, b, cost::AbstractCostModel, ::Type{A})
    mtx = AlignmentMatrix{Int}(length(a), length(b))
    return distance!(mtx, a, b, cost, A)
end

function distance!(mtx::AlignmentMatrix, a, b, cost::AbstractCostModel, ::Type{NaiveDP})
    m = length(a)
    n = length(b)
    fill_matrix!(mtx, a, 1, m, b, 1, n, cost, NaiveDP)
    return mtx[end,end]
end

function distance!(vec::AlignmentVector, a, b, cost::AbstractCostModel, ::Type{NaiveDP})
    m = length(a)
    n = length(b)
    fill_vector!(vec, a, 1, m, b, 1, n, cost, NaiveDP)
    return vec[end]
end

function distance!(mtx::AlignmentMatrix, a, b, cost::AbstractCostModel, ::Type{ShortDetourDP})
    m = length(a)
    n = length(b)
    t = 1
    while true
        try
            fill_matrix!(mtx, a, 1, m, b, 1, n, t, cost, ShortDetourDP)
        catch ex
            if isa(ex, AbberationError)
                t *= 2
                continue
            end
            rethrow()
        end
        break
    end
    return mtx[end,end]
end
