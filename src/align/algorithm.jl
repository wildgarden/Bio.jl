# matrix used for dynamic programming
type AlignmentMatrix{T<:Real}
    nrows::Int
    ncols::Int
    matrix::Matrix{T}
    function AlignmentMatrix{T}(::Type{T}, nrows, ncols)
        new(nrows, ncols, Array{T}(nrows + 1, ncols + 1))
    end
end

# index is 0-based
getindex(m::AlignmentMatrix, i::Integer, j::Integer) = m.matrix[i+1,j+1]
setindex!(m::AlignmentMatrix, x, i::Integer, j::Integer) = m.matrix[i+1,j+1] = x

function size(m::AlignmentMatrix, d::Integer)
    if d == 1
        return m.nrows
    elseif d == 2
        return m.ncols
    end
    return 1
end

function fill_matrix!(mtx::AlignmentMatrix, a, b, cost::AbstractCostModel)
    m = length(a)
    n = length(b)
    @assert m ≤ size(mtx, 1)
    @assert n ≤ size(mtx, 2)
    mtx[0,0] = 0
    for i in 1:m
        mtx[i,0] = mtx[i-1,0] + cost[a[i],GAP]
    end
    for j in 1:n
        mtx[0,j] = mtx[0,j-1] + cost[GAP,b[j]]
    end
    for j in 1:n
        for i in 1:m
            mtx[i,j] = min(
                mtx[i-1,j-1] + cost[a[i],b[j]],
                mtx[i-1,j  ] + cost[a[i],GAP ],
                mtx[i,  j-1] + cost[GAP, b[j]]
            )
        end
    end
    return mtx
end

immutable AbberationError <: Exception; end

function fill_matrix!{T}(mtx::AlignmentMatrix{T}, a, b, t::T, cost::AbstractCostModel)
    m = length(a)
    n = length(b)
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
        mtx[i,0] = mtx[i-1,0] + cost[a[i],GAP]
    end
    for j in 1:d+x
        mtx[0,j] = mtx[0,j-1] + cost[GAP,b[j]]
    end
    for j in 1:n
        l = max(1, j - (d + x))
        u = min(j + x, m)
        min_cost = typemax(T)
        for i in l:u
            c = mtx[i-1,j-1] + cost[a[i],b[j]]
            if i != l
                c = min(c, mtx[i-1,j] + cost[a[i],GAP])
            end
            if i != u
                c = min(c, mtx[i,j-1] + cost[GAP,b[j]])
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
