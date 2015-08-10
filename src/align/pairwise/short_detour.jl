# Short detour algorithm
# ----------------------
#
# Ukkonen's short detour algorithm.
#
# * type: global-global alignment
# * complexity: O(d*m), where d = edit distance
# * space: O(m*n) (TODO: improve)


type ShortDetour{M<:AbstractCostModel,T<:Real} <: PairwiseAlignmentAlgorithm
    model::M
    matrix::DPMatrix{T}
end

function call{T}(::Type{ShortDetour}, model::AbstractCostModel{T})
    ShortDetour(model, DPMatrix{T}())
end

immutable AbberationError <: Exception; end

function distance!(sd::ShortDetour, a, p, m, b, q, n)
    t = 1
    while true
        try
            dp!(sd, a, p, m, b, q, n, t)
        catch ex
            if isa(ex, AbberationError)
                # double threshold of aberration
                t *= 2
                continue
            end
            rethrow()
        end
        break
    end
    return sd.matrix[end,end]
end

function dp!{M,T}(sd::ShortDetour{M,T}, a, p, m, b, q, n, t::T)
    cost = sd.model
    mtx = sd.matrix
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
            if i != u || j - 1 + x ≥ m
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
    return sd
end
