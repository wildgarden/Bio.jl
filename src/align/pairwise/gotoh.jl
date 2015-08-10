# Gotoh's algorithm
# -----------------
#
# Gotoh's alignment algorithm for affine gap scores.
#
# Gotoh, Osamu. "An improved algorithm for matching biological sequences." Journal of molecular biology 162.3 (1982): 705-708.

type Gotoh{M<:AffineScoreModel,T<:Real} <: PairwiseAlignmentAlgorithm
    model::M
    sm::DPMatrix{T}
    sg::DPMatrix{T}
end

function call{T,C}(::Type{Gotoh}, model::AffineScoreModel{C,T})
    Gotoh(model, DPMatrix{T}(), DPMatrix{T}())
end

function score!(go::Gotoh, a, p, m, b, q, n)
    dp!(go, a, p, m, b, q, n)
    return max(go.sm[end,end], go.sg[end,end])
end

function dp!{M,T}(go::Gotoh{M,T}, a, p, m, b, q, n)
    score = go.model
    α = score.α
    β = score.β
    sm = go.sm
    sg = go.sg
    fitsize!(sm, m, n)
    fitsize!(sg, m, n)
    sm[0,0] = 0
    sg[0,0] = 0
    for i in 1:m
        sm[i,0] = typemin(T)
        sg[i,0] = -α - β * (i - 1)
    end
    for j in 1:n
        sm[0,j] = typemin(T)
        sg[0,j] = -α - β * (j - 1)
    end
    for j in 1:n
        for i in 1:m
            sm[i,j] = max(
                sm[i-1,j-1] + score[a[i],b[j]],
                sg[i-1,j-1] + score[a[i],b[j]]
            )
            sg[i,j] = max(
                sm[i-1,j] - α,
                sg[i-1,j] - β,
                sm[i,j-1] - α,
                sg[i,j-1] - β
            )
        end
    end
    return go
end
