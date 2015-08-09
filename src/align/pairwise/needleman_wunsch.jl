# Needleman-Wunsch algorithm
# --------------------------
#
# Naive Needleman-Wunsch algorithm.
#
# * type: global-global alignment
# * complexity: O(m*n)
# * space: O(m*n) for alignment, O(m) for score/distance


type NeedlemanWunsch{T} <: PairwiseAlignmentAlgorithm
    matrix::DPMatrix{T}
end

function call{T}(::Type{NeedlemanWunsch{T}})
    NeedlemanWunsch(DPMatrix{T}())
end

function call{T}(::Type{NeedlemanWunsch}, ::Union{AbstractScoreModel{T},AbstractCostModel{T}})
    NeedlemanWunsch(DPMatrix{T}())
end

function score!(nw::NeedlemanWunsch, a, p, m, b, q, n, score::AbstractScoreModel, linear_space=true)
    if linear_space
        linear_dp!(nw, a, p, m, b, q, n, score)
    else
        dp!(nw, a, p, m, b, q, n, score)
    end
    return nw.matrix[end,end]
end

function distance!(nw::NeedlemanWunsch, a, p, m, b, q, n, cost::AbstractCostModel, linear_space=true)
    if linear_space
        linear_dp!(nw, a, p, m, b, q, n, cost)
    else
        dp!(nw, a, p, m, b, q, n, cost)
    end
    return nw.matrix[end,end]
end

function dp!(nw::NeedlemanWunsch, a, p, m, b, q, n, score::AbstractScoreModel)
    mtx = nw.matrix
    fitsize!(mtx, m, n)
    mtx[0,0] = 0
    for i in 1:m
        mtx[i,0] = mtx[i-1,0] + score[a[i+p-1],GAP]
    end
    for j in 1:n
        mtx[0,j] = mtx[0,j-1] + score[GAP,b[j+q-1]]
        for i in 1:m
            mtx[i,j] = max(
                mtx[i-1,j-1] + score[a[i+p-1],b[j+q-1]],
                mtx[i-1,j  ] + score[a[i+p-1],GAP     ],
                mtx[i,  j-1] + score[GAP,     b[j+q-1]]
            )
        end
    end
    return nw
end

function linear_dp!(nw::NeedlemanWunsch, a, p, m, b, q, n, score::AbstractScoreModel)
    # TODO: inefficient when m >> n
    vec = nw.matrix
    fitsize!(vec, m)
    vec[0] = 0
    for i in 1:m
        vec[i] = vec[i-1] + score[a[i+p-1],GAP]
    end
    for j in 1:n
        diag = vec[0]
        vec[0] += score[GAP,b[j+q-1]]
        for i in 1:m
            tmp = vec[i]
            vec[i] = max(
                diag     + score[a[i+p-1],b[j+q-1]],
                vec[i-1] + score[a[i+p-1],GAP     ],
                vec[i  ] + score[GAP,     b[j+q-1]]
            )
            diag = tmp
        end
    end
    return nw
end

function dp!(nw::NeedlemanWunsch, a, p, m, b, q, n, cost::AbstractCostModel)
    mtx = nw.matrix
    fitsize!(mtx, m, n)
    mtx[0,0] = 0
    for i in 1:m
        mtx[i,0] = mtx[i-1,0] + cost[a[i+p-1],GAP]
    end
    for j in 1:n
        mtx[0,j] = mtx[0,j-1] + cost[GAP,b[j+q-1]]
        for i in 1:m
            mtx[i,j] = min(
                mtx[i-1,j-1] + cost[a[i+p-1],b[j+q-1]],
                mtx[i-1,j  ] + cost[a[i+p-1],GAP     ],
                mtx[i,  j-1] + cost[GAP,     b[j+q-1]]
            )
        end
    end
    return nw
end

function linear_dp!(nw::NeedlemanWunsch, a, p, m, b, q, n, cost::AbstractCostModel)
    # TODO: inefficient when m >> n
    vec = nw.matrix
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
    return nw
end
