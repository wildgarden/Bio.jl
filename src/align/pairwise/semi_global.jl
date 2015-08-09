# Semi-global alignment
# ---------------------
#
# A slightly modified version of Needleman-Wunsch algorthim.
# TODO: more specific name?
#
# * type: global-local alignment
# * complexity: O(m*n)
# * space: O(m*n)

type SemiGlobal{T} <: PairwiseAlignmentAlgorithm
    matrix::DPMatrix{T}
    min_cost::T
    function SemiGlobal(matrix)
        new(matrix)
    end
end

function call{T}(::Type{SemiGlobal{T}})
    SemiGlobal{T}(DPMatrix{T}())
end

function call{T}(::Type{SemiGlobal}, ::AbstractCostModel{T})
    SemiGlobal{T}(DPMatrix{T}())
end

function distance!(sg::SemiGlobal, a, p, m, b, q, n, cost::AbstractCostModel)
    dp!(sg, a, p, m, b, q, n, cost)
    return sg.min_cost
end

function dp!(sg::SemiGlobal, a, p, m, b, q, n, cost::AbstractCostModel)
    mtx = sg.matrix
    fitsize!(mtx, m, n)
    mtx[0,0] = 0
    for i in 1:m
        mtx[i,0] = mtx[i-1,0] + cost[a[i+p-1],GAP]
    end
    min_cost = mtx[m,0]
    for j in 1:n
        mtx[0,j] = 0
        for i in 1:m
            mtx[i,j] = min(
                mtx[i-1,j-1] + cost[a[i+p-1],b[j+q-1]],
                mtx[i-1,j  ] + cost[a[i+p-1],GAP     ],
                mtx[i,  j-1] + cost[GAP,     b[j+q-1]]
            )
        end
        min_cost = min(mtx[m,j], min_cost)
    end
    sg.min_cost = min_cost
    return sg
end
