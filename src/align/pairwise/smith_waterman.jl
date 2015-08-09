# Smith-Waterman algorithm
# ------------------------
#
# Naive Smith-Waterman algorithm.
#
# type: local-local alignment
# complexity: O(m*n)
# space: O(m*n)


type SmithWaterman{T} <: PairwiseAlignmentAlgorithm
    matrix::DPMatrix{T}
    max_score::T
    max_score_loc::Tuple{Int,Int}
    function SmithWaterman(matrix)
        new(matrix)
    end
end

function call{T}(::Type{SmithWaterman{T}})
    SmithWaterman{T}(DPMatrix{T}())
end

function call{T}(::Type{SmithWaterman}, ::AbstractScoreModel{T})
    SmithWaterman{T}(DPMatrix{T}())
end

function score!(sw::SmithWaterman, a, p, m, b, q, n, score)
    align!(sw, a, p, m, b, q, n, score)
    return sw.max_score
end

function align!(sw::SmithWaterman, a, p, m, b, q, n, score::AbstractScoreModel)
    mtx = sw.matrix
    fitsize!(mtx, m, n)
    mtx[0,0] = 0
    max_score = mtx[0,0]
    max_score_loc = (0, 0)
    # assuming the GAP score is non-positiev
    for i in 1:m
        mtx[i,0] = 0
    end
    for j in 1:n
        mtx[0,j] = 0
        for i in 1:m
            mtx[i,j] = max(
                0,
                mtx[i-1,j-1] + score[a[i+p-1],b[j+q-1]],
                mtx[i-1,j  ] + score[a[i+p-1],GAP     ],
                mtx[i,  j-1] + score[GAP,     b[j+q-1]]
            )
            if mtx[i,j] > max_score
                max_score = mtx[i,j]
                max_score_loc = (i, j)
            end
        end
    end
    sw.max_score = max_score
    sw.max_score_loc = max_score_loc
    return sw
end
