# Needleman-Wunsch algorithm
# --------------------------
#
# Naive Needleman-Wunsch algorithm.
#
# * type: global-global alignment
# * complexity: O(m*n)
# * space: O(m*n) for alignment, O(m) for score/distance
#
# references:
# * Needleman, Saul B., and Christian D. Wunsch. "A general method applicable to the search for similarities in the amino acid sequence of two proteins." Journal of molecular biology 48.3 (1970): 443-453.
# * Sankoff, David. "Matching sequences under deletion/insertion constraints." Proceedings of the National Academy of Sciences 69.1 (1972): 4-6.


type NeedlemanWunsch{M<:Union{AbstractScoreModel,AbstractCostModel},T<:Real} <: PairwiseAlignmentAlgorithm
    model::M
    matrix::DPMatrix{T}
end

function call{T}(::Type{NeedlemanWunsch}, model::Union{AbstractScoreModel{T},AbstractCostModel{T}})
    NeedlemanWunsch(model, DPMatrix{T}())
end

function score!(nw::NeedlemanWunsch, a, p, m, b, q, n, linear_space=true)
    if linear_space
        linear_dp!(nw, a, p, m, b, q, n)
    else
        dp!(nw, a, p, m, b, q, n)
    end
    return nw.matrix[end,end]
end

function distance!(nw::NeedlemanWunsch, a, p, m, b, q, n, linear_space=true)
    if linear_space
        linear_dp!(nw, a, p, m, b, q, n)
    else
        dp!(nw, a, p, m, b, q, n)
    end
    return nw.matrix[end,end]
end

@generated function dp!{M}(nw::NeedlemanWunsch{M}, a, p, m, b, q, n)
    if M <: AbstractScoreModel
        f = :max
    elseif M <: AbstractCostModel
        f = :min
    else
        @assert false
    end
    quote
        model = nw.model
        mtx = nw.matrix
        fitsize!(mtx, m, n)
        mtx[0,0] = 0
        for i in 1:m
            mtx[i,0] = mtx[i-1,0] + model[a[i+p-1],GAP]
        end
        for j in 1:n
            mtx[0,j] = mtx[0,j-1] + model[GAP,b[j+q-1]]
            for i in 1:m
                mtx[i,j] = $f(
                    mtx[i-1,j-1] + model[a[i+p-1],b[j+q-1]],
                    mtx[i-1,j  ] + model[a[i+p-1],GAP     ],
                    mtx[i,  j-1] + model[GAP,     b[j+q-1]]
                )
            end
        end
        return nw
    end
end

@generated function linear_dp!{M}(nw::NeedlemanWunsch{M}, a, p, m, b, q, n)
    if M <: AbstractScoreModel
        f = :max
    elseif M <: AbstractCostModel
        f = :min
    else
        @assert false
    end
    quote
        model = nw.model
        vec = nw.matrix
        fitsize!(vec, m)
        vec[0] = 0
        for i in 1:m
            vec[i] = vec[i-1] + model[a[i+p-1],GAP]
        end
        for j in 1:n
            diag = vec[0]
            vec[0] += model[GAP,b[j+q-1]]
            for i in 1:m
                tmp = vec[i]
                vec[i] = $f(
                    diag     + model[a[i+p-1],b[j+q-1]],
                    vec[i-1] + model[a[i+p-1],GAP     ],
                    vec[i  ] + model[GAP,     b[j+q-1]]
                )
                diag = tmp
            end
        end
        return nw
    end
end
