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


type NeedlemanWunsch{M<:Union{AbstractLinearGapModel,AbstractCostModel},T<:Real, S} <: PairwiseAlignmentAlgorithm
    model::M
    matrix::DPMatrix{T}
    trace_matrix::DPMatrix{S}
end

function call{T}(::Type{NeedlemanWunsch}, model::Union{AbstractLinearGapModel{T},AbstractCostModel{T}})
    NeedlemanWunsch(model, DPMatrix{T}(), DPMatrix{Array{Array{Int64,1}}}())
end

# Conventions of variables:
#   a: sequence A
#   b: sequence B
#   m: aligned length of sequence A
#   n: aligned length of sequence B
#   p: start position of sequence A
#   q: start position of sequence B
# This means a[p:p+m-1] and b[q:q+n-1] are aligned to each other.

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
    if M <: AbstractLinearGapModel
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

        mtx_path = nw.trace_matrix
        fitsize!(mtx_path, m, n)

        fill!(mtx_path, collect([]))
        for i in 0:m
            for j in 0:n
                mtx_path[i,j,0] = []
            end
        end
        push!(mtx_path[0,0,0],[0,0])

        for i in 1:m
            mtx[i,0] = mtx[i-1,0] + model[a[i+p-1],GAP]
            push!(mtx_path[i,0], [i-1,0])
        end

        for j in 1:n
            mtx[0,j] = mtx[0,j-1] + model[GAP,b[j+q-1]]
            push!(mtx_path[0,j], [0,j-1])
            for i in 1:m
                diagonal = mtx[i-1,j-1] + model[a[i+p-1],b[j+q-1]]
                top      = mtx[i-1,j  ] + model[a[i+p-1],GAP     ]
                left     = mtx[i,  j-1] + model[GAP,     b[j+q-1]]

                mtx[i,j] = $f(diagonal, top, left)

                for (idx, val) in enumerate([diagonal, top, left])
                    if val == mtx[i, j]
                        if idx == 1 #diag
                            push!(mtx_path[i,j], [i-1,j-1])
                        elseif idx == 2 # top
                            push!(mtx_path[i,j], [i-1,j])
                        else # left
                            push!(mtx_path[i,j], [i,j-1])
                        end
                    end
                end
            end
        end

        return nw
    end
end


"""
finds one optimal alignment for the given sequences
please note: for the given sequences more than optimal alignment may exist. This function returns only one.
"""
function backtrackSingle(nw::NeedlemanWunsch, seq_a, p::Int, seq_b, q::Int)
    align_a = ""
    align_b = ""

    i = length(seq_a)
    j = length(seq_b)

    while i > 0 || j > 0
        if i > 0 && j > 0 && nw.matrix[i, j] == (nw.matrix[i - 1, j - 1] + nw.model[seq_a[i + p - 1], seq_b[j + q - 1]])
            align_a *= string(seq_a[i + p - 1])
            align_b *= string(seq_b[j + q - 1])
            i -= 1
            j -= 1
        elseif i > 0 && nw.matrix[i, j] == (nw.matrix[i - 1, j] + nw.model[seq_a[i + p - 1], GAP])
            align_a *= string(seq_a[i + p - 1])
            align_b *= string("-")
            i -= 1
        else
            align_a *= string("-")
            align_b *= string(seq_b[j + q - 1])
            j -= 1
        end
    end

    return (reverse(align_a), reverse(align_b))
end


"""
Recursive backtracking returns all optimal alignments.
Parameter `result` contains all optimals alignments when function returns. Other params are unchanged.
"""
function backtrackAll!(nw::NeedlemanWunsch, result, row, col, seq_a_orig, seq_b_orig, seq_a_aligned = "", seq_b_aligned = "")
    if row == 0 && col == 0
        push!(result,(reverse(seq_a_aligned), reverse(seq_b_aligned)))
    else
        if length(nw.trace_matrix[row,col]) > 0
            for (i,j) in nw.trace_matrix[row,col]
                if i == (row - 1) && j == col # top
                    seq_a = string(seq_a_aligned, seq_a_orig[row])
                    seq_b = string(seq_b_aligned, '-')
                elseif i == (row - 1) && j == (col -1) # diagonal
                    seq_a = string(seq_a_aligned, seq_a_orig[row])
                    seq_b = string(seq_b_aligned, seq_b_orig[col])
                else # left
                    seq_a = string(seq_a_aligned, '-')
                    seq_b = string(seq_b_aligned, seq_b_orig[col])
                end

                backtrackAll!(nw, result, i, j, seq_a_orig, seq_b_orig, seq_a, seq_b)
            end
        else
            error("empty path in $row, $col")
        end
    end
end

@generated function linear_dp!{M}(nw::NeedlemanWunsch{M}, a, p, m, b, q, n)
    if M <: AbstractLinearGapModel
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
