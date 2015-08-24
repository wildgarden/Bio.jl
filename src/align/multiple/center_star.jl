# Center-Star-Method algorithm
# ------------------------
#
# type: approximation multiple sequence alignment
# sum-of-pairs distance is at most twice that of the optimal multiple alignment
# complexity: O(k^2 * n^2)
#
# reference: "Gusfield, Dan. "Efficient methods for multiple sequence alignment with guaranteed error bounds." Bulletin of mathematical biology 55.1 (1993): 141-154.

type CenterStar{M <: Union{AbstractLinearGapModel, AbstractCostModel}, T <: MultipleSequenceAlignment} <: MultipleAlignmentAlgorithm
    model::M
    alignment::T
end

function call{T}(::Type{CenterStar}, model::Union{AbstractLinearGapModel{T},AbstractCostModel{T}})
    CenterStar(model, MultipleSequenceAlignment())
end


"""
merges a sequence into a multiple alignment

param: seq_align_pair[1] must be the reference sequence
"""
function mergeSequence!(msa::MultipleSequenceAlignment, seq_align_pair::Tuple)
    msa = msa.alignment

    if isempty(msa)
        @assert length(seq_align_pair[1]) == length(seq_align_pair[2])
        push!(msa, seq_align_pair[1])
        push!(msa, seq_align_pair[2])
        # nothing more to do
        return
    end

    seq = seq_align_pair[2]

    pos = 1
    while length(msa[1]) >= pos && length(seq_align_pair[1]) >= pos
       # compare center sequences
        char_msa = msa[1][pos]
        char_seq = seq_align_pair[1][pos]
        if char_msa != char_seq
            if char_msa == '-' # GAP
                seq = insertGap(seq, pos)
            elseif char_seq == '-' # GAP
                insertGap!(msa, pos)
            end
        end
        pos += 1
    end

    # fill up to full length
    seq_len = length(seq)
    msa_len = length(msa[1])
    max_len = max(msa_len, seq_len)

    if  msa_len < max_len
        insertGap!(msa, max_len)
    elseif seq_len < max_len
        diff = max_len - seq_len
        seq = string(seq, "-"^diff)
    end

    push!(msa, seq)
end


"extends the given sequences by inserting a gap at the specified position"
function insertGap!(sequences::Array{ASCIIString, 1}, pos)
    for (i, seq) in enumerate(sequences)
        sequences[i] = insertGap(seq, pos)
    end
end


"inserts a gap at the given position"
function insertGap(sequence::String, pos)
    return string(sequence[1:pos - 1], '-', sequence[pos:end])
end


"summarizes the distances of a symmetric matrix"
function summarizedistances(distances)
    costs = []
    for i in 1:size(distances, 2)
        dist = 0
        for j in 1:size(distances, 2)
            if j < i
               dist += distances[j,i]
            else
                dist += distances[i,j]
            end
        end
        push!(costs, dist)
    end
    return costs
end


"""
Aligns multiple sequences.
Currently only cost-models are supported.
"""
function align(cs::CenterStar, sequences::Array)
    @assert typeof(cs.model) <: AbstractCostModel

    # Basic algorithm
    # 1. Find D(S_i, S_j) for all i, j
    # 2. Find center sequence with minimum distance ( sum from i to k over D(S_c, S_i))
    # 3. for every S_i element {S} - S_c , choose an optimal alignment between S_i and S_c
    # 4. Introduce spaces into S_c so that the multiple alignment M satisfies the alignments found in step 3

    # Step 1.
    msa          = cs.alignment
    len_seqs     = length(sequences)
    mtx_pw_align = Array{PairwiseAlignmentAlgorithm}(len_seqs, len_seqs)
    distances    = zeros(Real, len_seqs, len_seqs)

    # consider only the diagonal matrix, due to symmetry
    for i in 1:len_seqs
        for j in (i+1):len_seqs
            nw    = NeedlemanWunsch(cs.model)
            seq_a = sequences[i]
            seq_b = sequences[j]
            distances[i, j] = distance!(nw, seq_a, 1, length(seq_a), seq_b, 1, length(seq_b), false)
            mtx_pw_align[i, j] = nw
        end
    end

    # Step 2.
    costs = summarizedistances(distances)
    (mincost, idx_center) = findmin(costs)

    # Step 3.
    for col in 1:size(mtx_pw_align, 2)
        if col != idx_center
            optimal_alignment = chooseAlignment(mtx_pw_align, sequences, idx_center, col)
            mergeSequence!(msa, optimal_alignment) # step 4.
        end
    end

    return msa
end

"chooses an optimal pairwise alignment with respect to an symmetric alignment matrix"
function chooseAlignment(mtx, sequences, row, col)
    if col < row
        nw = mtx[col, row]
    else
        nw = mtx[row, col]
    end

    return backtrackSingle(nw, sequences[row], 1, sequences[col], 1)
end
