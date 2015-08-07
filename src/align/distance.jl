# `a` and `b` are something like a sequence
function distance{A<:AlignmentAlgorithm}(a, b, cost::AbstractCostModel=UnitCost, alg::Type{A}=NaiveDP)
    return distance(a, b, cost, alg)
end

function distance{A<:AlignmentAlgorithm}(a, b, cost::AbstractCostModel, ::Type{A})
    mtx = AlignmentMatrix{Int}(length(a), length(b))
    return distance!(mtx, a, b, cost, A)
end

function distance!(mtx::AlignmentMatrix, a, b, cost::AbstractCostModel, ::Type{NaiveDP})
    fill_matrix!(mtx, a, b, cost, NaiveDP)
    return mtx[end,end]
end

function distance!(mtx::AlignmentMatrix, a, b, cost::AbstractCostModel, ::Type{ShortDetourDP})
    t = 1
    while true
        try
            fill_matrix!(mtx, a, b, t, cost, ShortDetourDP)
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
