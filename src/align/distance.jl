const default_cost = UnitCost

function distance(a, b, cost::AbstractCostModel=default_cost)
    mtx = AlignmentMatrix{Int}(Int, length(a), length(b))
    distance!(mtx, a, b, cost)
    return mtx[end,end]
end

function distance!(mtx::AlignmentMatrix, a, b, cost::AbstractCostModel=default_cost)
    fill_matrix!(mtx, a, b, cost)
    return mtx[end,end]
end

function distance2(a, b, cost::AbstractCostModel=default_cost)
    mtx = AlignmentMatrix{Int}(Int, length(a), length(b))
    distance2!(mtx)
    return mtx[end,end]
end

function distance2!(mtx::AlignmentMatrix, a, b, cost::AbstractCostModel=default_cost)
    t = 1
    while true
        try
            fill_matrix!(mtx, a, b, t, cost)
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
