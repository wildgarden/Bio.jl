function distance(a, b, cost::AbstractCostModel=UnitCost)
    m = length(a)
    n = length(b)
    mtx = AlignmentMatrix{Int}(Int, m, n)
    fill_matrix!(mtx, a, b, cost)
    return mtx[m,n]
end

function distance2(a, b, cost::AbstractCostModel=UnitCost)
    m = length(a)
    n = length(b)
    mtx = AlignmentMatrix{Int}(Int, m, n)
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
    return mtx[m,n]
end
