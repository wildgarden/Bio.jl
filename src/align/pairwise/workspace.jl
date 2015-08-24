# 0-based matrix for dynamic programming
type DPMatrix{T} <: AbstractArray{T,2}
    nrows::Int
    ncols::Int
    data::Array{T,2}
    function DPMatrix{T}(::Type{T})
        return new(0, 0)
    end
end

function call{T}(::Type{DPMatrix{T}})
    return DPMatrix{T}(T)
end

function call{T}(::Type{DPMatrix{T}}, nrows::Integer, ncols::Integer)
    m = DPMatrix{T}(T)
    resize!(m, nrows, ncols)
end

size(m::DPMatrix) = (m.nrows, m.ncols)

getindex(m::DPMatrix, i::Integer, j::Integer, k::Integer) = m.data[i+1,j+1,k+1]
setindex!(m::DPMatrix, x, i::Integer, j::Integer, k::Integer) = m.data[i+1,j+1,k+1] = x

getindex(m::DPMatrix, i::Integer, j::Integer) = m.data[i+1,j+1]
setindex!(m::DPMatrix, x, i::Integer, j::Integer) = m.data[i+1,j+1] = x

getindex(m::DPMatrix, i::Integer) = m.data[i+1,1]
setindex!(m::DPMatrix, x, i::Integer) = m.data[i+1,1] = x

function resize!{T}(m::DPMatrix{T}, nrows, ncols)
    m.data = Array{T}(nrows + 1, ncols + 1)
    m.nrows = nrows
    m.ncols = ncols
    return m
end

function fitsize!(m::DPMatrix, nrows, ncols=0)
    if m.nrows < nrows || m.ncols < ncols
        resize!(m, nrows, ncols)
    else
        m.nrows = nrows
        m.ncols = ncols
    end
    return m
end
