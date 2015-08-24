module Bio

abstract FileFormat

include("ragel.jl")
include("seq/seq.jl")
include("align/align.jl")
include("intervals/intervals.jl")

end # module
