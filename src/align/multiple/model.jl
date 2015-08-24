abstract Alignment

# Stores sequences of multiple alignments
type MultipleSequenceAlignment
    alignment::Array{ASCIIString, 1}

    function MultipleSequenceAlignment(alignment::Array{ASCIIString, 1})
        return new(alignment)
    end

    function MultipleSequenceAlignment()
        return new(Array{ASCIIString, 1}())
    end
end

# look-up-table for color codes
char_to_color = Dict("A" => :yellow, "C" => :blue, "G" => :green, "T" => :red, "U" => :red, "-" => :magenta)

"shows the sequences of a multiple alignment"
function Base.show(io::IO, msa::MultipleSequenceAlignment)
    rows = length(msa.alignment)
    cols = isempty(msa.alignment) ? "0" : length(msa.alignment[1])
    write(io, "Multiple alignment with $rows rows and $cols cols", "\n")
    for seq in msa.alignment
        for c in seq
            color = get(char_to_color, string(c), :cyan)
            print_with_color(color, io, string(c))
        end
        write(io, "\n") # todo remove when typeof sequence is changed to SeqRecord
    end
end

getindex(m::MultipleSequenceAlignment, i::Integer) = m.alignment[i]
