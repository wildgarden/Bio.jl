

immutable BED <: FileFormat end

"Metadata for BED interval records"
immutable BEDMetadata
    name::Nullable{ASCIIString}
    score::Nullable{Int}
    thick_first::Nullable{Int}
    thick_last::Nullable{Int};
    item_rgb::Nullable{RGB{Float32}};
    block_count::Nullable{Int}
    block_sizes::Nullable{Vector{Int}}
    block_firsts::Nullable{Vector{Int}}

    function BEDMetadata(name=Nullable{ASCIIString}(),
                         score=Nullable{Int}(),
                         thick_first=Nullable{Int}(),
                         thick_last=Nullable{Int}(),
                         item_rgb=Nullable{Float32}(),
                         block_count=Nullable{Int}(),
                         block_sizes=Nullable{Vector{Int}}(),
                         block_firsts=Nullable{Vector{Int}}())
        return new(name, score, thick_first, thick_last, item_rgb,
                   block_count, block_sizes, block_firsts)
    end
end


function (==)(a::BEDMetadata, b::BEDMetadata)
    for name in fieldnames(BEDMetadata)
        aval = getfield(a, name)
        bval = getfield(b, name)
        if !((isnull(aval) && isnull(bval)) || get(aval) == get(bval))
            return false
        end
    end
    return true
end


"An `Interval` with associated metadata from a BED file"
typealias BEDInterval Interval{BEDMetadata}


module BEDParserImpl

import Bio.Intervals: Strand, STRAND_NA, BEDInterval, BEDMetadata
import Bio.Ragel
using Switch, Compat, Color
export BEDParser, takevalue!


%%{
    machine bed;

    action yield {
        yield = true
        # // fbreak causes will cause the pushmark action for the next seqname
        # // to be skipped, so we do it here
        Ragel.@pushmark!
        fbreak;
    }

    action count_line { input.state.linenum += 1 }
    action pushmark { Ragel.@pushmark! }

    action seqname     { input.seqname      = Ragel.@asciistring_from_mark! }
    action first       { input.first        = 1 + Ragel.@int64_from_mark! }
    action last        { input.last         = Ragel.@int64_from_mark! }
    action name        { input.name         = Nullable{ASCIIString}(Ragel.@asciistring_from_mark!) }
    action score       { input.score        = Ragel.@int64_from_mark! }
    action strand      { input.strand       = convert(Strand, Ragel.@char) }
    action thick_first { input.thick_first  = Ragel.@int64_from_mark! }
    action thick_last  { input.thick_last   = Ragel.@int64_from_mark! }
    action item_rgb_r  { input.red = input.green = input.blue = (Ragel.@int64_from_mark!) / 255.0 }
    action item_rgb_g  { input.green        = (Ragel.@int64_from_mark!) / 255.0 }
    action item_rgb_b  { input.blue         = (Ragel.@int64_from_mark!) / 255.0 }
    action item_rgb    { input.item_rgb     = RGB{Float32}(input.red, input.green, input.blue ) }
    action block_count { input.block_count  = Ragel.@int64_from_mark! }
    action block_size {
        if isnull(input.block_sizes)
            input.block_sizes = Array(Int, 0)
        end
        push!(get(input.block_sizes), (Ragel.@int64_from_mark!))
    }
    action block_first {
        if isnull(input.block_firsts)
            input.block_firsts = Array(Int, 0)
        end
        push!(get(input.block_firsts), (Ragel.@int64_from_mark!))
    }

    newline      = '\r'? '\n'     >count_line;
    hspace       = [ \t\v];
    blankline    = hspace* newline;

    seqname      = [ -~]*   >pushmark     %seqname;
    first        = digit+   >pushmark     %first;
    last         = digit+   >pushmark     %last;
    name         = [ -~]*   >pushmark     %name;
    score        = digit+   >pushmark     %score;
    strand       = [+\-\.?] >strand;
    thick_first  = digit+   >pushmark     %thick_first;
    thick_last   = digit+   >pushmark     %thick_last;

    item_rgb_r   = digit+   >pushmark     %item_rgb_r;
    item_rgb_g   = digit+   >pushmark     %item_rgb_g;
    item_rgb_b   = digit+   >pushmark     %item_rgb_b;
    item_rgb     = item_rgb_r (hspace* ',' hspace* item_rgb_g hspace* ',' hspace* item_rgb_b)? %item_rgb;

    block_count  = digit+   >pushmark    %block_count;

    block_size   = digit+   >pushmark    %block_size;
    block_sizes  = block_size (',' block_size)* ','?;

    block_first  = digit+   >pushmark    %block_first;
    block_firsts = block_first (',' block_first)* ','?;

    bed_entry = seqname '\t' first '\t' last (
                    '\t' name ( '\t' score ( '\t' strand ( '\t' thick_first (
                    '\t' thick_last ( '\t' item_rgb ( '\t' block_count (
                    '\t' block_sizes ( '\t' block_firsts )? )? )? )? )? )? )? )?
                    )?
                newline blankline*;

    main := blankline* (bed_entry %yield)*;
}%%


%% write data;


type BEDParser
    state::Ragel.State

    # intermediate values when parsing
    seqname::ASCIIString
    first::Int64
    last::Int64
    strand::Strand

    red::Float32
    green::Float32
    blue::Float32
    name::Nullable{ASCIIString}
    score::Nullable{Int}
    thick_first::Nullable{Int}
    thick_last::Nullable{Int};
    item_rgb::Nullable{RGB{Float32}};
    block_count::Nullable{Int}
    block_sizes::Nullable{Vector{Int}}
    block_firsts::Nullable{Vector{Int}}

    function BEDParser(input::Union(IO, String, Vector{Uint8}),
                       memory_map::Bool=false)
        %% write init;

        return new(Ragel.State(cs, input, memory_map),
                   "", 0, 0, STRAND_NA, 0.0, 0.0, 0.0,
                   Nullable{ASCIIString}(), Nullable{Int}(), Nullable{Int}(),
                   Nullable{Int}(), Nullable{RGB{Float32}}(), Nullable{Int}(),
                   Nullable{Vector{Int}}(), Nullable{Vector{Int}}())
    end
end


function Ragel.ragelstate(parser::BEDParser)
    return parser.state
end


function accept_state!(input::BEDParser, output::BEDInterval)
    output = input.nextitem
    input.nextitem = BEDInterval()
end


function takevalue!(input::BEDParser)
    value = BEDInterval(input.seqname, input.first, input.last, input.strand,
                        BEDMetadata(input.name, input.score, input.thick_first,
                                    input.thick_last, input.item_rgb,
                                    input.block_count, input.block_sizes,
                                    input.block_firsts))
    input.strand = STRAND_NA
    input.name = Nullable{ASCIIString}()
    input.score = Nullable{Int}()
    input.thick_first = Nullable{Int}()
    input.thick_last = Nullable{Int}()
    input.item_rgb = Nullable{RGB{Float32}}()
    input.block_count = Nullable{Int}()
    input.block_sizes = Nullable{Vector{Int}}()
    input.block_firsts = Nullable{Vector{Int}}()

    return value
end


Ragel.@generate_read_fuction("bed", BEDParser, BEDInterval,
    begin
        @inbounds begin
            %% write exec;
        end
    end,
    begin
        # TODO: If I'm going to do actual destructive parsing, I
        # need to be able to do some setup here.
        accept_state!(input, output)
    end)



end # module BEDParserImpl


# This inexplicably doesn't work, which is why I qualify BEDParser below.
#using BEDParserImpl

"An iterator over entries in a BED file or stream"
type BEDIterator <: IntervalStream{BEDMetadata}
    parser::BEDParserImpl.BEDParser
    nextitem::Nullable{BEDInterval}
end


"""
Parse a BED file.

# Arguments
  * `filename::String`: Path of the BED file.
  * `memory_map::Bool`: If true, attempt to memory map the file on supported
    platforms. (Default: `false`)

# Returns
An iterator over `BEDInterval`s contained in the file.
"""
function read(filename::String, ::Type{BED}; memory_map::Bool=false)
    it = BEDIterator(BEDParserImpl.BEDParser(filename, memory_map),
                       Nullable{BEDInterval}())
    return it
end


"""
Parse a BED file.

# Arguments
  * `input::IO`: Input stream containing BED data.

# Returns
An iterator over `BEDInterval`s contained in the file.
"""
function read(input::IO, ::Type{BED})
    return BEDIterator(BEDParserImpl.BEDParser(input), Nullable{BEDInterval}())
end


function start(it::BEDIterator)
    advance!(it)
    return nothing
end


function advance!(it::BEDIterator)
    isdone = !BEDParserImpl.advance!(it.parser)
    if isdone
        it.nextitem = Nullable{BEDInterval}()
    else
        it.nextitem = BEDParserImpl.takevalue!(it.parser)
    end
end


function next(it::BEDIterator, state::Nothing)
    item = get(it.nextitem)
    advance!(it)
    return item, nothing
end


function done(it::BEDIterator, state::Nothing)
    return isnull(it.nextitem)
end
