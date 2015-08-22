#=
convertes .PIR files created with modeller to FASTA files

idea taken from https://www.biostars.org/p/131001/

For addition info see "PIR parsing in Biopython"
http://biopython.org/DIST/docs/api/Bio.SeqIO.PirIO-pysrc.html
=#

"""
reads a given PIR inputfile, converts its content to FASTA-format
and writes it to the specified outputfile
"""
function read_PIR_and_save_as_FASTA(inputfile::String, outputfile::String)
    content = open(readall, inputfile)
    fasta = convert_PIR_to_FASTA(content)
    open(outputfile, "w+") do fh
        write(fh, fasta)
    end
end

"""
regex that deletes the line containing additional information in a PIR-file
and appends the deleted information to the FASTA-HEADER surrounded by `PIR=()`
"""
function convert_PIR_to_FASTA(pir)
    return replace(pir, r"(>.*)\n[ ]*(.*)(\n)", s"\1;PIR=(\2)\3")
end
