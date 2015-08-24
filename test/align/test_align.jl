module TestAlign

using FactCheck
using Bio
using Bio.Align

facts("Alignments") do
    context("Global") do
        context("Pairwise alignment") do
            context("NeedlemanWunsch") do
                nw = NeedlemanWunsch(UnitScore)
                @fact align(nw, "GCATGCU", "GATTACA") --> ("GCA-TGCU", "G-ATTACA")

                expected = [
                    ("AABB","-A-B"),
                    ("AABB","A--B"),
                    ("AABB","-AB-"),
                    ("AABB","A-B-")
                ]
                @fact alignAll("AABB", "AB") --> expected
            end
        end
        context("Multiple") do
            context("Center-Star-Method") do
                sequences = [
                    "CCTGCTGCAG",
                    "GATGTGCCG",
                    "GATGTGCAG",
                    "CCGCTAGCAG",
                    "CCTGTAGG"
                ]

                msa = [
                    "CCTGCT-GCAG",
                    "GATG-T-GCCG",
                    "GATG-T-GCAG",
                    "CC-GCTAGCAG",
                    "CCTG-T--AGG",
                ]
                expected = Bio.Align.MultipleSequenceAlignment(msa)
                result = align(sequences)
                @fact result.alignment --> expected.alignment
            end
        end
    end
end

end # TestAlign
