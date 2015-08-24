module TestSeq

using FactCheck
include("../../src/helper/pir.jl")

facts("PIR") do
    context("Conversion") do
        context("convert PIR to FASTA-Format") do
            pir = ">P1;5fd1\n  structureX:5fd1:1    :A:106  :A:ferredoxin:Azotobacter vinelandii: 1.90: 0.19\nAFVVTDNCIKCKYTDCVEVCPVDCFYEGPNFLVIHPDECIDCALCEPECPAQAIFSEDEVPEDMQEFIQLNAELA\nEVWPNITEKKDPLPDAEDWDGVKGKLQHLER*"
            expected = ">P1;5fd1;PIR=(structureX:5fd1:1    :A:106  :A:ferredoxin:Azotobacter vinelandii: 1.90: 0.19)\nAFVVTDNCIKCKYTDCVEVCPVDCFYEGPNFLVIHPDECIDCALCEPECPAQAIFSEDEVPEDMQEFIQLNAELA\nEVWPNITEKKDPLPDAEDWDGVKGKLQHLER*"
            @fact convert_PIR_to_FASTA(pir) --> expected
        end
    end
end

end
