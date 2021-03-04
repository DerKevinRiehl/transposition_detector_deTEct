############################################################################
##### Transposition Event Detector deTEct - part of Transposon Ultimate ####
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

# Methods
def helpExplanations():
    print("Transposition Event Detector deTEct (v1.0, 2021, Kevin Riehl)")
    print("")
    print("USAGE")
    print("  transposition_deTEct [-help] -seqHeadTXT <FILE> -transpGFF3 <FILE> -assmFasta <FILE> -svTool <FILE> -svFile <FILE> -outParsedFile <FILE> -outResultFile <FILE>")
    print("")
    print("(mandatory):")
    print("  -seqHeadTXT <String>")
    print("    Sequence head names, TXT file (produced by reasonaTE)")
    print("  -transpGFF3 <String>")
    print("    Transposon annotation file, GFF3 file (produced by reasonaTE)")
    print("  -assmFasta <String>")
    print("    Assembly file of reference genome, FASTA file")
    print("  -svTool <String>")
    print("    Structural variant detection tool: \"pbsv\", \"sniffles\"")
    print("  -svFile <String>")
    print("    Structural variant detection output file, VCF file")
    print("  -outParsedFile <String>")
    print("    Target file for VCF parsed outputs")
    print("  -outResultFile <String>")
    print("    Target file for final results with transposition events")
    print("")