############################################################################
##### Transposition Event Detector deTEct - part of Transposon Ultimate ####
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

# Methods
def parseSniffles_SVs(sequenceHeadsFile, inputVcfFile, outputGffFile):
    # Load Sequence Heads
    sequenceDictA = {}
    sequenceDictB = {}
    f = open(sequenceHeadsFile, "r")
    line = f.readline()
    while line!="":
        parts = line.replace("\n","").replace(">","").split("\t")
        sequenceDictA[parts[0]] = parts[1]
        sequenceDictB[parts[1]] = parts[0]
        line = f.readline()
    f.close()
    # Convert Structural Variants to Annotations
    fW = open(outputGffFile,"w+")
    fR = open(inputVcfFile, "r")
    line = fR.readline()
    while line!="":
        if(line.startswith("#")):
            line = fR.readline()
            continue
        parts = line.replace("\n","").split("\t")
        chrom = parts[0]
        start = str(int(parts[1])+1)
        end = str(int(parts[7].split(";END=")[1].split(";")[0])+1)
        info = parts[7]
        for key in sequenceDictB:
            info = info.replace(key,sequenceDictB[key])
            info = info.replace(key.upper(),sequenceDictB[key].upper())
        if("SVTYPE=DEL" in line):
            fW.write(sequenceDictB[chrom]+"\t"+"SVIM"+"\t"+"deletion"+"\t"+start+"\t"+end+"\t"+"."+"\t"+"+"+"\t"+"."+"\t"+info)
        if("SVTYPE=DUP" in line):
            fW.write(sequenceDictB[chrom]+"\t"+"SVIM"+"\t"+"duplication"+"\t"+start+"\t"+end+"\t"+"."+"\t"+"+"+"\t"+"."+"\t"+info)
        if("SVTYPE=INV" in line):
            fW.write(sequenceDictB[chrom]+"\t"+"SVIM"+"\t"+"inversion"+"\t"+start+"\t"+end+"\t"+"."+"\t"+"+"+"\t"+"."+"\t"+info)
        if("SVTYPE=INVDUP" in line):
            fW.write(sequenceDictB[chrom]+"\t"+"SVIM"+"\t"+"inversion_duplicate_unknownBound"+"\t"+start+"\t"+end+"\t"+"."+"\t"+"+"+"\t"+"."+"\t"+info)
        if("SVTYPE=TRA" in line):
            fW.write(sequenceDictB[chrom]+"\t"+"SVIM"+"\t"+"translocation"+"\t"+start+"\t"+end+"\t"+"."+"\t"+"+"+"\t"+"."+"\t"+info)
        if("SVTYPE=INS" in line):
            fW.write(sequenceDictB[chrom]+"\t"+"SVIM"+"\t"+"insertion"+"\t"+start+"\t"+end+"\t"+"."+"\t"+"+"+"\t"+"."+"\t"+info)
        fW.write("\n")
        line = fR.readline()
    fR.close()
    fW.close()
  
#import os
#folder = "C:/Users/kevin/Desktop/detectionTest"
#sequenceHeadsFile     = os.path.join(folder,"PRJEB28388Annotations","sequence_heads.txt")
#inputVcfFile          = os.path.join(folder,"RICCIO_EC_amares_Sniffles_MINIMAP2.vcf")
#inputBreakendsBEDFile = os.path.join(folder,"candidates_NGLMR","candidates_breakends.bed")
#outputGffFile         = os.path.join(folder,"RICCIO_EC_amares_Sniffles_MINIMAP2.vcf.gff3")
#parseSniffles_SVs(sequenceHeadsFile, inputVcfFile, outputGffFile)
#  