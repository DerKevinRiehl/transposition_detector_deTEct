############################################################################
##### Transposition Event Detector deTEct - part of Transposon Ultimate ####
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

# Methods
def parseSVIM_SVs(sequenceHeadsFile, inputVcfFile, inputBreakendsBEDFile, outputGffFile):
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
        quality = parts[5]
        if("svim.INV" in line):
            end = str(int(parts[7].split(";END=")[1].split(";")[0])+1)
            fW.write(sequenceDictB[chrom]+"\t"+"SVIM"+"\t"+"inversion"+"\t"+start+"\t"+end+"\t"+quality+"\t"+"+"+"\t"+"."+"\t"+parts[2])
        if("svim.DEL" in line):
            end = str(int(parts[7].split(";END=")[1].split(";")[0])+1)
            fW.write(sequenceDictB[chrom]+"\t"+"SVIM"+"\t"+"deletion"+"\t"+start+"\t"+end+"\t"+quality+"\t"+"+"+"\t"+"."+"\t"+parts[2])
        if("svim.INS" in line):
            end = str(int(start)+int(parts[7].split(";SVLEN=")[1].split(";")[0])+1)
            fW.write(sequenceDictB[chrom]+"\t"+"SVIM"+"\t"+"insertion"+"\t"+start+"\t"+end+"\t"+quality+"\t"+"+"+"\t"+"."+"\t"+parts[2])
        if("svim.BND" in line):
            line = fR.readline()
            continue
        if("svim.DUP_TANDEM" in line):
            end = str(int(parts[7].split(";END=")[1].split(";")[0])+1)
            quality = ""
            fW.write(sequenceDictB[chrom]+"\t"+"SVIM"+"\t"+"duplication_tandem"+"\t"+start+"\t"+end+"\t"+quality+"\t"+"+"+"\t"+"."+"\t"+parts[2])
        if("svim.DUP_INT" in line):
            end = str(int(parts[7].split(";END=")[1].split(";")[0])+1)
            fW.write(sequenceDictB[chrom]+"\t"+"SVIM"+"\t"+"duplication_interspersed"+"\t"+start+"\t"+end+"\t"+quality+"\t"+"+"+"\t"+"."+"\t"+parts[2])
        fW.write("\n")
        line = fR.readline()
    fR.close()
    fR = open(inputBreakendsBEDFile, "r")
    line = fR.readline()
    while line!="":
        parts = line.replace("\n","").split("\t")
        chrom = parts[0]
        start = str(int(parts[1])+1)
        end   = str(int(parts[2])+1)
        info  = parts[3]+parts[5]
        for key in sequenceDictB:
            info = info.replace(key,sequenceDictB[key])
        fW.write(sequenceDictB[chrom]+"\t"+"SVIM"+"\t"+"translocation"+"\t"+start+"\t"+end+"\t"+"."+"\t"+"+"+"\t"+"."+"\t"+info)
        fW.write("\n")
        line = fR.readline()
    fR.close()
    fW.close()

#import os
#folder = "C:/Users/kevin/Desktop/detectionTest"
#sequenceHeadsFile     = os.path.join(folder,"PRJEB28388Annotations","sequence_heads.txt")
#inputVcfFile          = os.path.join(folder,"RICCIO_EC_amares_SVIM_MINIMAP2.vcf")
#inputBreakendsBEDFile = os.path.join(folder,"candidates_SVIM_MINIMAP2","candidates_breakends.bed")
#outputGffFile         = os.path.join(folder,"RICCIO_EC_amares_SVIM_MINIMAP2.vcf.gff3")
#parseSVIM_SVs(sequenceHeadsFile, inputVcfFile, inputBreakendsBEDFile, outputGffFile)