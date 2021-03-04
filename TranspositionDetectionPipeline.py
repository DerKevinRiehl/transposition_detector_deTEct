############################################################################
##### Transposition Event Detector deTEct - part of Transposon Ultimate ####
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

# Imports
import ast 
from operator import itemgetter

# Methods
def getChromosomeLengths(genomeFile):
    lengths = {}
    f = open(genomeFile,"r")
    line = f.readline().replace("\n","")
    seq = 0
    tot = 0
    chrom = ""
    while line!="":
        if(line.startswith(">")):
            if(chrom!="" and chrom!="M"):
                lengths[chrom] = seq
                tot += seq
            chrom = line.replace(">","").replace("\n","")
            seq = 0
        else:        
            seq += len(line.replace("\n",""))
        line = f.readline().replace("\n","")
    f.close()
    lengths[chrom] = seq
    tot += seq
    return lengths,tot

def isIntersection(a,b,c,d, tol):
    if(a-tol<=c):
        if(b+tol>=c):
            return 1
    else:
        if(a-tol<=d):
            return 1
    return 0

def getOverlap(a,b,c,d):
    if(a<c):
        return min(d-c,b-c)
    else:
        return min(b-a,d-a)
    
def detectTranspositionEvents(seqHeadsTXT, transpGFF3, assemblyFASTA, svGFF3, outFileGFF3):
    # Load Sequence Heads
    sequenceDictA = {}
    sequenceDictB = {}
    f = open(seqHeadsTXT, "r")
    line = f.readline()
    while line!="":
        parts = line.replace("\n","").replace(">","").split("\t")
        sequenceDictA[parts[0]] = "chr"+parts[1]
        sequenceDictB["chr"+parts[1]] = parts[0]
        line = f.readline()
    f.close()
    # Load Lengths
    chromLenghts, totalLength = getChromosomeLengths(assemblyFASTA)
    # Load Transposon Annotations
    annotations = {}
    f = open(transpGFF3, "r")
    line = f.readline()
    while line!="":
        parts = line.replace("\n","").split("\t")
        chrom = parts[0]
        if(chrom not in annotations):
            annotations[chrom] = list()
        start = int(parts[3])
        end   = int(parts[4])
        info  = parts[8].split("transposon=")[1].split(";")[0]
        clasL = parts[8].split("class=")[1].split(";")[0]
        annotations[chrom].append([start,end,info,clasL])
        line = f.readline()
    f.close()
    for chrom in annotations:
        annotations[chrom] = sorted(annotations[chrom], key=itemgetter(0))
    # Filter1: Filter Structural Variants by min and max length
    filterLenMIN = 50
    filterLenMAX = totalLength*0.01
    fW = open(svGFF3+".filtered.gff3","w+")
    fR = open(svGFF3, "r")
    line = fR.readline()
    while line!="":
        line = line.replace("\n","")
        if(line==""):
            continue
        parts = line.replace("\n","").split("\t")
        start = int(parts[3])
        end   = int(parts[4])
        if(abs(end-start)>filterLenMIN and abs(end-start)<filterLenMAX):
            fW.write(line)
            fW.write("\n")
        line = fR.readline()
    fR.close()
    fW.close()
    # Filter2: Remove Duplicates in SVs
    svAnnotations = {}
    fR = open(svGFF3+".filtered.gff3", "r")
    line = fR.readline()
    while line!="":
        line = line.replace("\n","")
        if(line==""):
            continue
        parts = line.replace("\n","").split("\t")
        chrom = parts[0]
        start = int(parts[3])
        end   = int(parts[4])
        typ   = parts[2]
        if(chrom not in svAnnotations):
            svAnnotations[chrom] = list()
        svAnnotations[chrom].append([start,end,typ,line])
        line = fR.readline()
    fR.close()
    for chrom in svAnnotations:
        svAnnotations[chrom] = sorted(svAnnotations[chrom], key=itemgetter(0))
    svAnnotations2 = {}
    for chrom in svAnnotations:
        iteration = 0
        run = True
        while run:
            iteration+=1
            svAnnotations2[chrom] = list()
            duplFilter = list()
            for ix1 in range(0,len(svAnnotations[chrom])):
                if(ix1%1000==0):
                    print("\t","...","\t",chrom,"\t",ix1,"\t",len(svAnnotations[chrom]))
                ann1 = svAnnotations[chrom][ix1]
                if(ann1[3] in duplFilter):
                    continue
                leng  = abs(ann1[1]-ann1[0])+1
                foundIx2 = -1
                for ix2 in range(ix1+1,len(svAnnotations[chrom])):
                    ann2 = svAnnotations[chrom][ix2]
                    if(ann2[3] in duplFilter):
                        continue
                    if(ann2[0]>ann1[1]):
                        break
                    if(ann1[2]==ann2[2]):
                        if(isIntersection(ann1[0],ann1[1],ann2[0],ann2[1],0)):
                            duplFilter.append(ann2[3])
                            foundIx2 = ix2
                            break
                if(foundIx2==-1):
                    svAnnotations2[chrom].append(svAnnotations[chrom][ix1])
                else:
                    svAnnotations[chrom][ix1][1] = svAnnotations[chrom][foundIx2][1]
                    line = svAnnotations[chrom][ix1][3]
                    parts = line.split("\t")
                    parts[4] = str(svAnnotations[chrom][ix1][1])
                    parts[8] = "merged;"+parts[8]
                    svAnnotations[chrom][ix1][3] = "\t".join(parts)
                    svAnnotations2[chrom].append(svAnnotations[chrom][ix1])
            svAnnotations[chrom] = svAnnotations2[chrom]
            print(chrom,"\t","Iteration ",iteration,"\tdropped\t",len(duplFilter))
            if(len(duplFilter)==0):
                run = False
    fW = open(svGFF3+".filtered2.gff3","w+")
    chromosomes = list(svAnnotations.keys())
    chromosomes.sort()
    for chrom in chromosomes:
        for ann in svAnnotations[chrom]:
            fW.write(ann[3])
            fW.write("\n")
    fW.close()
    # Matching: Match Transposons to Structural Variants
    tolerance = 0.1
    fW = open(svGFF3+".matches.gff3","w+")
    fR = open(svGFF3+".filtered2.gff3", "r")
    line = fR.readline()
    while line!="":
        line = line.replace("\n","")
        if(line==""):
            continue
        parts = line.replace("\n","").split("\t")
        chrom = parts[0] 
        start = int(parts[3])
        end   = int(parts[4])
        leng  = abs(end-start)+1
        info  = parts[8]
        candidates = list()
        for ann in annotations[chrom]:
            if(isIntersection(start,end,ann[0],ann[1], leng*tolerance)):
                candidates.append(ann[2])
        if(len(candidates)!=0):
            parts[8] = str(candidates)+";"+info
            fW.write("\t".join(parts))
            fW.write("\n")
        line = fR.readline()
    fR.close()
    fW.close()
    # Filter3: Filter Matches by similar length to Transposons and Structural Variants 
    lenTol = 0.5
    fW = open(outFileGFF3,"w+")
    fR = open(svGFF3+".matches.gff3", "r")
    line = fR.readline()
    while line!="":
        line = line.replace("\n","")
        if(line==""):
            continue
        parts = line.replace("\n","").split("\t")
        chrom = parts[0] 
        start = int(parts[3])
        end   = int(parts[4])
        candidates = ast.literal_eval(parts[8].split(";")[0]) 
        newCandidates = list()
        leng1  = abs(end-start)+1
        candC = -1
        for c in candidates:
            cAnn = -1
            for ann in annotations[chrom]:
                if(ann[2]==c):
                    cAnn = ann
                    candC = ann
                    break
            leng2 = abs(ann[1]-ann[0])+1
            lengDiff = abs(leng2-leng1)
            if(lengDiff<lenTol*leng1 and lengDiff<lenTol*leng2):
                newCandidates.append(c)
        if(len(newCandidates)==1):
            parts[8] = "Transposon="+str(newCandidates[0].replace(" ",""))+";Class="+candC[3]+";"+";".join(parts[8].split(";")[1:])
            fW.write("\t".join(parts))
            fW.write("\n")
        elif(len(newCandidates)>1):
            maxOverlap = -1
            maxId = 0
            for c in newCandidates:
                cAnn = -1
                for ann in annotations[chrom]:
                    if(ann[2] in c):
                        cAnn = ann
                        break
                overlap = getOverlap(start,end,cAnn[0],cAnn[1])
                if(overlap>maxOverlap):
                    maxOverlap = overlap
                    maxId = cAnn[2]
            if(maxId!=0):
                parts[8] = "Transposon="+str(maxId.replace(" ",""))+";Class="+cAnn[3]+";"+";".join(parts[8].split(";")[1:])
                fW.write("\t".join(parts))
                fW.write("\n")
            else:
                parts[8] = "Transposon="+str(newCandidates[0].replace(" ",""))+";Class="+cAnn[3]+";multiple_transposons;"+";".join(parts[8].split(";")[1:])
                fW.write("\t".join(parts))
                fW.write("\n")
        line = fR.readline()
    fR.close()
    fW.close()

#folder         = "C:/Users/kevin/Desktop/detectionTest"
#variantFile   = "RICCIO_EC_amares_SVIM_MINIMAP2.vcf" # "RICCIO_EC_amares_SVIM.vcf" # "RICCIO_EC_amares_Sniffles.vcf"
#seqHeadsTXT   = os.path.join(folder,"PRJEB28388Annotations","sequence_heads.txt")
#transpGFF3    = os.path.join(folder,"PRJEB28388Annotations","FinalAnnotations_Transposons.gff3")
#assemblyFASTA = os.path.join(folder,"PRJEB28388Annotations","sequence.fasta")
#svGFF3        = os.path.join(folder,variantFile+".gff3")
#outFileGFF3   = os.path.join(folder,variantFile+".transpositionEvents.gff3")
#detectTranspositionEvents(seqHeadsTXT, transpGFF3, assemblyFASTA, svGFF3, outFileGFF3)
    
#folder         = "C:/Users/kevin/Desktop/detectionTest"
#variantFile   = "RICCIO_EC_amares_PBSV_PBMM2.var.vcf" # "RICCIO_EC_amares_SVIM.vcf" # "RICCIO_EC_amares_Sniffles.vcf"
#seqHeadsTXT   = os.path.join(folder,"PRJEB28388Annotations","sequence_heads.txt")
#transpGFF3    = os.path.join(folder,"PRJEB28388Annotations","FinalAnnotations_Transposons.gff3")
#assemblyFASTA = os.path.join(folder,"PRJEB28388Annotations","sequence.fasta")
#svGFF3        = os.path.join(folder,variantFile+".gff3")
#outFileGFF3   = os.path.join(folder,variantFile+".transpositionEvents.gff3")
#detectTranspositionEvents(seqHeadsTXT, transpGFF3, assemblyFASTA, svGFF3, outFileGFF3)