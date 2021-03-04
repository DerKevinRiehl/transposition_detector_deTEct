############################################################################
##### Transposition Event Detector deTEct - part of Transposon Ultimate ####
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

# Imports
import sys
import os
from ParserPBSV import parsePBSV_SVs
#from ParserSVIM import parseSVIM_SVs
from ParserSniffles import parseSniffles_SVs
from TranspositionDetectionPipeline import detectTranspositionEvents
from TranspositionDetector_Help import helpExplanations

# Methods
def getArgument(args, title):
    for i in range(0, len(args)):
        if(args[i].startswith("-"+title)):
            if(i<len(args)-1):
                return args[i+1]
            else:
                return ""
    return ""

## Main Code
args        = sys.argv[1:]   
if("-help" in args or "-h" in args):
    helpExplanations()
else:
    seqHeadFile = getArgument(args, "seqHeadTXT")
    transpAnnot = getArgument(args, "transpGFF3")
    assemblyFas = getArgument(args, "assmFasta")
    svTool      = getArgument(args, "svTool")
    svFile      = getArgument(args, "svFile")
#    bedFile     = getArgument(args, "bedFile")
    ouFile1     = getArgument(args, "outParsedFile")
    ouFile2     = getArgument(args, "outResultFile")
    
    # Check inputs
    if(seqHeadFile==""):
        print("ERROR: please provide valid sequence heads file (seqHeadFile)!")
        sys.exit(-1)
    else:
        if(not os.path.isfile(seqHeadFile)):
            print("ERROR: provided seqHeadFile \""+seqHeadFile+"\" could not be found!")
            sys.exit(-1)
    if(transpAnnot==""):
        print("ERROR: please provide valid transposon annotation (transpGFF3)!")
        sys.exit(-1)
    else:
        if(not os.path.isfile(transpAnnot)):
            print("ERROR: provided transpGFF3 \""+transpAnnot+"\" could not be found!")
            sys.exit(-1)
    if(assemblyFas==""):
        print("ERROR: please provide valid reference genome assembly (assmFasta)!")
        sys.exit(-1)
    else:
        if(not os.path.isfile(assemblyFas)):
            print("ERROR: provided assmFasta \""+assemblyFas+"\" could not be found!")
            sys.exit(-1)
    if(svFile==""):
        print("ERROR: please provide valid structural variant file (svFile)!")
        sys.exit(-1)
    else:
        if(not os.path.isfile(svFile)):
            print("ERROR: provided svFile \""+svFile+"\" could not be found!")
            sys.exit(-1)
    if(ouFile1==""):
        print("ERROR: please provide desired outputfile for VCF parser (ouFile1)!")
        sys.exit(-1)
    if(ouFile2==""):
        print("ERROR: please provide desired outputfile for transposition events (ouFile2)!")
        sys.exit(-1)
            
    # Parse Structural Variant VCF File
    if(svTool=="pbsv"):
        parsePBSV_SVs(seqHeadFile, svFile, ouFile1)
#    if(svTool=="svim"):
#        parseSVIM_SVs(seqHeadFile, svFile, bedFile, ouFile1)
    elif(svTool=="sniffles"):
        parseSniffles_SVs(seqHeadFile, svFile, ouFile1)
    else:
        print("ERROR: no valid svTool \""+svTool+"\" could be found... (valid options \"pbsv\" and \"sniffles\")")
        sys.exit(-1)
        
    # Run Pipeline
    detectTranspositionEvents(seqHeadFile, transpAnnot, assemblyFas, ouFile1, ouFile2)