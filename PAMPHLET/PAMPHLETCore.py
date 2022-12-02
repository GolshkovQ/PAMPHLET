#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
PAMPHLET - PAM Prediction HomoLogous Enhancement Toolkit
Author: Chen Qi, Baitao Li, Lei Huang
University of Chinese Academy of Sciences, College of Life Sciences, Beijing, China
BGI Research, Shenzhen, China
Email: qichen@genomics.cn; libaitao@genomics.cn; huanglei@genomics.cn
'''

import sys
running_python3 = False
if sys.version_info > (3, 0):
    running_python3 = True

import os
import time


print('''

██████   █████  ███    ███ ██████  ██   ██ ██      ███████ ████████ 
██   ██ ██   ██ ████  ████ ██   ██ ██   ██ ██      ██         ██    
██████  ███████ ██ ████ ██ ██████  ███████ ██      █████      ██    
██      ██   ██ ██  ██  ██ ██      ██   ██ ██      ██         ██    
██      ██   ██ ██      ██ ██      ██   ██ ███████ ███████    ██   

''')

print("PAMPHLET")
print("PAM Prediction HomoLogous Enhancement Toolkit")
print("Author: Chen Qi, Baitao Li, Lei Huang")
print("University of Chinese Academy of Sciences, College of Life Sciences, Beijing, China")
print("BGI Research, Shenzhen, China")

print("===============================================================")
print("Start loading resources...")
resourcesTime = time.time()
runningTime = time.time()

from PAMPHLET import PAMPHLETParams
from PAMPHLET import PAMPHLETReviser
from PAMPHLET import PAMPHLETResources
from PAMPHLET import OnlineResources

print("Loading resources finished. Time used: %s seconds" % (time.time() - resourcesTime))

def main():
    print("===============================================================")
    print("Initializing parameters...")

    arg_parser = PAMPHLETParams.get_arg_parser()
    opts = arg_parser.parse_args()

    SpacerFile = opts.spacer
    RepeatSeq = opts.repeat
    ProteinFile = opts.protein
    OutDir = opts.outdir
    Oriention = opts.orientation
    UniqueMode = opts.unique
    ProteinBlastResult = opts.proteinblast
    ProteinLength = opts.proteinlen
    SpacerBlastDB = opts.spacerdb
    FlankLength = opts.flanklen
    RefMode = opts.refmode
    FreqMode = opts.freqmode

    ### FILE CHECKS

    PAMPHLETParams.check_protein_argument(RepeatSeq,ProteinFile,ProteinBlastResult,ProteinLength)
    PAMPHLETParams.check_spacer_sequences(SpacerFile)
    PAMPHLETParams.check_environment()
    PAMPHLETParams.check_outdir(OutDir)

    ### MODE CHECKS AND SPACER REVISE

    print("Initializing parameters finished. Time used: %s seconds" % (time.time() - runningTime))
    print("===============================================================")

    if ProteinFile != None and ProteinBlastResult == None:
        print("Using protein sequences to revise spacers...")
        PAMPHLETParams.check_protein_files(ProteinFile)
        finalspacer = PAMPHLETReviser.main(SpacerFile,RepeatSeq,ProteinFile,OutDir,UniqueMode,ProteinBlastResult,ProteinLength)
        print("Using protein sequences to revise spacers finished. Time used: %s seconds" % (time.time() - runningTime))
        print("===============================================================")
        if finalspacer == False:
            finalspacer = SpacerFile
    elif ProteinFile == None and ProteinBlastResult != None:
        print("Using protein blast results to revise spacers...")
        finalspacer = PAMPHLETReviser.main(SpacerFile,RepeatSeq,ProteinBlastResult,OutDir,UniqueMode,ProteinBlastResult,ProteinLength,RefMode)
        if finalspacer == False:
            finalspacer = SpacerFile
    else:
        print("No protein sequences provided, skip spacer revise step.")
        print("===============================================================")
        finalspacer = SpacerFile

    print("Start running spacer blast...")
    print("This step may take a long time, please wait patiently...")
    print("===============================================================")

    ### RENAME SPACER ID
    SpacerIDLenList = PAMPHLETResources.rename_spacer_id(finalspacer,finalspacer+".label",Oriention)
    finalspacer = finalspacer + ".tmp"

    ### SETTING SPACER BLAST DIRECTORY
    SpacerBlastDir = os.path.join(OutDir,"SpacerBlast/")
    PAMPHLETParams.check_directory(SpacerBlastDir)
    SpacerInputDir = os.path.join(SpacerBlastDir,"Input/")
    PAMPHLETParams.check_directory(SpacerInputDir)
    SpacerBlastOutputDir = os.path.join(SpacerBlastDir,"RawOutput/")
    PAMPHLETParams.check_directory(SpacerBlastOutputDir)

    ### CHECK SPACER LENGTH, IF LENGTH > 20, CUT TO 20
    SpacerNumbers = os.popen("grep -c '>' %s" % finalspacer).read().strip()
    if int(SpacerNumbers) > 20:
        PAMPHLETResources.split_spacer_file(finalspacer,SpacerInputDir,int(SpacerNumbers))
    else:
        os.system("cp %s %s" % (finalspacer,SpacerInputDir))

    ### BLAST SPACER AGAINST REPEAT DATABASE
    for InputSpacerFile in os.listdir(SpacerInputDir):
        InputSpacerFile = os.path.join(SpacerInputDir,InputSpacerFile)
        OutputSpacerFile = os.path.join(SpacerBlastOutputDir,os.path.basename(InputSpacerFile))
        OnlineResources.run_spacer_blast(InputSpacerFile,OutputSpacerFile,SpacerBlastDB)

    ### MERGE BLAST OUTPUT IF THERE HAVE MULTIPLE BLAST OUTPUT
    MergedSpacerBlastOutput = os.path.join(SpacerBlastDir,"MergedSpacerBlastOutput.txt")
    if len(os.listdir(SpacerBlastOutputDir)) > 1:
        PAMPHLETResources.merge_blast_output(SpacerBlastOutputDir,MergedSpacerBlastOutput)
    else:
        os.system("cp %s %s" % (os.path.join(SpacerBlastOutputDir,os.listdir(SpacerBlastOutputDir)[0]),MergedSpacerBlastOutput))

    ### GET SIGNIFICANT SPACER BLAST OUTPUT AND NON-SIGNIFICANT HIT SPACER ID
    SignificantSpacerBlastOutput = os.path.join(SpacerBlastDir,"SignificantSpacerBlastOutputRaw.txt")
    InSignificantSpacerID = PAMPHLETResources.get_significant_spacer_blast_output(MergedSpacerBlastOutput,SignificantSpacerBlastOutput,SpacerIDLenList)
    
    ### SET FINAL SPACER HIT FILE
    FinalSpacerHitFile = os.path.join(OutDir,"FinalSpacerHit.txt")

    ### CHECK INSIGNIFICANT ID NUMBERS, IF IS 0, CONTINUE TO SEQDUMP DOWNLOAD; ELSE, REVISE SPACER LENGTH.
    if len(InSignificantSpacerID) != 0:

        print("First round spacer blast finished. Time used: %s seconds" % (time.time() - runningTime))
        print("Start revising spacer length for non-significant spacer blast hits...")
        print("This step may take a long time, please wait patiently...")
        print("Revised spacers: "+"; ".join(InSignificantSpacerID))
        print("===============================================================")

        ## REVISE SPACER BLAST DIRECTORY
        RevisedLengthSpacerFile = os.path.join(SpacerBlastDir,"RevisedLengthSpacerFile.txt")
        PAMPHLETResources.revise_spacer_length(finalspacer,RevisedLengthSpacerFile,InSignificantSpacerID)
        RevisedSpacerBlastDir = os.path.join(OutDir,"RevisedSpacerBlast/")
        PAMPHLETParams.check_directory(RevisedSpacerBlastDir)
        RevisedSpacerInputDir = os.path.join(RevisedSpacerBlastDir,"Input/")
        PAMPHLETParams.check_directory(RevisedSpacerInputDir)
        RevisedSpacerBlastOutputDir = os.path.join(RevisedSpacerBlastDir,"RawOutput/")
        PAMPHLETParams.check_directory(RevisedSpacerBlastOutputDir)

        ### CHECK SPACER LENGTH, IF LENGTH > 20, CUT TO 20
        if len(InSignificantSpacerID) >= 20:
            PAMPHLETResources.split_spacer_file(RevisedLengthSpacerFile,RevisedSpacerInputDir,len(InSignificantSpacerID))
        else:
            os.system("cp %s %s" % (RevisedLengthSpacerFile,RevisedSpacerInputDir))
        
        ### BLAST SPACER AGAINST REPEAT DATABASE
        for InputSpacerFile in os.listdir(RevisedSpacerInputDir):
            InputSpacerFile = os.path.join(RevisedSpacerInputDir,InputSpacerFile)
            OutputSpacerFile = os.path.join(RevisedSpacerBlastOutputDir,os.path.basename(InputSpacerFile))
            OnlineResources.run_spacer_blast(InputSpacerFile,OutputSpacerFile,SpacerBlastDB)
        
        ### MERGE BLAST OUTPUT IF THERE HAVE MULTIPLE BLAST OUTPUT
        if len(os.listdir(RevisedSpacerBlastOutputDir)) > 1:
            print("Merging blast output...")
            RevisedSpacerBlastOutput = os.path.join(RevisedSpacerBlastDir,"MergedSpacerBlastOutput.txt")
            PAMPHLETResources.merge_blast_output(RevisedSpacerBlastOutputDir,RevisedSpacerBlastOutput)
        else:
            RevisedSpacerBlastOutput = os.path.join(RevisedSpacerBlastDir,"MergedSpacerBlastOutput.txt")
            os.system("cp %s %s" % (os.path.join(RevisedSpacerBlastOutputDir,os.listdir(RevisedSpacerBlastOutputDir)[0]),RevisedSpacerBlastOutput))

        ### GET SIGNIFICANT SPACER BLAST OUTPUT AND NON-SIGNIFICANT HIT SPACER ID
        print("Getting significant spacer blast output...")
        RevisedSpacerLen = PAMPHLETResources.build_revised_spacer_len_dict(SpacerIDLenList,InSignificantSpacerID)
        SignificantRevisedSpacerBlastOutput = os.path.join(RevisedSpacerBlastDir,"SignificantRevisedSpacerBlastOutputRaw.txt")
        FalseSpacerIDList = PAMPHLETResources.get_significant_spacer_blast_output(RevisedSpacerBlastOutput,SignificantRevisedSpacerBlastOutput,RevisedSpacerLen)

        print("Revised false spacer: "+"; ".join(FalseSpacerIDList))
        print("Getting final spacer blast output...")

        ### CHECK SIGNIFICANT SPACER BLAST OUTPUT, IF IS 0, CONTINUE TO SEQDUMP DOWNLOAD; ELSE, REVISE SPACER LENGTH.
        if os.path.getsize(SignificantRevisedSpacerBlastOutput) == 0:
            print("WARNING: No significant spacer blast hit for revised spacer length.")
            if os.path.getsize(SignificantSpacerBlastOutput) == 0:
                print("No significant spacer blast hit for both original spacer and revised spacer .")
                print("Cannot find spacer PAM sequence.")
                sys.exit(1)
            else:
                print("Use original spacer blast hit.")
                os.system("mv %s %s" % (SignificantSpacerBlastOutput,FinalSpacerHitFile))
        else:
            PAMPHLETResources.merge_spacer_blast_output(SignificantSpacerBlastOutput,SignificantRevisedSpacerBlastOutput,FinalSpacerHitFile)
    
    else:
        print("Use original spacer blast hit.")
        print("Spacer blast finished. Time used: %s seconds" % (time.time() - runningTime))
        print("===============================================================")
        os.system("mv %s %s" % (SignificantSpacerBlastOutput,FinalSpacerHitFile))

    ### GET GENOME SEQDUMP FILE
    print("Start downloading genome seqdump file...")
    GenomeSeqdumpFile = os.path.join(OutDir,"GenomeSeqdump.txt")
    OnlineResources.get_genome_seqdump_files(FinalSpacerHitFile,GenomeSeqdumpFile)
    
    print("Genome seqdump file download finished. Time used: %s seconds" % (time.time() - runningTime))
    print("===============================================================")

    ### GET FLANK SEQUENCE BASED ON SIGNIFICANT SPACER BLAST OUTPUT

    print("Extract flank sequence and calculate PAM...")

    UpStreamFlankSeq = os.path.join(OutDir,"UpStreamFlankSeq.txt")
    DownStreamFlankSeq = os.path.join(OutDir,"DownStreamFlankSeq.txt")
    tempMinCEDDir = os.path.join(OutDir,"SpacerSeqdumpMinCED/")
    PAMPHLETParams.check_directory(tempMinCEDDir)
    UpstreamDict, DownstreamDict = PAMPHLETResources.get_flank_seq(FinalSpacerHitFile,GenomeSeqdumpFile,UpStreamFlankSeq,DownStreamFlankSeq,FlankLength,tempMinCEDDir)

    if len(UpstreamDict) == 0 or len(DownstreamDict) == 0:
        print("Cannot find upstream or downstream flank sequence.")
        sys.exit(1)

    ### CONVERT SEQUENCE DICT TO FREQUENCY LIST
    UpstreamFreqList = PAMPHLETResources.convert_seq_to_freq(UpstreamDict,FlankLength,FreqMode)
    DownstreamFreqList = PAMPHLETResources.convert_seq_to_freq(DownstreamDict,FlankLength,FreqMode)

    ### WRITE FREQDICT TO FILE, BUILD MATRIX
    UpstreamFreqFile = os.path.join(OutDir,"UpstreamFreq.txt")
    DownstreamFreqFile = os.path.join(OutDir,"DownstreamFreq.txt")
    PAMPHLETResources.write_freq_dict_to_file(UpstreamFreqList,UpstreamFreqFile)
    PAMPHLETResources.write_freq_dict_to_file(DownstreamFreqList,DownstreamFreqFile)

    ### NOW, DRAW WEBLOGO AND SAVE TO FILE
    UpstreamLogoFile = os.path.join(OutDir,"UpstreamLogo.png")
    DownstreamLogoFile = os.path.join(OutDir,"DownstreamLogo.png")
    PAMPHLETResources.draw_weblogo(UpstreamFreqFile,UpstreamLogoFile,"Upstream",FlankLength)
    PAMPHLETResources.draw_weblogo(DownstreamFreqFile,DownstreamLogoFile,"Downstream",FlankLength)

    ### PAMPHLET FINISHED
    print("PAMPHLET finished.")
    print("Total time used: %s seconds" % (time.time() - runningTime))
    print("===============================================================")

if __name__ == "__main__":
    main()
