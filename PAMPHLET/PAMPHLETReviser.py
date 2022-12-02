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
from Bio import SeqIO

from PAMPHLET import PAMPHLETParams
from PAMPHLET import OnlineResources
from PAMPHLET import PAMPHLETResources

def main(spacerFile,repeatseq,proteinFile,outDir,inunique,inblast,inlen,refmode):

    if inblast is None:

        ### Now, reading protein length and check revise status
        proteinLength = len(str(SeqIO.read(proteinFile,'fasta').seq))
        OnlineBlastDir = os.path.join(outDir,'OnlineBlast/')
        OnlineBlastResult = os.path.join(OnlineBlastDir,'OnlineBlastResult.txt')
        PAMPHLETParams.check_directory(OnlineBlastDir)

        ### Now, running protein blastp to predict related genomes
        OnlineResources.run_online_blastp(proteinFile,OnlineBlastResult)
        if PAMPHLETParams.check_null_files(OnlineBlastResult):
            print("No related protein found. Spacer revise failed.")
            return False
    
    else:

        OnlineBlastResult = inblast
        proteinLength = inlen
        OnlineBlastDir = os.path.join(outDir,'OnlineBlast/')
        PAMPHLETParams.check_directory(OnlineBlastDir)
    
    ### Get significant hits and return putative source tax name
    SignificantHitResult = os.path.join(OnlineBlastDir,"SignificantHits.txt")
    PutativeSourceTaxBox = PAMPHLETResources.get_significant_blastp_hits(OnlineBlastResult,SignificantHitResult,proteinLength)
    if PAMPHLETParams.check_null_files(SignificantHitResult):
        print("No significant hits found. Spacer revise failed.")
        return False

    ### Get seq dump file
    SeqDumpFile = os.path.join(OnlineBlastDir,"SeqDump.txt")
    refseqInformation = OnlineResources.get_seq_dump_file(SignificantHitResult,SeqDumpFile)

    if PAMPHLETParams.check_null_files(SeqDumpFile):
        print("No seq dump file found. Spacer revise failed.")
        return False

    ### Build sequenceID to taxID dictionary
    ReferenceInfo = PAMPHLETResources.get_seqid2seqtaxid(refseqInformation,PutativeSourceTaxBox,refmode)

    ### Now, extract the related contig
    RelatedSeqDumpDir = os.path.join(outDir,'RelatedSeqDump/')
    PAMPHLETParams.check_directory(RelatedSeqDumpDir)
    PAMPHLETResources.extract_related_seqdump(SeqDumpFile,RelatedSeqDumpDir,ReferenceInfo)

    ### Now, run the MinCED to predict the related contig CRISPR array
    MinCEDDir = os.path.join(outDir,'MinCED/')
    PAMPHLETParams.check_directory(MinCEDDir)
    PAMPHLETResources.run_minced(RelatedSeqDumpDir,MinCEDDir)

    ### Check minCED result, if all failed, return False;
    if not PAMPHLETParams.check_minced_files(MinCEDDir):
        print("MinCED failed. Spacer revise failed.")
        return False

    ### Now, select the correct DRs based on the MinCED result and repeat file
    SelectedSpacerFilesDict = PAMPHLETResources.select_correct_dr(MinCEDDir,repeatseq)

    ### Now, check the selected spacer file dict, if the dict is empty, return False
    if len(SelectedSpacerFilesDict) == 0:
        print("No correct DR found. Spacer revise failed.")
        return False
    
    ### Now, revise the spacer file
    RevisedSpacerFile = os.path.join(outDir,'RevisedSpacer.txt')
    PAMPHLETResources.revise_spacer_file(spacerFile,MinCEDDir,SelectedSpacerFilesDict,RevisedSpacerFile,inunique)

    ### Return revised spacer file
    return RevisedSpacerFile

if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])