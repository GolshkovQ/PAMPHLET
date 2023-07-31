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
import os
import requests
import time
import func_timeout
from func_timeout import func_set_timeout
from urllib.parse import quote
from pyfaidx import Fasta

def run_online_blastp(infile,outfile):
    ### This online blastp function, which is based on NCBI BLAST, is used to predict the homolog protein-related genome sequences.
    ### Search against Bacteria database, which taxid = 2

    ErrorMessage = "CPU usage limit was exceeded"

    upquery = open(infile,'r').read()

    UpLoadQuery = quote(upquery)
    arguments = "CMD=Put&PROGRAM=blastp&DATABASE=nr&ENTREZ_QUERY=txid2[ORGN]&EXPECT=0.00001&QUERY=" + UpLoadQuery

    r = requests.put("https://blast.ncbi.nlm.nih.gov/Blast.cgi?"+arguments)
    RID = r.text.split("RID = ")[1].split()[0]
    print("PROTEIN blastp PROGRAM RID: "+RID)

    rStatus = requests.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID="+RID)
    Status = rStatus.text.split("Status=")[1].split()[0]

    while Status != "READY":
        time.sleep(30)
        rStatus = requests.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID="+RID)
        Status = rStatus.text.split("Status=")[1].split()[0]
        print("Refresh in 30s, status: "+Status)
    
    rResult = requests.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_TYPE=Tabular&RID="+RID)

    ### Check rResult error message, if error message is "CPU usage limit was exceeded", then re-run the blastp program.
    if ErrorMessage in rResult.text:
        print("NCBI blastp CPU usage limit was exceeded, re-run blastp program.")
        run_online_blastp(infile,outfile)
    
    else:
        with open(outfile+".temp",'w') as fot:
            fot.write(rResult.text)
        fot.close()
        flag = 0

        with open(outfile+".temp",'r') as fotr, open(outfile,'w') as foc:
            for fotrecord in fotr.readlines():
                if flag == 1:
                    if fotrecord.strip() != "</PRE>":
                        foc.write(fotrecord)
                if "hits found" in fotrecord:
                    flag = 1
                if fotrecord.strip() == "</PRE>":
                    flag = 0

        os.system("rm "+outfile+".temp")
        return True

def run_spacer_blast(infile,outfile,indb,blastmode):
    ### This online blastn function, which is based on NCBI BLAST, is used to predict the putative protospacer sources. 
    ### Search against viruses database, which taxid = 10239

    UpLoadQuery = quote(open(infile,'r').read())
    if indb == 'phage':
        if blastmode == 'relax':
            arguments = "CMD=Put&PROGRAM=blastn&MEGABLAST=on&NUCL_REWARD=1&WORD_SIZE=16&GAPCOSTS=5 2&NUCL_PENALTY=-1&DATABASE=nt&ENTREZ_QUERY=txid10239[ORGN]&QUERY=" + UpLoadQuery
        else:
            arguments = "CMD=Put&PROGRAM=blastn&MEGABLAST=on&DATABASE=nt&ENTREZ_QUERY=txid10239[ORGN]&QUERY=" + UpLoadQuery
    else:
        if blastmode == 'relax':
            arguments = "CMD=Put&PROGRAM=blastn&MEGABLAST=on&NUCL_REWARD=1&WORD_SIZE=16&GAPCOSTS=5 2&NUCL_PENALTY=-1&DATABASE=nt&ENTREZ_QUERY=NOT%20txid2759[ORGN]&QUERY=" + UpLoadQuery
        else:
            arguments = "CMD=Put&PROGRAM=blastn&MEGABLAST=on&DATABASE=nt&ENTREZ_QUERY=NOT%20txid2759[ORGN]&QUERY=" + UpLoadQuery

    r = requests.put("https://blast.ncbi.nlm.nih.gov/Blast.cgi?"+arguments)
    RID = r.text.split("RID = ")[1].split()[0]
    while RID == "":
        print("RID extraction failed, re-searching...")
        time.sleep(2)
        RID = r.text.split("RID = ")[1].split()[0]

    print("SPACER BLAST PROGRAM RID: "+RID)

    rStatus = requests.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID="+RID)
    Status = rStatus.text.split("Status=")[1].split()[0]

    while Status != "READY":
        time.sleep(30)
        rStatus = requests.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID="+RID)
        Status = rStatus.text.split("Status=")[1].split()[0]
        Status = rStatus.text.split("Status=")[1].split()[0]
        print("Refresh in 30s, status: "+Status)

    rResult = requests.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_TYPE=Tabular&RID="+RID)

    with open(outfile+".temp",'w') as fot:
        fot.write(rResult.text)
    fot.close()
    flag = 0

    with open(outfile+".temp",'r') as fotr, open(outfile,'w') as foc:
        for fotrecord in fotr.readlines():
            if flag == 1:
                if fotrecord.strip() != "</PRE>" and not fotrecord.startswith("#"):
                    foc.write(fotrecord)
            if "hits found" in fotrecord:
                flag = 1
            if fotrecord.strip() == "</PRE>":
                flag = 0

    os.system("rm "+outfile+".temp")
    return True

### Set timeout for get_session_post_text function, time out limit is 10s
@func_set_timeout(10)
def get_session_post_text(inid):
    urls = "https://www.ncbi.nlm.nih.gov/ipg/"+inid
    session = requests.Session()
    r5 = session.get(urls)
    url = "https://www.ncbi.nlm.nih.gov/sviewer/ipg/ipg.cgi?"
    post_data = {'query':'omitHeader%3Dfalse%26wt%3Dxml%26indent%3Dtrue%26rows%3D50%26start%3D0%26q%3D'+ inid +'%26sort%3Dpriority_i%2520asc'}
    post_data_1 = {'db':'ipg','solr':'1','id':inid}
    session.post(url,data= post_data_1).text
    r6 = session.post(url,data=post_data).text
    session.cookies.clear()
    return r6

@func_set_timeout(20)
def get_protein_sequences(inurl):
    proteintext = requests.get(inurl).text
    return proteintext

@func_set_timeout(20)
def get_genome_fasta(inurl):
    genometext = requests.get(inurl).text
    return genometext

def get_seq_dump_file(insig,outdump,inmaxsize):
    infodict = {}
    tempdict = {}
    duplicateBox = []
    genomedupbox = []

    genomepage = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id="
    proteinpage = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id="
    with open(insig,'r') as fa, open(outdump,'w') as fb:
        for faline in fa.readlines():
            if len(duplicateBox) == inmaxsize:
                break
            sbjctID = faline.split("\t")[1]
            protein2fasta = proteinpage + sbjctID + "&rettype=fasta"
            retry_counter = 0
            while retry_counter < 3:
                try:
                    proteinfasta = get_protein_sequences(protein2fasta)
                    break
                except:
                    retry_counter += 1
                    print("Retry in: "+sbjctID+", retry time: "+str(retry_counter))
            try:
                proteinsequence = "".join(proteinfasta.split("\n")[1:])
            except:
                ### Into the next loop
                print("Unknown error in: "+sbjctID)
                continue
            if not proteinsequence in duplicateBox:
                URLCorrectFlag = True
                retry_counter = 0
                ### Set timeout for get_session_post_text function, time out limit is 10s
                while retry_counter < 3:
                    try:
                        r6 = get_session_post_text(sbjctID)
                        break
                    except:
                        retry_counter += 1
                        print("Retry in: "+sbjctID+", retry time: "+str(retry_counter))
                genomerefheader = "<str name=\"cds_accver_s\">"
                genome_tax_header = "<str name=\"org_s\">"
                orfstartheader = "<long name=\"cds_start_l\">"
                orfendheader = "<long name=\"cds_stop_l\">"
                try:
                    refid = r6.split(genomerefheader)[1].split("<")[0]
                    genometaxid = r6.split(genome_tax_header)[1].split("<")[0].split("|")[1]
                    orfstart = r6.split(orfstartheader)[1].split("<")[0]
                    orfend = r6.split(orfendheader)[1].split("<")[0]
                except:
                    URLCorrectFlag = False
                if URLCorrectFlag:
                    duplicateBox.append(proteinsequence)
                    print("Putative homolog protein : "+sbjctID+", source: "+refid+", taxid: "+genometaxid+", orf start: "+orfstart+", orf end: "+orfend)
                    if not sbjctID in infodict.keys() and not refid in genomedupbox:
                        infodict[sbjctID] = [refid,genometaxid,int(orfstart),int(orfend)]
                        genomedupbox.append(refid)
                        tempdict[refid] = sbjctID
                        gi2fasta = genomepage+refid+"&rettype=fasta"
                        ### retry for 3 times
                        retry_counter = 0
                        while retry_counter < 3:
                            try:
                                genomefasta = get_genome_fasta(gi2fasta)
                                if ">" in genomefasta and not "Error" in genomefasta and not "error" in genomefasta:
                                    fb.write(genomefasta[:-1])
                                    fb.flush()
                                    break
                            except:
                                retry_counter += 1
                                print("Time out in: "+sbjctID+", retry in "+str(retry_counter)+" time(s)...")
                                continue
                else:
                    print("Failed ID convert for: "+sbjctID)
        fa.close()
        fb.close()

        ### Now, check fasta file
        with open(outdump,'r') as fa, open(outdump+".tmp",'w') as fb:
            for faline in fa.readlines():
                if faline.strip() == "" or "error" in faline or "Error" in faline:
                    continue
                else:
                    fb.write(faline)
        fa.close()
        fb.close()

        os.system("mv "+outdump+".tmp "+outdump)

        refseqs = Fasta(outdump)
        for refseq in refseqs.keys():
            if refseq.split(".")[0] in tempdict.keys():
                targetsbj = tempdict[refseq.split(".")[0]]
                ### Now, revise the infodict refid
                infodict[targetsbj][0] = refseq

    return infodict

@func_set_timeout(10)
def get_genome_request_text(inurl):
    r = requests.get(inurl)
    return r.text

def get_genome_seqdump_files(inhit,outfile):
    seqdumplimit = 500
    seqdumpcounter = 0
    if seqdumpcounter >= 298:
        print("Seqdump numbers is reaching the limit, if you need a larger number, please modify the code.")
    genomepage = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id="
    duplicatebox = []
    with open(inhit,'r') as fa, open(outfile,'w') as fb:
        for faline in fa.readlines():
            sbjctid = faline.split()[1]
            if not sbjctid in duplicatebox:
                if seqdumpcounter >= seqdumplimit:
                    break
                duplicatebox.append(sbjctid)
                gi2fasta = genomepage + sbjctid + "&rettype=fasta"
                retry_counter = 0
                ### if false, retry 3 times
                while retry_counter < 3:
                    try:
                        genomefasta = get_genome_request_text(gi2fasta)
                        if ">" in genomefasta and not "Error" in genomefasta and not "error" in genomefasta:
                            seqdumpcounter += 1
                            fb.write(genomefasta[:-1])
                            break
                        else:
                            continue
                    except:
                        retry_counter += 1
                        print("Time out in: "+sbjctid+", retry for "+str(retry_counter)+" time(s)...")
                        continue
            else:
                continue
    fa.close()
    fb.close()
    return True

def main(infile,outfile,moder):
    if moder == "blastp":
        run_online_blastp(infile,outfile)
    else:
        run_spacer_blast(infile,outfile)

if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2],sys.argv[3])
