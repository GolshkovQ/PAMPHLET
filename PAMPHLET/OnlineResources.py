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

def run_spacer_blast(infile,outfile,indb):
    ### This online blastn function, which is based on NCBI BLAST, is used to predict the putative protospacer sources. 
    ### Search against viruses database, which taxid = 10239

    UpLoadQuery = quote(open(infile,'r').read())
    if indb == 'phage':
        arguments = "CMD=Put&PROGRAM=blastn&MEGABLAST=on&DATABASE=nt&ENTREZ_QUERY=txid10239[ORGN]&QUERY=" + UpLoadQuery
    else:
        arguments = "CMD=Put&PROGRAM=blastn&MEGABLAST=on&DATABASE=nt&ENTREZ_QUERY=NOT%20txid2759[ORGN]&QUERY=" + UpLoadQuery

    r = requests.put("https://blast.ncbi.nlm.nih.gov/Blast.cgi?"+arguments)
    RID = r.text.split("RID = ")[1].split()[0]

    print("SPACER BLAST PROGRAM RID: "+RID)

    rStatus = requests.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID="+RID)
    Status = rStatus.text.split("Status=")[1].split()[0]

    while Status != "READY":
        time.sleep(30)
        rStatus = requests.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID="+RID)
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

def get_seq_dump_file(insig,outdump):
    infodict = {}
    tempdict = {}
    duplicateBox = []
    genomedupbox = []

    genomepage = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id="
    proteinpage = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id="
    with open(insig,'r') as fa, open(outdump,'w') as fb:
        for faline in fa.readlines():
            if len(duplicateBox) == 10:
                break
            sbjctID = faline.split("\t")[1]
            protein2fasta = proteinpage + sbjctID + "&rettype=fasta"
            proteinfasta = requests.get(protein2fasta).text
            try:
                proteinsequence = proteinfasta.split("\n")[1]
            except:
                ### Into the next loop
                print("Unknown error in: "+sbjctID)
                continue
            if not proteinsequence in duplicateBox:
                URLCorrectFlag = True
                url5 = "https://www.ncbi.nlm.nih.gov/ipg/"+sbjctID
                session = requests.Session()
                r5 = session.get(url5)
                url = "https://www.ncbi.nlm.nih.gov/sviewer/ipg/ipg.cgi?"
                post_data = {'query':'omitHeader%3Dfalse%26wt%3Dxml%26indent%3Dtrue%26rows%3D50%26start%3D0%26q%3D'+ sbjctID +'%26sort%3Dpriority_i%2520asc'}
                post_data_1 = {'db':'ipg','solr':'1','id':sbjctID}
                session.post(url,data= post_data_1).text
                r6 = session.post(url,data=post_data).text
                session.cookies.clear()
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
                        genomefasta = requests.get(gi2fasta)
                        if ">" in genomefasta.text and not "Error" in genomefasta.text and not "error" in genomefasta.text:
                            fb.write(genomefasta.text[:-1])
                            fb.flush()
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
                genomefasta = requests.get(gi2fasta)
                if ">" in genomefasta.text and not "Error" in genomefasta.text and not "error" in genomefasta.text:
                    seqdumpcounter += 1
                    fb.write(genomefasta.text[:-1])
                    fb.flush()
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
