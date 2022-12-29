#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
PAMPHLET - PAM Prediction HomoLogous Enhancement Toolkit
Author: Chen Qi, Baitao Li, Lei Huang
University of Chinese Academy of Sciences, College of Life Sciences, Beijing, China
BGI Research, Shenzhen, China
Email: qichen@genomics.cn; libaitao@genomics.cn; huanglei@genomics.cn
'''

from pyfaidx import Fasta
import math
import os
from Bio import pairwise2
from Bio import SeqIO

def reverse_complete(seq):
    trantab = str.maketrans('ATCGatcg','TAGCtagc')
    return seq.translate(trantab)[::-1]

def near_crispr_array(indict,inid,instart,inend):
    ### indict format:
    ### {[contigid]:[[CRISPR_1_start,CRISPR_1_end],[CRISPR_2_start,CRISPR_2_end]...],...}
    if inid in indict:
        for units in indict[inid]:
            startpos = units[0]
            endpos = units[1]
            ### Inside CRISPR array
            if instart <= endpos and instart >= startpos:
                return True
            if inend <= endpos and instart >= endpos:
                return True
            if abs(inend-endpos) <= 50 or abs(inend-startpos) <= 50:
                return True
            if abs(instart-endpos) <= 50 or abs(instart-startpos) <= 50:
                return True
    return False

def get_significant_blastp_hits(blastres,significantres,proteinlen,incov,inident):
    ### This step, get significant hits
    ### The significant hits are defined as the hits with evalue <= 1e-5 & qcov >= 0.9 & pident >= 90.0

    linecount = 0
    ### Only check first 5 lines as the top hits candidates
    TopHit = []

    with open(blastres,'r') as fa, open(significantres,'w') as fb:
        for records in fa.readlines():
            if records.strip() != "" and not records.startswith("#"):
                linecount += 1
                Evalue = float(records.split("\t")[10])
                HitLength = int(records.split("\t")[3])
                QueryCoverage = HitLength/proteinlen
                QueryIdentity = float(records.split("\t")[2])
                SubjectID = records.split("\t")[1]

                if Evalue <= 1e-5 and QueryCoverage >= incov and QueryIdentity >= inident*100:   # if True:
                    fb.write(records)
                if Evalue <= 1e-30 and QueryCoverage >= 0.98 and QueryIdentity >= 98.0 and linecount <= 5: # if True:
                    TopHit.append(SubjectID)

    fa.close()
    fb.close()
    return TopHit

def get_seqid2seqtaxid(infodict,tophitbox,inmode):
    ### This step, get seqid2seqtaxid
    ### If already have tophitID, the tempdict select will under no-reference mode; else, the tempdict select will under reference mode, the reference is the tophit
    
    tempdict = {}
    finaldict = {}

    ### Build tophitbox_tax
    tophitbox_tax = []
    for wps in tophitbox:
        if wps in infodict:
            tophitbox_tax.append(infodict[wps][1])

    if len(tophitbox_tax) != 0 and inmode == "r":
        ### Now, use the reference mode, select the tophit
        for wpid, wpinfo in infodict.items():
            genomeid = wpinfo[0]
            taxid = wpinfo[1]
            if taxid in tophitbox_tax:
                finaldict[wpid] = wpinfo
    elif inmode == "nr":
        identifybox = []
        ### Now, use the non-reference mode, stat all tax duplication, then select the max one, if have same max, then select the first one
        for wpid, wpinfo in infodict.items():
            genomeid = wpinfo[0]
            taxid = wpinfo[1]
            if genomeid not in tempdict:
                tempdict[taxid] = 0
            tempdict[taxid] += 1
        maxvalue = max(tempdict.values())
        for taxid, value in tempdict.items():
            if value == maxvalue:
                identifybox.append(taxid)
        for wpid, wpinfo in infodict.items():
            genomeid = wpinfo[0]
            taxid = wpinfo[1]
            if taxid in identifybox:
                finaldict[wpid] = wpinfo
    elif inmode == "a":
        ### Now, use all mode, select all
        for wpid, wpinfo in infodict.items():
            genomeid = wpinfo[0]
            taxid = wpinfo[1]
            finaldict[wpid] = wpinfo

    return finaldict

def extract_related_seqdump(infile,outdir,refdict):
    refseqs = Fasta(infile)

    ### Now, build the hit position dict based on inhit
    duplicatebox = []
    for wpid, wpinfo in refdict.items():
        genomeid = wpinfo[0]
        start = int(wpinfo[2])
        end = int(wpinfo[3])
        if start >= end:
            start, end = end, start
        if not genomeid in duplicatebox and genomeid in refseqs.keys():
            with open(os.path.join(outdir,genomeid + ".fa"),'a') as fa:
                contiglen = len(str(refseqs[genomeid][::]))
                if start <= 20000:
                    start = 1
                else:
                    start = start - 20000
                if end + 20000 >= contiglen:
                    end = contiglen
                else:
                    end = end + 20000
                fa.write(">"+genomeid+"\n"+str(refseqs[genomeid][start:end])+"\n")
            fa.close()

    return True

def run_minced(indir,outdir):
    for fastafile in os.listdir(indir):
        if fastafile.endswith(".fa"):
            CMD = "minced -spacers "+os.path.join(indir,fastafile)+" "+\
                os.path.join(outdir,fastafile.replace(".fa",".minced"))+" "+\
                os.path.join(outdir,fastafile.replace(".fa",".gff"))
            print("Now running minced for "+fastafile+", command: "+CMD)
            os.system(CMD)
    return True

def select_correct_dr(indir,tempseq,inident):
    finaldict = {}
    for mincedfile in os.listdir(indir):
        if mincedfile.endswith(".gff"):
            with open(os.path.join(indir,mincedfile),'r') as fa:
                for faline in fa.readlines():
                    if not faline.startswith("#"):
                        rptunitseq = faline.strip().split("\t")[8].split("rpt_unit_seq=")[1].upper().replace("U","T")
                        rptunitid = faline.strip().split("\t")[8].split("ID=")[1].split(";")[0].replace("CRISPR","CRISPR_")
                        ### Now, check the similarity between the rptunitseq and the tempseq
                        topalign = pairwise2.align.globalms(rptunitseq,tempseq.replace("U","T"),10,0,-10,-0.5,one_alignment_only=True)
                        formatalign = pairwise2.format_alignment(*topalign[0])
                        similarity = formatalign.split("\n")[1].count("|")/len(topalign[0][0].lstrip("-").rstrip("-"))
                        if similarity >= inident:
                            finaldict[mincedfile.replace(".gff","_spacers.fa")] = ["P",rptunitid]
                        else:
                            rcrptunitseq = reverse_complete(rptunitseq)
                            topalign = pairwise2.align.globalms(rcrptunitseq,tempseq.replace("U","T"),10,0,-10,-0.5,one_alignment_only=True)
                            formatalign = pairwise2.format_alignment(*topalign[0])
                            similarity = formatalign.split("\n")[1].count("|")/len(topalign[0][0].lstrip("-").rstrip("-"))
                            if similarity >= inident:
                                finaldict[mincedfile.replace(".gff","_spacers.fa")] = ["N",rptunitid]
            fa.close()
    return finaldict

def revise_spacer_file(infile,spacerloc,spacerdict,outfile,inmode):
    if inmode:
        uniquebox = []
        with open(outfile,'w') as fb:
            for temprecords in SeqIO.parse(infile,'fasta'):
                if not temprecords.seq in uniquebox:
                    fb.write(">"+str(temprecords.id)+"\n"+str(temprecords.seq)+"\n")
                    uniquebox.append(str(temprecords.seq))
            for targetFiles in spacerdict.keys():
                if spacerdict[targetFiles][0] == "P":
                    spacerfile = os.path.join(spacerloc,targetFiles)
                    for reviserecords in SeqIO.parse(spacerfile,'fasta'):
                        targetCRISPRID = "_".join(str(reviserecords.id).replace(targetFiles.replace("_spacers.fa",""),"").split("_")[:-2])[1:]
                        if targetCRISPRID == spacerdict[targetFiles][1] and not str(reviserecords.seq) in uniquebox:
                            fb.write(">"+str(reviserecords.id)+"\n"+str(reviserecords.seq)+"\n")
                            uniquebox.append(str(reviserecords.seq))
                else:
                    spacerfile = os.path.join(spacerloc,targetFiles)
                    for reviserecords in SeqIO.parse(spacerfile,'fasta'):
                        targetCRISPRID = "_".join(str(reviserecords.id).replace(targetFiles.replace("_spacers.fa",""),"").split("_")[:-2])[1:]
                        if targetCRISPRID == spacerdict[targetFiles][1] and not reverse_complete(str(reviserecords.seq)) in uniquebox:
                            fb.write(">"+str(reviserecords.id)+"\n"+reverse_complete(str(reviserecords.seq))+"\n")
                            uniquebox.append(reverse_complete(str(reviserecords.seq)))
        fb.close()
    else:
        with open(outfile,'w') as fb:
            for temprecords in SeqIO.parse(infile,'fasta'):
                fb.write(">"+str(temprecords.id)+"\n"+str(temprecords.seq)+"\n")
            for targetFiles in spacerdict.keys():
                if spacerdict[targetFiles][0] == "P":
                    spacerfile = os.path.join(spacerloc,targetFiles)
                    for reviserecords in SeqIO.parse(spacerfile,'fasta'):
                        targetCRISPRID = "_".join(str(reviserecords.id).replace(targetFiles.replace("_spacers.fa",""),"").split("_")[:-2])[1:]
                        if targetCRISPRID == spacerdict[targetFiles][1]:
                            fb.write(">"+str(reviserecords.id)+"\n"+str(reviserecords.seq)+"\n")
                else:
                    spacerfile = os.path.join(spacerloc,targetFiles)
                    for reviserecords in SeqIO.parse(spacerfile,'fasta'):
                        targetCRISPRID = "_".join(str(reviserecords.id).replace(targetFiles.replace("_spacers.fa",""),"").split("_")[:-2])[1:]
                        if targetCRISPRID == spacerdict[targetFiles][1]:
                            fb.write(">"+str(reviserecords.id)+"\n"+reverse_complete(str(reviserecords.seq))+"\n")
        fb.close()
    return True

def split_spacer_file(infile,outdir,seqnum):
    splitpartnum = math.ceil(float(seqnum/15))
    ### Use seqkit to split the file
    CMD = "seqkit split -p "+str(splitpartnum)+" "+infile+" --quiet"
    os.system(CMD)
    splitteddir = os.path.join(infile+".split")
    os.system("mv "+os.path.join(splitteddir,"*")+" "+outdir)
    os.system("rm -r "+splitteddir)
    return True

def rename_spacer_id(infile,inlabel,inorientation):
    finalSpacerBox = {}
    spacercounter = 1
    with open(infile+".tmp",'w') as fa, open(inlabel,'w') as fb:
        for spacers in SeqIO.parse(infile,'fasta'):
            if inorientation == "positive":
                fa.write(">SpacerID"+str(spacercounter)+"\n"+str(spacers.seq)+"\n")
                fb.write(">SpacerID"+str(spacercounter)+"\t"+str(spacers.id)+"\n")
            else:
                fa.write(">SpacerID"+str(spacercounter)+"\n"+reverse_complete(str(spacers.seq))+"\n")
                fb.write(">SpacerID"+str(spacercounter)+"\t"+str(spacers.id)+"\n")
            finalSpacerBox["SpacerID"+str(spacercounter)] = len(str(spacers.seq))
            spacercounter += 1
    fa.close()
    fb.close()
    return finalSpacerBox

def merge_blast_output(indir,outfile):
    with open(outfile,'w') as fa:
        for files in os.listdir(indir):
            with open(os.path.join(indir,files),'r') as fb:
                for fblines in fb.readlines():
                    ### Check blank line
                    if not fblines.strip() == "":
                        fa.write(fblines)
            fb.close()
    fa.close()
    return True

def get_significant_spacer_blast_output(inraw,insig,indict):
    significant_box = {}
    insignificant_box = []
    with open(inraw,'r') as fa, open(insig,'w') as fb:
        for faline in fa.readlines():
            if not faline.strip() == "":
                queryid = faline.strip().split("\t")[0]
                mismatches = int(faline.strip().split("\t")[4])
                gapnumbers = int(faline.strip().split("\t")[5])
                evalue = float(faline.strip().split("\t")[-2])
                querystartpos = int(faline.strip().split("\t")[6])
                queryendpos = int(faline.strip().split("\t")[7])
                qlen = indict[queryid]
                qcovs = (queryendpos-querystartpos+1)/qlen
                if evalue <= 1 and qcovs == 1 and mismatches <= 3 and gapnumbers == 0:
                    if not queryid in significant_box:
                        significant_box[queryid] = 0
                    significant_box[queryid] += 1
                    fb.write(faline)
                else:
                    ### This step, check the spacer alignment region, if only start/end 1nt not aligned, then it is still significant. But the mismatch numbers should be less than 2.
                    if querystartpos == 2 or queryendpos == qlen-1:
                        if mismatches <= 2 and gapnumbers == 0:
                            if not queryid in significant_box:
                                significant_box[queryid] = 0
                            significant_box[queryid] += 1
                            fb.write(faline)
    for spacerid in indict.keys():
        if not spacerid in significant_box or significant_box[spacerid] == 0:
            insignificant_box.append(spacerid)
    return insignificant_box

def revise_spacer_length(infile,outfile,inid):
    refseqs = Fasta(infile)
    with open(outfile,'w') as fa:
        for seqid in inid:
            fa.write(">"+seqid+"\n"+str(refseqs[seqid][1:-1])+"\n")
    fa.close()
    return True

def merge_spacer_blast_output(fileA,fileB,outfile):
    ### fileA = significantHit, fileB = revisedSignificantHit
    duplicateBox = []
    with open(fileA,'r') as fa, open(fileB,'r') as fb, open(outfile,'w') as fc:
        for faline in fa.readlines():
            if not faline.strip() == "":
                if int(faline.strip().split("\t")[8]) > int(faline.strip().split("\t")[9]):
                    complexesC = faline.strip().split("\t")[0] + "-" + faline.strip().split("\t")[1] + "-" +str(int(faline.strip().split("\t")[9])+1) + "-" + str(int(faline.strip().split("\t")[8])-1)
                    complexesL = faline.strip().split("\t")[0] + "-" + faline.strip().split("\t")[1] + "-" +str(int(faline.strip().split("\t")[9])+1) + "-" + str(int(faline.strip().split("\t")[8]))
                    complexesR = faline.strip().split("\t")[0] + "-" + faline.strip().split("\t")[1] + "-" +str(int(faline.strip().split("\t")[9])) + "-" + str(int(faline.strip().split("\t")[8])-1)
                else:
                    complexesC = faline.strip().split("\t")[0] + "-" + faline.strip().split("\t")[1] + "-" +str(int(faline.strip().split("\t")[8])+1) + "-" + str(int(faline.strip().split("\t")[9])-1)
                    complexesL = faline.strip().split("\t")[0] + "-" + faline.strip().split("\t")[1] + "-" +str(int(faline.strip().split("\t")[8])+1) + "-" + str(int(faline.strip().split("\t")[9]))
                    complexesR = faline.strip().split("\t")[0] + "-" + faline.strip().split("\t")[1] + "-" +str(int(faline.strip().split("\t")[8])) + "-" + str(int(faline.strip().split("\t")[9])-1)
                duplicateBox.append(complexesC)
                duplicateBox.append(complexesL)
                duplicateBox.append(complexesR)
                fc.write(faline)
        for fbline in fb.readlines():
            if not fbline.strip() == "":
                if int(fbline.strip().split("\t")[8]) > int(fbline.strip().split("\t")[9]):
                    revisedcomplexes = fbline.strip().split("\t")[0] + "-" + fbline.strip().split("\t")[1] + "-" +str(fbline.strip().split("\t")[8])
                else:
                    revisedcomplexes = fbline.strip().split("\t")[0] + "-" + fbline.strip().split("\t")[1] + "-" +str(fbline.strip().split("\t")[9])
                revisedcomplexes = fbline.strip().split("\t")[0] + "-" + fbline.strip().split("\t")[1] + "-" +str(fbline.strip().split("\t")[8]) + "-" + str(fbline.strip().split("\t")[9])
                if not revisedcomplexes in duplicateBox:
                    fc.write(fbline)
    fa.close()
    fb.close()
    fc.close()
    return True

def build_revised_spacer_len_dict(indict,inlist):
    finaldict = {}
    for units in inlist:
        if units in indict.keys():
            finaldict[units] = indict[units]-2
    return finaldict

def get_flank_seq(inhit,ingenome,inupstream,indownstream,flanklen,tempdir):
    refseqs = Fasta(ingenome)
    upstreamDict = {}
    downstreamDict = {}

    ### Run MinCED vs. ingenome, generate CRISPR gff files and build CRISPR gff dict
    CRISPRDict = {}
    CMD = "minced "+ingenome+" "+os.path.join(tempdir,os.path.basename(ingenome)+".crispr")+" "+os.path.join(tempdir,os.path.basename(ingenome)+".gff")
    os.system(CMD)
    print("Running MinCED to generate CRISPR gff file, command: "+CMD)
    with open(os.path.join(tempdir,os.path.basename(ingenome)+".gff"),'r') as fg:
        for fgline in fg.readlines():
            if not fgline.startswith("#"):
                contigID = fgline.split("\t")[0]
                CRISPRstart = int(fgline.split("\t")[3])
                CRISPRend = int(fgline.split("\t")[4])
                if not contigID in CRISPRDict:
                    CRISPRDict[contigID] = []
                CRISPRDict[contigID].append([CRISPRstart,CRISPRend])
    
    ### Now, start extract flanking sequence, if sbjctstartpos & sbjctendpos inside/near 50bp of a CRISPR array, then drop this hit.
    with open(inhit,'r') as fa, open(inupstream,'w') as fb, open(indownstream,'w') as fc:
        for faline in fa.readlines():
            queryid = faline.strip().split("\t")[0]
            sbjctid = faline.strip().split("\t")[1]
            sbjctstartpos = int(faline.strip().split("\t")[8])
            sbjctendpos = int(faline.strip().split("\t")[9])
            if sbjctid in refseqs.keys():
                ### Check the strand
                if not near_crispr_array(CRISPRDict,sbjctid,sbjctstartpos,sbjctendpos):
                    if sbjctstartpos < sbjctendpos:
                        strand = "+"
                    else:
                        strand = "-"
                    
                    ### Check the genome length, if not enough to extract flank (flanklen) nt sequence, then skip
                    if strand == "+":
                        if sbjctstartpos > flanklen and len(str(refseqs[sbjctid][::])) - sbjctendpos > flanklen:
                            upstreamSeqs = str(refseqs[sbjctid][sbjctstartpos-flanklen-1:sbjctstartpos-1])
                            downstreamSeqs = str(refseqs[sbjctid][sbjctendpos:sbjctendpos+flanklen])
                            fb.write(">"+queryid+"#"+sbjctid+"#upstream"+"\n"+str(upstreamSeqs)+"\n")
                            fc.write(">"+queryid+"#"+sbjctid+"#downstream"+"\n"+str(downstreamSeqs)+"\n")
                            if not queryid in upstreamDict.keys():
                                upstreamDict[queryid] = []
                            upstreamDict[queryid].append(upstreamSeqs)
                            if not queryid in downstreamDict.keys():
                                downstreamDict[queryid] = []
                            downstreamDict[queryid].append(downstreamSeqs)
                    else:
                        if sbjctendpos > flanklen and len(str(refseqs[sbjctid][::])) - sbjctstartpos > flanklen:
                            upstreamSeqs = reverse_complete(str(refseqs[sbjctid][sbjctstartpos:sbjctstartpos+flanklen]))
                            downstreamSeqs = reverse_complete(str(refseqs[sbjctid][sbjctendpos-flanklen-1:sbjctendpos-1]))
                            fb.write(">"+queryid+"#"+sbjctid+"#upstream"+"\n"+str(upstreamSeqs)+"\n")
                            fc.write(">"+queryid+"#"+sbjctid+"#downstream"+"\n"+str(downstreamSeqs)+"\n")
                            if not queryid in upstreamDict.keys():
                                upstreamDict[queryid] = []
                            upstreamDict[queryid].append(upstreamSeqs)
                            if not queryid in downstreamDict.keys():
                                downstreamDict[queryid] = []
                            downstreamDict[queryid].append(downstreamSeqs)
    fa.close()
    fb.close()
    fc.close()
    return upstreamDict,downstreamDict

def convert_seq_to_freq(indict,flanklen,freqmode):
    ### This function, will return a list, each element is a dictionary and key is A/T/C/G and the value is the score.
    ### The score is the sum of each position of the base frequency.

    ### This step, will build a dictionary, key is the spacer id and the value is a list, which contains the frequency of each position.
    ### The list is a list of dictionary, each dictionary is the frequency of each base at each position.
    freqDictTemp = {}
    freqFinal = []
    for seqid in indict.keys():
        seqIDfreqBox = []
        for base_position in range(0,flanklen):
            base_freq_dict = {}
            counterA = 0
            counterT = 0
            counterG = 0
            counterC = 0
            for seqs in indict[seqid]:
                if seqs[base_position] == "A":
                    counterA += 1
                elif seqs[base_position] == "T":
                    counterT += 1
                elif seqs[base_position] == "G":
                    counterG += 1
                elif seqs[base_position] == "C":
                    counterC += 1
            freqA = counterA/len(indict[seqid])
            freqT = counterT/len(indict[seqid])
            freqG = counterG/len(indict[seqid])
            freqC = counterC/len(indict[seqid])
            base_freq_dict['A'] = freqA
            base_freq_dict['T'] = freqT
            base_freq_dict['G'] = freqG
            base_freq_dict['C'] = freqC
            seqIDfreqBox.append(base_freq_dict)
        freqDictTemp[seqid] = seqIDfreqBox
    
    ### Now, calculate each base position each base score, the score is the sum of all the same base frequency
    for pos in range(0,flanklen):
        TotalCount = 0
        sumAfreq = 0
        sumTfreq = 0
        sumGfreq = 0
        sumCfreq = 0
        for seqid in freqDictTemp.keys():
            TotalCount += 1
            sumAfreq += freqDictTemp[seqid][pos]['A']
            sumTfreq += freqDictTemp[seqid][pos]['T']
            sumGfreq += freqDictTemp[seqid][pos]['G']
            sumCfreq += freqDictTemp[seqid][pos]['C']

        if freqmode == "sigmoid":
            ### Now, calculate the average frequency of each base at each position, normalized by sigmoid function
            sigmoid_normalized_A = sumAfreq*(1/(1+math.exp(-TotalCount)))
            sigmoid_normalized_T = sumTfreq*(1/(1+math.exp(-TotalCount)))
            sigmoid_normalized_G = sumGfreq*(1/(1+math.exp(-TotalCount)))
            sigmoid_normalized_C = sumCfreq*(1/(1+math.exp(-TotalCount)))
            freqFinal.append({'A':sigmoid_normalized_A,'T':sigmoid_normalized_T,'G':sigmoid_normalized_G,'C':sigmoid_normalized_C})
        else:
            freqFinal.append({'A':sumAfreq,'T':sumTfreq,'G':sumGfreq,'C':sumCfreq})

    return freqFinal

def write_freq_dict_to_file(indict,infile):
    with open(infile,'w') as fa:
        fa.write("PO\tA\tT\tG\tC\n")
        for positions in range(0,len(indict)):
            fa.write(str(positions+1)+"\t"+str(indict[positions]['A'])+"\t"+str(indict[positions]['T'])+"\t"+str(indict[positions]['G'])+"\t"+str(indict[positions]['C'])+"\n")
    fa.close()
    return True

def draw_weblogo(infile,outpic,flag,flankLen):
    if flag == "Upstream":
        xannotation = ",".join([str(i) for i in range(-1*flankLen,0)])
        weblogotitle = "Upstream"
    else:
        xannotation = ",".join([str(i) for i in range(1,flankLen+1)])
        weblogotitle = "Downstream"
    CMD = "weblogo -f "+infile+" -o "+outpic+" -F jpeg --title "+weblogotitle+" --size large --annotate "+xannotation+" --resolution 600 --color blue C \'C\' --color red T \'T\' --color green A \'A\' --color orange G \'G\'"
    os.system(CMD)

def main():
    print("Resources test.")

if __name__ == '__main__':
    main()
