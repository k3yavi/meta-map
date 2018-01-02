from __future__ import print_function
import pysam
import pandas as pd
import numpy as np
import sys
import os
import math
import click

mil = 1000000
hund = 100

def get_ref_id(aln):
    # get id for the mapped reference
    try:
        return aln.reference_name.split("|")[1]
    except:
        return aln.reference_name.split(".")[0]

def print_details(qId, rId, aln, id2phlm, lvl="sing"):
    if lvl == "sing":
        print ("ERROR", file=sys.stderr)
        print("SAM Algns", file=sys.stderr)
        print (aln.query_name, aln.reference_name, file=sys.stderr)
        print("SUBSETTING", file=sys.stderr)
        print (qId, rId, file=sys.stderr)
        print ("PHYLA", file=sys.stderr)
        print (id2phlm[qId], id2phlm[rId], file=sys.stderr, end="\n")
    else:
        print ("ERROR", file=sys.stderr)
        print ("PHYLA", file=sys.stderr)
        print (qId, file=sys.stderr, end="\t")
        for an in aln:
            print (an, file=sys.stderr, end="\t")
        print("\n", file=sys.stderr)


def write_stats(totCount, singCount, totReads, roseCount, euCount, orphanCount, skipCount, TP, FP, FN, TN, cwd):
    mmCount = round(totCount - singCount)
    unmapCount = totReads - totCount
    FN = totReads - TP - FP - TN - roseCount - euCount


    sen = TP/float(TP+FN)
    spec = TN/float(TN+FP)
    ppv = TP/float(TP+FP)
    npv = TN/float(TN+FN)
    mcc = ((TP*TN)-(FP*FN)) / math.sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)  )
    stText = "\n\n" + \
    "====================================================================================\n" + \
    "Total Number of reads: {0} ({1:.2f}M)\n".format(totReads, totReads/mil)+ \
    "Number of Unmapped reads: {0} ({1:.2f}M, {2:.2f}%)\n".format(unmapCount, unmapCount/mil, unmapCount*hund/totReads)+ \
    "Number of Mapped reads {0}({1:.2f}M, {2:.2f}%)\n".format(totCount, totCount/mil, totCount*hund/totReads)+ \
    "\n\n"+ \
    "============================ OUT OF MAPPED READS ===================================\n"+ \
    "Number of Singly Mapped reads: {0} ({1:.2f}M, {2:.2f}%)\n".format(singCount, singCount/mil, singCount*hund/totReads)+ \
    "Number of Multimapped reads: {0} ({1:.2f}M, {2:.2f}%)\n".format(mmCount, mmCount/mil, mmCount*hund/totReads)+ \
    "Number of Mapped but Skipped reads: {0} ({1:.2f}M, {2:.2f}%)\n".format(skipCount, skipCount/mil, skipCount*hund/totReads)+ \
    "Number of Orphaned (Ignored)ALIGNMENTS (Should be significantly low): {}\n".format(orphanCount)+ \
    "====================================================================================\n"+ \
    "\n ===================== \n ASSUMPTIONS \n =====================\n"+ \
    "1: Any Multi-mapped read has the Original Phyla in ATLEAST 1 alignment\n"+ \
    "2: Don't know what to do with Rose Sequence Ignoring for now\n"+ \
    "3: No reference of Eukaryotes added\n"+ \
    "OVERALL: Atmost 10% Reads could have been mapped more.\n"+ \
    "Eukaryotes Counts: {0}({1:.2f}%)\n".format(euCount, euCount*hund/totReads)+ \
    "Rose Counts: {0}({1:.2f}%)\n".format(roseCount, roseCount*hund/totReads)+ \
    "====================================================================================\n"+ \
    "\n ===================== \n ACCURACY METRIC \n =====================\n"+ \
    "Number of True positives(TP) reads: {0} ({1:.2f}M, {2:.2f}%)\n".format(TP, TP/mil, TP*hund/totReads)+ \
    "Number of False Negatives(FN) reads: {0} ({1:.2f}M, {2:.2f}%)\n".format(FN, FN/mil, FN*hund/totReads)+ \
    "Number of False positives(FP) reads: {0} ({1:.2f}M, {2:.2f}%)\n".format(FP, FP/mil, FP*hund/totReads)+ \
    "Number of True Negatives(TN) reads: {0} ({1:.2f}M, {2:.2f}%)\n".format(TN, TN/mil, TN*hund/totReads)+ \
    "sens:{0:.2f} spec:{1:.2f} prec:{2:.2f} mcc:{3:.2f}\n".format(sen, spec, ppv, mcc)+ \
    "====================================================================================\n\n\n\n"

    filename = cwd + "/report.txt"
    with open(filename, 'w') as f:
        f.write(stText)


@click.command()
@click.option('--sam',  help='path for the sam file generated by pufmap')
@click.option('--fq',  help='path for the fastq file')
@click.option('--level',  help='all/some alignments to use for comparison', default="all")
def get_stats(sam, fq, level):
    cwd = os.getcwd()
    totReads = 0
    TN = 0
    FN = 0
    reads_list = []
    roseCount = 0
    euCount = 0

    with open("/mnt/scratch2/avi/meta-map/reads/meta/s2.tsv") as f:
        data = pd.read_table(f).set_index('EMBL ID').drop(["Strain/species details"], axis=1).to_dict()["Phylum"]
    id2phlm = {}
    for k,v in data.items():
        nk = k.split(".")[0]
        id2phlm[nk] = v
        # one of the entry got two identifiers
        if (nk == "CM000636"):
            id2phlm["CP006835"] = v
        elif v == "Rhizobium_Bradyrhizobium":
            id2phlm[nk] = "Proteobacteria"
        elif v == "Pathogens":
            id2phlm[nk] = "Proteobacteria"
    del data

    with open(fq) as f:
        for line in f:
            # counting total reads
            totReads += 1
            # Progress Monitoring
            if(totReads % mil == 0):
                print ("\r Done reading {} Million reads from fastq.".format(int(round(totReads)/1000000)), end="")
                sys.stdout.flush()

            #extracting relevant part of read
            read = line.strip().replace("/1","").replace("@","")

            if "Random" in read:
                TN += 1
            if "Eukaryotes" in read:
                euCount += 1
            if "Rose" in read:
                roseCount += 1

            #making a list of read id
            reads_list.append(read)

            # skip next 4 lines
            for _ in range(3):
                f.next()

    if len(reads_list) != len(set(reads_list)):
        print ("ERROR: Repeating reads found")
        return 0

    print ("Done reading fastq.")
    with pysam.AlignmentFile(sam) as f:
        TP = 0
        FP = 0
        totCount = 0.0
        singCount = 0
        orphanCount = 0
        skipCount = 0.0
        for aln in f:
            #get mate of the read
            mate_aln = f.next()

            # count total Number of reads
            totCount += 1

            # get number of alignments
            n_alns = aln.get_tag('NH')

            #ignoring Rose Sequence
            if "Rose" in aln.query_name or "Eukaryotes" in aln.query_name:
                skipCount += 1.0/n_alns
                continue

            # Ignoring Orphan alignments for now
            if(aln.reference_name != mate_aln.reference_name):
                orphanCount += 1
                print ("WARNING: ORPHANS Detected statistics Needs to be re-evaluated")
                continue

            # Progress Monitoring
            if(round(totCount) % mil == 0):
                print ("\r Done reading {} Million reads from SAM.".format(int(round(totCount)/1000000)), end="")
                sys.stdout.flush()

            # get ground truth id
            qId = aln.query_name.split('-')[0]
            if "|" in qId:
                qId = qId.split("|")[1]
            elif "_" in qId:
                qId = qId.split("_")[0]

            # for singly mapped reads only
            if n_alns == 1:
                # Increment the single count
                singCount += 1

                # get id for the mapped reference
                rId = get_ref_id(aln)

                # compare the labels
                try:
                    if qId == rId or id2phlm[qId]==id2phlm[rId]:
                        # get TP count
                        TP += 1
                    else:
                        print_details(qId, rId, aln, id2phlm)
                        # get FP count
                        FP += 1
                except:
                    print ("In Single mapping")
                    print_details(qId, rId, aln, id2phlm)


            #handling MultiMapped Reads
            else:
                # list of all alignments
                algns = [get_ref_id(aln)]

                # iterate over all alignments
                for _ in range(1, n_alns):
                    aln = f.next()
                    mate_aln = f.next()

                    # Ignoring Orphan alignments for now
                    if(aln.reference_name != mate_aln.reference_name):
                        orphanCount += 1
                    else:
                        algns.append(get_ref_id(aln))

                # SIMULATED FALSE POSITIVES
                if "Random" in aln.query_name:
                    TN -= 1
                    FP += 1
                    continue

                flag = False
                if level == "some":
                    for rId in algns:
                        try:
                            if qId == rId or id2phlm[qId]==id2phlm[rId]:
                                flag = True
                                break
                        except:
                            print ("In Multimapping")
                            print_details(qId, rId, aln, id2phlm)
                elif level == "all":
                    plist = set([])
                    try:
                        qId = id2phlm[qId]
                        for rId in algns:
                            plist.add( id2phlm[rId] )
                            if len(plist) > 1:
                                break
                    except:
                        print ("In Multimapping")
                        print_details(qId, rId, aln, id2phlm)

                    if(len(plist) == 1 and list(plist)[0] == qId):
                        flag=True

                if flag:
                    TP += 1
                else:
                    #print_details(qId, rId, plist, id2phlm, "multi")
                    FP += 1

    write_stats(totCount, singCount, totReads, roseCount, euCount, orphanCount, skipCount, TP, FP, FN, TN, cwd)
    print ("\noutput written to {}".format(os.getcwd()+"/reports.txt"))

if __name__=="__main__":
    get_stats()
