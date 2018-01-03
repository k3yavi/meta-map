from __future__ import print_function
import pysam
import pandas as pd
import os
import math
import click
import sys

mil = 1000000
hund = 100

def populate_phylum_dict(pname):
    with open(pname) as f:
        data = pd.read_table(f).set_index('EMBL ID').drop(["Strain/species details"], axis=1).to_dict()["Phylum"]
    id2phlm = {}
    for k,v in data.items():
        nk = k.split(".")[0]
        id2phlm[nk] = v
        if (nk == "CM000636"):
            id2phlm["CP006835"] = v
        elif v == "Rhizobium_Bradyrhizobium":
            id2phlm[nk] = "Proteobacteria"
        elif v == "Pathogens":
            id2phlm[nk] = "Proteobacteria"
    id2phlm["Rose"] = "Rose"
    id2phlm["Eukaryotes"] = "Eukaryotes"
    del data
    return id2phlm

def get_ref_id(aln):
    # get id for the mapped reference
    rname = aln.reference_name

    euList = ["Arabidopsis", "Human", "Lizard", "Chicken", "Eagle", "Turtle", "Yeast"]
    for eu in euList:
        if eu in rname:
            return "Eukaryotes"

    if "Rose" in rname:
        return "Rose"
    try:
        return aln.reference_name.split("|")[1]
    except:
        return aln.reference_name.split(".")[0]

def get_query_id(aln):
    qname = aln.query_name

    # get ground truth id
    if "Eukaryotes" in qname:
        return "Eukaryotes"
    elif "Rose" in qname:
        return "Rose"

    qId = qname.split('-')[0]
    if "|" in qId:
        qId = qId.split("|")[1]
    elif "_" in qId:
        qId = qId.split("_")[0]

    return qId

def print_details(qId, rIds, aln):
    print ("ALIGNMENT", file=sys.stderr)
    print (aln, file=sys.stderr )
    print ("PHYLA", file=sys.stderr)
    print ("QUERY:\t" + qId + "\tMAPPINGS:", end="\t", file=sys.stderr)
    for rId in rIds:
        print ( rId , end="\t", file=sys.stderr)
    print ("\n", file=sys.stderr)

def parse_fq(rname):
    totReads = 0
    TN = 0
    reads_list = []
    with open(rname) as f:
        for line in f:
            # counting total reads
            totReads += 1
            # Progress Monitoring
            if(totReads % mil == 0):
                print ("\r Done reading {} Million reads from FASTQ.".format(int(round(totReads)/1000000)), end="")
                sys.stdout.flush()

            #extracting relevant part of read
            read = line.strip().replace("/1","").replace("@","")

            if "Random" in read:
                TN += 1

            #making a list of read id
            reads_list.append(read)

            # skip next 4 lines
            for _ in range(3):
                f.next()

    if len(reads_list) != len(set(reads_list)):
        print ("ERROR: Repeating reads found")
        exit(1)
    return totReads, TN, reads_list

def print_stats(singCount, totCount, totReads, TP, FP, TN, orphanCount):
    mmCount = round(totCount - singCount)
    unmapCount = totReads - totCount
    FN = totReads - TP - FP - TN
    cwd = os.getcwd()

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
    "Number of Orphaned (Ignored)ALIGNMENTS (Should be significantly low): {}\n".format(orphanCount)+ \
    "====================================================================================\n"+ \
    "\n ===================== \n ACCURACY METRIC \n =====================\n"+ \
    "Number of True positives(TP) reads: {0} ({1:.2f}M, {2:.2f}%)\n".format(TP, TP/mil, TP*hund/totReads)+ \
    "Number of False Negatives(FN) reads: {0} ({1:.2f}M, {2:.2f}%)\n".format(FN, FN/mil, FN*hund/totReads)+ \
    "Number of False positives(FP) reads: {0} ({1:.2f}M, {2:.2f}%)\n".format(FP, FP/mil, FP*hund/totReads)+ \
    "Number of True Negatives(TN) reads: {0} ({1:.2f}M, {2:.2f}%)\n".format(TN, TN/mil, TN*hund/totReads)+ \
	"\n ===================== \n SECONDARY ACCURACY METRIC \n =====================\n"+ \
    "Senstivity: {}\n".format(sen)+ \
    "Specificity: {}\n".format(spec)+ \
    "Precision: {}\n".format(ppv)+ \
    "Neg Pred. Value: {}\n".format(npv)+ \
    "MCC: {}\n".format(mcc)+ \
    "====================================================================================\n\n\n\n"

    filename = cwd + "/report.txt"
    with open(filename, 'w') as f:
        f.write(stText)


def perform_counting(fname, totReads, TN, reads_list, id2phlm):
    with pysam.AlignmentFile(fname) as f:
        TP = 0
        FP = 0
        totCount = 0.0
        singCount = 0
        orphanCount = 0
        for aln in f:
            #get mate of the read
            mate_aln = f.next()

            # count total Number of reads
            totCount += 1

            # get number of alignments
            n_alns = aln.get_tag('NH')

            # for singly mapped reads only
            if n_alns == 1:
                # Increment the single count
                singCount += 1

            # Ignoring Orphan alignments for now
            if(aln.reference_name != mate_aln.reference_name):
                orphanCount += 1
                print ("WARNING: ORPHANS Detected statistics Neess to be re-evaluated")
                continue

            # Progress Monitoring
            if(round(totCount) % mil == 0):
                print ("\r Done reading {} Million reads from BAM....".format(int(round(totCount)/1000000)), end="")
                sys.stdout.flush()

            qId = get_query_id(aln)

            # list of all alignments
            rIds = [get_ref_id(aln)]

            # iterate over all alignments
            for _ in range(1, n_alns):
                aln = f.next()
                mate_aln = f.next()

                # Ignoring Orphan alignments for now
                if(aln.reference_name != mate_aln.reference_name):
                    orphanCount += 1
                else:
                    rIds.append(get_ref_id(aln))

            # skip the whole alignment list of it's a Random read
            if "Random" in aln.query_name:
                TN -= 1
                FP += 1
                continue

            plist = set([])
            try:
                qId_plm = id2phlm[qId]
                for rId in rIds:
                    rId_plm = id2phlm[rId]
                    plist.add( rId_plm )
                    if len(plist) > 1:
                        break
            except:
                print (qId, rIds)
                print_details(qId_plm, plist, aln)
                break

            if(len(plist) == 1 and list(plist)[0] == qId_plm):
                TP += 1
            else:
#                 print_details(qId_plm, plist, aln)
                FP += 1
    return singCount, totCount, TP, FP, TN, orphanCount

@click.command()
@click.option('--sam',  help='path for the sam file generated by pufmap')
@click.option('--fq',  help='path for the fastq file')
def get_stats(fq, sam):
    pname = "/mnt/scratch2/avi/meta-map/reads/meta/s2.tsv"

    # populate orgaism id to phylum dictionary
    id2phlm = populate_phylum_dict(pname)

    # parse the fastq file for TN calculations
    totReads, TN, reads_list = parse_fq(fq)

    # Parse the BAM and perform the counting
    singCount, totCount, TP, FP, TN, orphanCount = perform_counting(sam, totReads, TN, reads_list, id2phlm)

    # Calculate the stats and print it
    print_stats(singCount, totCount, totReads, TP, FP, TN, orphanCount)
    print ("\noutput written to {}".format(os.getcwd()+"/reports.txt"))

    return sen, spec, ppv, npv, mcc

if __name__=="__main__":
    get_stats()
