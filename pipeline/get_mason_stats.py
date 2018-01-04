from __future__ import print_function
import pandas as pd
import os
import click
import sys
from collections import defaultdict
import scipy.stats

mil = 1000000
hund = 100

def print_stats(singCount, totCount, orphanCount, phlmDict):
    mmCount = round(totCount - singCount)
    cwd = os.getcwd()

    stText = "\n\n" + \
    "====================================================================================\n" + \
        "Number of Mapped reads {0}({1:.2f}M)\n".format(totCount, totCount/mil)+ \
    "\n\n"+ \
    "============================ OUT OF MAPPED READS ===================================\n"+ \
    "Number of Singly Mapped reads: {0} ({1:.2f}M, {2:.2f}%)\n".format(singCount, singCount/mil, singCount*hund/totCount)+ \
    "Number of Multimapped reads: {0} ({1:.2f}M, {2:.2f}%)\n".format(mmCount, mmCount/mil, mmCount*hund/totCount)+ \
    "Number of Orphaned (Ignored)ALIGNMENTS (Should be significantly low): {}\n".format(orphanCount)+ \
    "====================================================================================\n\n\n\n"

    filename = cwd + "/report.txt"
    with open(filename, 'w') as f:
        f.write(stText)

    with open(cwd+"/dist.txt", 'w') as f:
        for k,v in phlmDict.items():
            f.write(k+","+str(v)+"\n")

def get_algn(aln):
    return aln.split()[2]

def perform_counting(fname, id2phlm):
    with open(fname) as f:
        totCount = 0.0
        singCount = 0
        orphanCount = 0
        phlmDict = defaultdict(int)
        for aln in f:
            toks = aln.strip().split()
            if toks[0][0] == "@":
                continue
            #get mate of the read
            #mate_aln = f.next()

            # get number of alignments
            n_alns = int(toks[-1].replace("NH:i:",""))

            # for singly mapped reads only
            if n_alns == 1:
                # Increment the single count
                singCount += 1
            elif n_alns == 0:
                continue

            # count total Number of reads
            totCount += 1

            rId = get_algn(aln)

            # list of all alignments
            rIds = set([ id2phlm[ rId  ] ])

            # Progress Monitoring
            if(round(totCount) % mil == 0):
                print ("\r Done reading {} Million reads from BAM....".format(int(round(totCount)/1000000)), end="")
                sys.stdout.flush()


            # iterate over all alignments
            for _ in range(1, n_alns):
                aln = f.next()

                rId = get_algn(aln)
                try:
                    rIds.add(id2phlm[ rId ])
                except:
                    print(aln)
                    print(n_alns)
                    exit(1)

            if len(rIds) == 1:
                phlmDict[list(rIds)[0]] += 1

    return singCount, totCount, orphanCount, phlmDict

def get_phlm_dict():
    id2phlm={}
    with open("/mnt/scratch2/avi/meta-map/kraken/meta/taxa.csv") as f:
            for line in f:
                toks = line.strip().split(",")
                val = toks[0]
                for i in range(1, len(toks)):
                    id2phlm[toks[i].strip(">")] = val
    return id2phlm

def get_correlation(phlmDict):
    with open("/mnt/scratch2/avi/meta-map/kraken/meta/truth_hc1.txt") as f:
        truth = pd.read_table(f).set_index("species").drop(["taxid","size","dataset"], 1)

    with open("/mnt/scratch2/avi/meta-map/kraken/meta/Huttenhower_HC1_Kraken.txt") as f:
        krak = pd.read_table(f, header=None).set_index(4).drop([0,2,3],1)

    print ("Truth Shape: # Reads(# phyla) {}({})".format(truth.sum().values[0], truth.shape[0]))
    print ("Kraken Shape: # Reads (# phyla) {}({})".format(krak.sum().values[0], krak.shape[0]))
    print ("Pufferfish Shape:  # Reads (# phyla) {}({})".format(sum(phlmDict.values()), len(phlmDict)))

    puff =  pd.DataFrame(phlmDict.items(), columns=['Name', 'Puff']).set_index("Name")

    trList = []
    puList = []
    krList = []
    for x in truth.index:
        trList.append( truth.loc[x].values[0] )
        krList.append( krak.loc[x].values[0] )
        flag = False
        for y,v in phlmDict.items():
            y = y.replace("_", " ")
            if x in y:
                puList.append(v)
                flag = True
                break
        if flag==False:
            puList.append(0)
    print ("Kraken v Truth")
    print (scipy.stats.spearmanr(krList, trList))
    print ("Pufferfish v Truth")
    print (scipy.stats.spearmanr(puList, trList))

@click.command()
@click.option('--sam',  help='path for the sam file generated by pufmap')
def get_stats(sam):
    #get taxa dict
    id2phlm = get_phlm_dict()

    ## Parse the BAM and perform the counting
    singCount, totCount, orphanCount, phlmDict = perform_counting(sam, id2phlm)
    #phlmDict = {}
    get_correlation(phlmDict)
    ## Calculate the stats and print it
    print_stats(singCount, totCount, orphanCount, phlmDict)
    print ("\noutput written to {}".format(os.getcwd()+"/report.txt"))

if __name__=="__main__":
    get_stats()