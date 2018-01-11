from __future__ import print_function
import pysam
import pandas as pd
import os
import math
import click
import sys
from collections import defaultdict
import random
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import seaborn as sns


hund = 100
hundk = 100000

class Taxa:
    def __init__(self, level):
        self.prune_level = level
        self.leaves = set([])
        self.pruning_nodes = set([])
        self.internal_nodes = set([1])
        self.num_nodes = 0
        self.tree = defaultdict(int)

    def insert(self, tid, pid, rank):
        if type(tid) == int and type(pid) == int:
            self.tree[tid] = pid
            self.num_nodes += 1
        else:
            print("ERROR: Non integer nodes for taxa tree")
            sys.exit(1)

        if rank == self.prune_level:
            self.pruning_nodes.add(tid)

        # Hackish online counting of internal and leaf nodes
        # Check if you have time
        if tid != 1 and tid != pid:
            if tid not in self.internal_nodes:
                self.leaves.add(tid)
            if pid in self.leaves:
                self.leaves.discard(pid)
            if pid not in self.internal_nodes:
                self.internal_nodes.add(pid)

    def get_parent(self, tid):
        if tid in self.tree:
            return self.tree[tid]
        else:
            return 1

    def get_path(self, tid):
        p_tid = tid
        path = [tid]
        while(p_tid != 1):
            p_tid = self.get_parent(p_tid)
            path.append( p_tid )
        return path

    def in_leaves(self, tid):
        return tid in self.leaves

# class for inverted parent -> child subtree
class Tree:
    def __init__(self, paths, intvs):
        self.mappings = defaultdict(set)
        self.coverage = defaultdict(set)
        self.leaves = []
        self.root = 1
        self.populate_tree(paths, intvs)

    def get_score(self, tid):
        try:
            return len(self.coverage[tid])
        except:
            print("ERROR: {} not found in subtree".format(tid))

    def get_coverage(self, intvs, node="int"):
        nums = set([])
        for intv in intvs:
            if node == "leaf":
                start = int(intv[0])
                end = start + int(intv[1])
                nums |= set(range(start, end))
            else:
                nums |= intv
        return nums

    def populate_tree(self, paths, intvs):
        nodes = set([])
        for idx, path in enumerate(paths):
            leaf = path[0]
            covg = self.get_coverage(intvs[idx], "leaf")

            # resolving duplicate seq 2 taxid conversion
            if leaf in self.leaves:
                if self.get_score(leaf) > len(covg):
                    continue
            else:
                self.leaves.append(leaf)
                for i in range(1, len(path)):
                    from_tid = path[i]
                    to_tid = path[i-1]
                    nodes.add(from_tid)
                    self.mappings[from_tid].add(to_tid)

            # update the new coverage
            self.coverage[leaf] = covg

        while(len(nodes) != 0):
            node = random.choice(list(nodes))
            children = self.mappings[node]
            child_intvs = []
            for child in children:
                if child not in self.coverage:
                    break
                else:
                    child_intvs.append(self.coverage[child])
            if len(child_intvs) == len(children):
                self.coverage[node] = self.get_coverage(child_intvs)
                nodes.discard(node)

def read_taxa(level):
    tf = "/mnt/scratch2/avi/meta-map/kraken/KrakenDB/taxonomy/nodes.dmp"
    # create taxa object
    taxa = Taxa(level)
    with open(tf) as f:
        for line in f:
            toks = line.rstrip("\t|\n").split("\t|\t")
            taxa.insert(int(toks[0]), int(toks[1]), toks[2])

    num_internal = float(len(taxa.internal_nodes))
    num_leaves = float(len(taxa.leaves))
    print ("Found Total # of -->")
    print ("Nodes: {}".format(taxa.num_nodes))
    print ("Internal Nodes: {0:.0f}({1:.3f}%)".format(num_internal, hund*num_internal/taxa.num_nodes))
    print ("Leaves: {0:.0f}({1:.3f}%)".format(num_leaves, hund*num_leaves/taxa.num_nodes))
    return taxa

def read_map():
    ref = "/mnt/scratch2/avi/meta-map/kraken/KrakenDB/seqid2taxid.map"
    with open(ref) as f:
        return pd.read_table(f, header=None).set_index(0).to_dict()[1]

def get_best_mapping(taxids, intvs, taxa):
    cov_threshold = 44
    n_maps = len(taxids)

    if(n_maps != len(intvs)):
        print("ERROR: number of intervals not consistent with # of refs")
        print("Exiting")
        sys.exit(1)

    #if n_maps == 1:
    #    return taxids[0]

    paths = []
    for taxid in taxids:
        paths.append( taxa.get_path(taxid) )

    # make a small tree with only relevant taxids
    sub_tree = Tree(paths, intvs)
    head = sub_tree.root

    if sub_tree.get_score(head) < cov_threshold:
        return 0

    while( head not in sub_tree.leaves ):
        children = list(sub_tree.mappings[head])
        # if only one child the just go down the tree
        if len(children) == 1:
            head = children[0]
            continue

        # get scores of all the children
        all_scores = []
        for child in children:
            all_scores.append( sub_tree.get_score(child) )

        # get the max score
        max_score = max(all_scores)

        if max_score < cov_threshold:
            break

        # if only one max score then move down the tree
        if all_scores.count(max_score) == 1:
            # get the taxid for max scoring child
            head = children[ all_scores.index(max_score) ]
            continue

        # if found ambiguous multiple paths then prune the tree
        #head = children[0]
        break

    if head not in taxa.pruning_nodes:
        # collapse this to required level
        path = taxa.get_path(head)
        for node in path:
            if node in taxa.pruning_nodes:
                head = node
                break
    return head

def perform_counting(sam, ref2tax, taxa):
    tax_count = defaultdict(int)
    read_count = 0
    not_found = []
    with open(sam) as f:
        for line in f:
            # progress monitor
            read_count += 1
            if(read_count%hundk == 0):
                print("\r {} Processed Reads".format(read_count), end="")
                sys.stdout.flush()

            rid, n_alns = line.strip().split()

            taxids = []
            intvs = []
            for _ in range(int(n_alns)):
                toks = f.next().strip().split()
                try:
                    ref_taxid = ref2tax[ toks[0] ]
                except:
                    not_found.append(toks[0])
                    continue
                n_intvs = int(toks[1])
                intv = []
                for i in range(2, 2*n_intvs+2, 2):
                    intv.append((toks[i], toks[i+1]))
                taxids.append(ref_taxid)
                intvs.append(intv)

            if len(taxids) == 0:
                continue

            key_taxid = get_best_mapping(taxids, intvs, taxa)
            tax_count[key_taxid] += 1

    print("\nTotal {} reference ids not found"
          " and missed {} alignments".format(len(set(not_found)),
                                             len(not_found)))
    return tax_count

def write_taxa_count(tax_count, tax_count_uniq, taxa, ds):
    cwd = os.getcwd()
    with open(cwd+"/"+ ds + "_report.txt", 'w') as f:
        f.write("Feature\tCount\n")
        for k,v in tax_count.items():
            f.write( str(k) +"\t"+str(v)+"\n")

    with open(cwd+"/"+ ds + "_uniq_report.txt", 'w') as f:
        f.write("Feature\tCount\n")
        for k,v in tax_count_uniq.items():
            f.write( str(k) +"\t"+str(v)+"\n")


def get_correlation(df_list, names):
    com_all = set(df_list[0].index)
    for df in df_list:
        com_all |= set(df.index)
    print("\n\nInital Shapes: ", end="")
    for idx,df in enumerate(df_list):
        print(names[idx]+":"+str(df.shape)+"\t", end="")
    print("Number of common species: {}".format(len(com_all)))
    ct = pd.concat([ x.loc[list(com_all)] for x in df_list ], axis=1).fillna(0)
    corr = ct.corr(method="spearman")["truth"]
    print("Corr: {}".format(names[1]))
    print(corr[1:].values[0])
    return ct, corr[1:].values[0]

def get_ards(ct_df, name):
    mard = np.mean( np.abs(ct_df["truth"] - ct_df[name]) / ( ct_df["truth"] + ct_df[name] ) )
    print("MARD: {}".format(name))
    print(mard)
    return mard

def print_correlation(tax_count, tax_count_unq, taxa, dataset):
    print("Reading truth")
    with open("/mnt/scratch2/avi/meta-map/kraken/meta/truth.txt") as f:
        truth = pd.read_table(f).set_index("taxid")
    truth = truth[ "Huttenhower_"+dataset == truth["dataset"] ]
    truth = truth.drop(["dataset", "species", "size"], 1)

    trDict = truth.to_dict()["counts"]
    new_kr = defaultdict(int)
    for tid,ct in trDict.items():
        if tid not in taxa.pruning_nodes:
            path = taxa.get_path(tid)
            for node in path:
                if node in taxa.pruning_nodes:
                    tid = node
                    break
        new_kr[tid] += ct

    truth = pd.DataFrame(new_kr.items()).set_index(0)
    truth.columns = ["truth"]

    print("Reading kraken")
    with open("/mnt/scratch2/avi/meta-map/kraken/krakOut/"+dataset+".rpt") as f:
        krak = pd.read_table(f, header=None).set_index(4).drop([0,1,3,5],1).drop([0])
    krDict = krak.to_dict()[2]
    new_kr = defaultdict(int)
    for tid,ct in krDict.items():
        if tid not in taxa.pruning_nodes:
            path = taxa.get_path(tid)
            for node in path:
                if node in taxa.pruning_nodes:
                    tid = node
                    break
        new_kr[tid] += ct

    krak = pd.DataFrame(new_kr.items()).set_index(0)
    krak.columns = ["kraken"]
    krDict = []
    new_kr = []

    print("Reading kraken")
    with open("/mnt/scratch2/avi/meta-map/kraken/krakOut/"+dataset+"_unfilt.rpt") as f:
        krak_unft = pd.read_table(f, header=None).set_index(4).drop([0,1,3,5],1).drop([0])
    krDict_unft = krak_unft.to_dict()[2]
    new_kr = defaultdict(int)
    for tid,ct in krDict_unft.items():
        if tid not in taxa.pruning_nodes:
            path = taxa.get_path(tid)
            for node in path:
                if node in taxa.pruning_nodes:
                    tid = node
                    break
        new_kr[tid] += ct
    krak_unft = pd.DataFrame(new_kr.items()).set_index(0)
    krak_unft.columns = ["kraken_unfiltered"]
    krDict_unft = []
    new_kr = []


    print("Reading Puff")
    puff = pd.DataFrame(tax_count.items()).set_index(0)
    puff.columns = ["MM"]

    print("Reading Puff")
    puff_unq = pd.DataFrame(tax_count_unq.items()).set_index(0)
    puff_unq.columns = ["Unique"]

    mards = []
    corrs = []
    ct_df, corr = get_correlation([truth, krak],
                                  ["truth", "kraken"])
    mard = get_ards(ct_df, "kraken")
    mards.append(mard)
    corrs.append(corr)

    ct_df, corr = get_correlation([truth, krak_unft],
                                  ["truth", "kraken-unfiltered"])
    mard = get_ards(ct_df, "kraken_unfiltered")
    mards.append(mard)
    corrs.append(corr)


    ct_df, corr = get_correlation([truth, puff],
                                  ["truth", "Puff_MM"])
    mard = get_ards(ct_df, "MM")
    mards.append(mard)
    corrs.append(corr)


    ct_df, corr = get_correlation([truth, puff_unq],
                                  ["truth", "Puff_Unq"])
    mard = get_ards(ct_df, "Unique")
    mards.append(mard)
    corrs.append(corr)

    return mards, corrs

def make_boxplot(mards, corrs, level):
    smzip = [[s,m] for sp,ma in zip(corrs, mards) for s,m in zip(sp, ma)]

    df = pd.DataFrame(smzip, columns=["Spearman", "mards"])
    df['method'] = pd.Series(["Kraken", "kraken-unfilter", "Puff_MM", "Puff_Unique"]*10)
    dd=pd.melt(df, id_vars=['method'], value_vars=["Spearman", "mards"], var_name='metric')
    sns.boxplot(x='method', y='value', data=dd, hue='metric')
    plt.show()
    plt.title(level)
    plt.savefig(level+".pdf")

@click.command()
@click.option('--level',  help='base level to get counts for')
@click.option('--report',  is_flag=True,  help='report counts or not', default=False)
@click.option('--plot', is_flag=True,  help='report counts or not')
@click.option('--ds',  help='specific dataset')
@click.option('--base',  help='new base path')
def run(level, report, ds, plot, base):
    if base == None:
        dir = "/mnt/scratch2/avi/meta-map/kraken/puff/dmps/"
    else:
        dir = base

    if ds == None:
        datasets = ["HC1", "HC2", "LC1", "LC2", "LC3", "LC4", "LC5", "LC6", "LC7", "LC8"]
    else:
        datasets = [ds]
    if level == "all":
        levels = ["species", "genus", "family", "order", "class", "phylum"]
    else:
        levels = [level]

    # read in taxonomy information from the nodes.dmp file
    taxa = read_taxa(level)

    # read in reference to taxa map
    ref2tax = read_map()

    for level in levels:
        mards = []
        corrs = []
        for ds in datasets:
            # do the counting operation
            tax_count = perform_counting(dir+ds+".dmp", ref2tax, taxa)
            tax_count_unq = perform_counting(dir+ds+"_unq.dmp", ref2tax, taxa)

            if report:
                # write the counted taxa
                write_taxa_count(tax_count, tax_count_unq, taxa, ds)

            # report the correlation
            mard, corr = print_correlation(tax_count, tax_count_unq, taxa, ds)

            mards.append(mard)
            corrs.append(corr)

        if plot:
            make_boxplot(mards, corrs, level)
        else:
            print (mards)
            print (corrs)

if __name__=="__main__":
    run()
