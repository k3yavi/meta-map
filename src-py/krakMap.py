from __future__ import print_function
import pysam
import pandas as pd
import os
import math
import click
import sys
from collections import defaultdict
import random

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
        try:
            return self.tree[tid]
        except:
            print ("ERROR: taxid not found in the tree")
            sys.exit(1)

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
            sys.exit(1)

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

        # if only one max score then move down the tree
        if all_scores.count(max_score) == 1:
            # get the taxid for max scoring child
            head = children[ all_scores.index(max_score) ]
            continue

        # if found ambiguous multiple paths then prune the tree
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

def write_taxa_count(tax_count, taxa):
    cwd = os.getcwd()
    with open(cwd+"/report.txt", 'w') as f:
        f.write("Feature\tCount\n")
        for k,v in tax_count.items():
            f.write( str(k) +"\t"+str(v)+"\n")

def get_correlation(df_list, names):
    com_all = set(df_list[0].index)
    for df in df_list:
        com_all &= set(df.index)
    print("\n\nInital Shapes: ", end="")
    for idx,df in enumerate(df_list):
        print(names[idx]+":"+str(df.shape)+"\t", end="")
    print("Number of common species: {}".format(len(com_all)))
    ct = pd.concat([ x.loc[list(com_all)] for x in df_list ], axis=1)
    print(ct.corr(method="spearman"))


def print_correlation(tax_count, taxa):
    print("Reading truth")
    with open("/mnt/scratch2/avi/meta-map/kraken/meta/truth_hc1.txt") as f:
        truth = pd.read_table(f).set_index("taxid").drop(["species",
                                                          "size",
                                                          "dataset"], 1)

    print("Reading kraken")
    with open("/mnt/scratch2/avi/meta-map/kraken/report/non_zero.txt") as f:
        krak = pd.read_table(f, header=None, sep=" ").set_index(1).drop([0])

    krDict = krak.to_dict()[0]
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
    krDict = []

    with open("krak.txt", 'w') as f:
        for k,v in new_kr.items():
            if k in truth.index:
                f.write(str(k)+"\t"+str(v)+"\n")


    print("Reading Puff")
    puff = pd.DataFrame(tax_count.items()).set_index(0)

    get_correlation([truth, krak, puff], ["truth", "kraken", "Puff"])
    get_correlation([truth, puff], ["truth", "Puff"])
    get_correlation([truth, krak], ["truth", "kraken"])
    get_correlation([krak, puff], ["kraken", "Puff"])


@click.command()
@click.option('--sam',  help='sam/sam-type pre-processed file')
@click.option('--level',  help='base level to get counts for')
def run(sam, level):
    # read in taxonomy information from the nodes.dmp file
    taxa = read_taxa(level)

    # read in reference to taxa map
    ref2tax = read_map()

    # do the counting operation
    tax_count = perform_counting(sam, ref2tax, taxa)

    # write the counted taxa
    write_taxa_count(tax_count, taxa)

    # report the correlation
    print_correlation(tax_count, taxa)

if __name__=="__main__":
    run()
