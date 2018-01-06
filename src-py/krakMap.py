from __future__ import print_function
import pysam
import pandas as pd
import os
import math
import click
import sys
from collections import defaultdict

hund = 100

class Taxa:
    def __init__(self):
        self.leaves = set([])
        self.internal_nodes = set([1])
        self.num_nodes = 0
        self.tree = defaultdict(int)

    def insert(self, tid, pid):
        if type(tid) == int and type(pid) == int:
            self.tree[tid] = pid
            self.num_nodes += 1
        else:
            print("ERROR: Non integer nodes for taxa tree")
            exit(1)

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
            exit(1)

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
    def __init__(self, paths):
        self.mappings, self.leaves = self.get_tree(paths)
        self.root = 1

    def get_tree(paths):
        mappings = defaultdict(set)
        leaves = []
        for path in paths:
            leaves.append(path[0])
            for i in range(1, len(path)):
                from_tid = path[i]
                to_tid = path[i-1]
                mappings[from_tid].add(to_tid)
        return mappings

def read_taxa():
    tf = "/mnt/scratch2/avi/meta-map/kraken/KrakenDB/taxonomy/nodes.dmp"
    # create taxa object
    taxa = Taxa()
    with open(tf) as f:
        for line in f:
            toks = line.rstrip("\t|\n").split("\t|\t")
            taxa.insert(int(toks[0]), int(toks[1]))

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

def get_mapping(taxids, intvs, taxa):
    n_maps = len(taxids)

    if(n_maps != len(intvs)):
        print("ERROR: number of intervals not consistent with # of refs")
        print("Exiting")
        exit(1)

    if n_maps == 1:
        if taxa.in_leaves(taxids[0]):
            return taxids[0]
        else:
            print("ERROR: taxid not from the leaves\n Exiting")
            exit(1)

    paths = []
    for taxid in taxids:
        if taxa.in_leaves(taxid):
            paths.append( taxa.get_path(taxid) )
        else:
            print("ERROR: taxid not from the leaves\n Exiting")
            exit(1)

    # make a small tree with only relevant taxids
    sub_tree = Tree(paths)
    head = sub_tree.root
    while( head not in sub_tree.leaves ):
        children = sub_tree.mappings[head]
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
    return head

def perform_counting(sam, ref2tax, taxa):
    tax_count = defaultdict(int)
    with open(sam) as f:
        for line in f:
            rid, n_alns = line.strip().split()
            taxids = []
            intvs = []
            for _ in range(int(n_alns)):
                toks = f.next().strip().split()
                ref_taxid = ref2tax[ toks[0] ]
                n_intvs = int(toks[1])
                intv = []
                for i in range(2, 2*n_intvs+2, 2):
                    intv.append((toks[i], toks[i+1]))
                taxids.append(ref_taxid)
                intvs.append(intv)
            key_taxid = get_mapping(taxids, intvs, taxa)
            tax_count[key_taxid] += 1
    return tax_count

@click.command()
@click.option('--sam',  help='sam/sam-type pre-processed file')
def run(sam):
    # read in taxonomy information from the nodes.dmp file
    taxa = read_taxa()

    # read in reference to taxa map
    ref2tax = read_map()

    # do the counting operation
    tax_count = perform_counting(sam, ref2tax, taxa)

    print(len(tax_count))

if __name__=="__main__":
    run()
