from argparse import ArgumentParser
import re
import os
import csv

MATCHALEAF = re.compile('[(,]\\s*n[0-9]+')

def get_tree_leaves(tree):
    leaves = []
    for lraw in re.findall(MATCHALEAF,tree):
        leaves.append('n'+lraw.split('n')[1])
    return leaves


def read_leaf_genomes(gnmdir,leaves):
    genomes = []
    for lf in leaves:
        filename = os.path.join(gnmdir,lf+'_GENOME.tsv')
        with open(filename) as f:
            genes_raw = [line.split() for line in f.readlines()[1::]]
            genes = []
            for g in genes_raw:
                g[0] = int(g[0])
                genes.append(tuple(g))
            genes.sort()
            rmpos = lambda x: '' if x=='+' else x
            genes_unimog=[rmpos(o)+fam for _,fam,o,_ in genes]
            genomes.append(">%s\n%s )"%(lf,' '.join(genes_unimog)))
    return genomes
                
        
            
def main():
    parser = ArgumentParser()
    parser.add_argument('extantTree')
    parser.add_argument('genomeDir')
    args=parser.parse_args()
    with open(args.extantTree) as tf:
        tree = ''.join(tf.readlines())
    leaves = get_tree_leaves(tree)
    for gnm in read_leaf_genomes(args.genomeDir,leaves):
        print(gnm)

main()
