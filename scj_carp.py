from argparse import ArgumentParser, FileType
from utils import *
import sys

def canonicize_split(split,colors):
    if 2*len(split) > len(colors):
        return sorted(list(set(colors).difference(set(split))))
    elif 2*len(split)<len(colors):
        return set(split)
    elif colors[0] in split:
        return set(split)
    else:
        return set(colors).difference(set(split))

def other_split_side(split,colors):
    return set(colors).difference(split)

def set_to_id(s):
    return '...'.join(sorted(list(s)))

def id_to_set(ids):
    return set(ids.split('...'))



def color_edges(mbpg,u,b):
    return [(u,x) for x in mbpg[u] if len(mbpg[u][x]['colors'].intersection(b))>0]




def calc_carp_splits(mbpg,colors):
    splits = {}
    isolated_contested=[]
    for u,v,cls in mbpg.edges(data='colors'):
        if len(cls)==len(colors):
            continue
        sid = set_to_id(canonicize_split(cls,colors))
        clsu = [mbpg[u][x]['colors'] for x in mbpg[u] if x!=v]
        clsv = [mbpg[v][x]['colors'] for x in mbpg[v] if x!=u]
        if len(clsu)==0 and len(clsv)!=0: 
            isolated_contested.append((u,v))
            continue
        if len(clsu)!=0 and len(clsv)==0:
            isolated_contested.append((v,u))
            continue
        if len(clsu)==0 and len(clsv)==0:
            continue
        contested = False
        for clso in clsu+clsv:
            if len(clso.intersection(cls)) > 0:
                contested=True
                break

        if contested:
            continue
        if not sid in splits:
            splits[sid] = 0
        splits[sid]+=1
    #splits_stupid_case = dict([(k,0) for k in splits])
    return splits
    '''
    for k in splits:
        X=id_to_set(k)
        Y=other_split_side(k,colors)
        print(X,Y,file=sys.stderr)
        for u in mbpg.nodes():
            for A,B in [(X,Y),(Y,X)]:
                esa = color_edges(mbpg,u,A)
                if len(esa)!=1:
                    #has to be uncontested in A
                    continue
                if TELO in [x[1][-1] for x in esa[0]]:
                    continue
                esb = color_edges(mbpg,u,B)
                if len(esb) <=1 or esa[0] not in esb:
                    #has to be contested in B
                    continue
                u,v = esa[0]
                
                if len(color_edges(mbpg,v,B))==1:
                    #edge that is uncontested in A, must be contested in B and have degree 1 on other side
                    splits_stupid_case[k]+=1
                    #print("Stupid case for ",X,Y,u,v)
    for k, val in splits_stupid_case.items():
        #print("Split changed by {}".format(val))
        splits[k]+=val
    '''
    


def calc_carp_splits_alt_alt(mbpg,colors):
    splits = {}
    carpi=calc_carp_index(mbpg)
    for u,v,cls in mbpg.edges(data='colors'):
        if len(cls)==len(colors):
            continue
        sid = set_to_id(canonicize_split(cls,colors))
        clsu = [mbpg[u][x]['colors'] for x in mbpg[u] if x!=v]
        clsv = [mbpg[v][x]['colors'] for x in mbpg[v] if x!=u]
        if len(clsu)==0 and len(clsv)==0:
            continue
        contested = False
        for clso in clsu+clsv:
            if len(clso.intersection(cls)) > 0:
                contested=True
                break

        if contested:
            continue
        if not sid in splits:
            splits[sid] = 0
        splits[sid]+=1
    for sstring in splits.keys():
        A=id_to_set(sstring)
        carpia=calc_carp_index(proj_to_color(mbpg,A))
        B = other_split_side(A,colors)
        carpib=calc_carp_index(proj_to_color(mbpg,B))
        print("C(%s): %d C(%s): %d res: %d"%(','.join(A),carpia,','.join(B),carpib,carpi-carpia-carpib))


def proj_to_color(mbpg,colors):
    mbpgc = nx.Graph()
    for v,cls in mbpg.nodes(data='colors'):
        mbpgc.add_node(v,colors=cls)
    for u,v,cls in mbpg.edges(data='colors'):
        ints = cls.intersection(colors)
        if len(ints)>0:
            mbpgc.add_edge(u,v,colors=ints)
    return mbpgc

def calc_carp_splits_alt(mbpg,colors,threshold=0):
    splits = {}
    for u,v,cls in mbpg.edges(data='colors'):
        if len(cls)==len(colors):
            continue
        sid = set_to_id(canonicize_split(cls,colors))
        clsu = [mbpg[u][x]['colors'] for x in mbpg[u] if x!=v]
        clsv = [mbpg[v][x]['colors'] for x in mbpg[v] if x!=u]
        if len(clsu)==0 and len(clsv)==0:
            continue
        contested = False
        for clso in clsu+clsv:
            if len(clso.intersection(cls)) > 0:
                contested=True
                break

        if contested:
            continue
        if not sid in splits:
            splits[sid] = 0
        splits[sid]+=1
    for k,v in list(splits.items()):
        if v < threshold:
            del splits[k]
    print_splits(splits)
    print("Calculate negative contributions",file=sys.stderr)
    for sstring in list(splits.keys()):
        A=id_to_set(sstring)
        B = other_split_side(A,colors)
        for u,v in mbpg.edges():
            if len(cls)==len(colors):
                continue
            clsu = [mbpg[u][x]['colors'] for x in mbpg[u] if x!=v]
            clsv = [mbpg[v][x]['colors'] for x in mbpg[v] if x!=u]
            if len(clsu)==0 and len(clsv)==0:
                continue
            contested_A = False
            contested_B=False
            for clso in clsu+clsv:
                if len(clso.intersection(A)) > 0:
                    contested_A=True
                if len(clso.intersection(B)) > 0:
                    contested_B=True
            if contested_A and contested_B:
                splits[sstring]-=1
            #if splits[sstring]<=0:
            #    del splits[sstring]
            #    break
    return splits
                

            
            
def calc_carp_index(mbpg,get_edge_partition=False):
    if get_edge_partition:
        uncontested=[]
        contested = []
    carp_index=0
    for u,v in mbpg.edges():
        clsu = [mbpg[u][x]['colors'] for x in mbpg[u] if x!=v]
        clsv = [mbpg[v][x]['colors'] for x in mbpg[v] if x!=u]
        if len(clsu)==0 and len(clsv)==0:
            if get_edge_partition:
                uncontested.append((u,v))
            continue
        if get_edge_partition:
            contested.append((u,v))
        carp_index+=1
    return carp_index if (not get_edge_partition) else (carp_index,contested,uncontested)

    

def print_splits(splits,top=-1):
    ranked = sorted([(s,v) for s,v in splits.items()],key=lambda x:x[1],reverse=True)
    if top>0:
        if top<1:
            n = round(len(ranked)*top)
        else:
            n=top
        ranked=ranked[0:n]
    for s,v in ranked:
        print('%d\t'%v+s.replace('...','\t'))


def print_mbpg(mbpg):
    for u,v,colors in mbpg.edges(data='colors'):
        print(u,v,colors)

def filter_tree(splits):
    tree = []
    for v,k in sorted([(v,k) for k,v in splits.items()],reverse=True):
        s = id_to_set(k)
        compatible=True
        for so,_ in tree:
            if len(s.intersection(so))!=0  and not (s.issubset(so) or so.issubset(s)):
                compatible=False
                break
        if compatible:
            tree.append((k,v))
    return tree

def predict_adjacencies(tree,mg):
    pass

def partition_adjacencies(mbpg,colorset):
    ambpg = proj_to_color(mbpg,colorset)
    _,contested,uncontested = calc_carp_index(ambpg,get_edge_partition=True)
    return contested,uncontested

def main():
    parser=ArgumentParser()
    parser.add_argument('unimog')
    parser.add_argument('--alt',action='store_true',help='Activate alternative split calculation')
    parser.add_argument("--core",action='store_true')
    parser.add_argument("--greedy-tree",action='store_true')
    parser.add_argument("--candidate-adjacencies",action='store_true')
    parser.add_argument("--write-measure",type=FileType('w'))
    parser.add_argument("--write-carp-adjacencies",type=FileType('w'))
    args=parser.parse_args()
    with open(args.unimog) as f:
        gnms = readGenomes(f.readlines())
        if args.core:
            core = get_core(gnms)
            gnms = project_to_gene_set(gnms,core)
        print("Genome sizes%s:"%(" after projection to core" if args.core else ""),file=sys.stderr)
        for x,chrs in gnms:
            tlen = sum([len(chr) for _,chr in chrs])
            print(x,tlen,file=sys.stderr)
        colors = [x for x,_ in gnms]
        mbpg = build_multi_bp_graph(gnms)
        if args.write_carp_adjacencies:
            ci,_,uc= calc_carp_index(mbpg,get_edge_partition=True)

            for u,v in uc:
                un,ux = u
                vn,vx = v
                print("{un} {ux}\t{vn} {vx}".format(un=un,ux=ux,vn=vn,vx=vx),file=args.write_carp_adjacencies)
        else:
            ci=calc_carp_index(mbpg)
        print("Carp index: %d"%ci, file=args.write_measure if args.write_measure else sys.stdout)

        #print_mbpg(mbpg)
        #if args.alt:
        #    splits=calc_carp_splits_alt_alt(mbpg,colors)
        #else:
        #    splits = calc_carp_splits(mbpg,colors)
        #    print_splits(splits)

main()
