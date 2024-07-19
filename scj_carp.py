from argparse import ArgumentParser
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
    return splits


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
                

            
            
def calc_carp_index(mbpg):
    carp_index=0
    for u,v in mbpg.edges():
        clsu = [mbpg[u][x]['colors'] for x in mbpg[u] if x!=v]
        clsv = [mbpg[v][x]['colors'] for x in mbpg[v] if x!=u]
        if len(clsu)==0 and len(clsv)==0:
            continue
        carp_index+=1
    return carp_index

    

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
def main():
    parser=ArgumentParser()
    parser.add_argument('unimog')
    parser.add_argument('--alt',action='store_true',help='Activate alternative split calculation')
    args=parser.parse_args()
    with open(args.unimog) as f:
        gnms = readGenomes(f.readlines())
        core = get_core(gnms)
        gnms = project_to_gene_set(gnms,core)
        print("Genome sizes after projection to core:",file=sys.stderr)
        for x,chrs in gnms:
            tlen = sum([len(chr) for _,chr in chrs])
            print(x,tlen,file=sys.stderr)
        colors = [x for x,_ in gnms]
        mbpg = build_multi_bp_graph(gnms)
        print("Carp index: %d"%calc_carp_index(mbpg))
        #print_mbpg(mbpg)
        if args.alt:
            splits=calc_carp_splits_alt_alt(mbpg,colors)
        else:
            splits = calc_carp_splits(mbpg,colors)
            print_splits(splits)

main()
