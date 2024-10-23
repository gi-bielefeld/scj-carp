from argparse import ArgumentParser, FileType
from utils import *
import sys

        
            
def calc_carp_index(mbpg,get_edge_partition=False):
    if get_edge_partition:
        uncontested=[]
        contested = []
    carp_index=0
    for u,v in mbpg.edges():
        if u == (TELO,TELO) or v == (TELO,TELO):
            uncontested.append((u,v))
            continue
        if u==v:
            contested.append((u,v))
            carp_index+=1
            continue
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



def print_mbpg(mbpg):
    for u,v,colors in mbpg.edges(data='colors'):
        print(u,v,colors)

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
        pgnmsz=0
        for x,chrs in gnms:
            tlen = sum([len(chr) for _,chr in chrs])
            pgnmsz+=tlen
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
        cin=ci/pgnmsz
        print("Carp index: %d,%d,%f"%(ci,pgnmsz,cin), file=args.write_measure if args.write_measure else sys.stdout)

if __name__ == '__main__':
    main()
