import utils as ut
from argparse import ArgumentParser
import networkx
parser =ArgumentParser()
parser.add_argument('genome')
args = parser.parse_args()
with open(args.genome) as f:
    gnms = ut.readGenomes(f.readlines())

graph = ut.build_multi_bp_graph(gnms)

for u,v in graph.edges():
    un,ux =  u
    vn,vx = v
    print("{un} {ux}\t{vn} {vx}".format(un=un,ux=ux,vn=vn,vx=vx))

    

