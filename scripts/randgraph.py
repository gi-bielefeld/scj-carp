import networkx as nx
import random as r
import numpy as np
import os


SAMPLES = 1


def is_head(i):
    return i%2==1


def is_tail(i):
    return i%2==0

def marker(i):
    m = (i+1)//2
    return m

def randgraph(size,conn_prob,self_conn_prob):
    rg = nx.erdos_renyi_graph(size, conn_prob)
    for v in rg.nodes():
        if r.random() < self_conn_prob:
            rg.add_edge(v,v)
        if rg.degree(v) == 0:
            rg.add_edge(v,0)
    return rg

def rand_gfa(nmarkers,conn_prob,self_conn_prob,sfr,sto,outfile):
    rg = randgraph(2*nmarkers+1,conn_prob,self_conn_prob)
    for i in range(1,nmarkers+1):
        print("S\t{}\t*\tLN:i:{}".format("S%d"%i,r.randrange(sfr,sto)),file=outfile)
    for u,v in rg.edges():
        if u==0 or v==0:
            continue
        a = marker(u)
        b = marker(v)
        ass = "S%d"%a
        bs = "S%d"%b
        orienta = "+" if is_head(u) else "-"
        orientb = "+" if is_tail(b) else "-"
        print(f"L\t{ass}\t{orienta}\t{bs}\t{orientb}\t0M",file=outfile)


def get_max_rfile_nr():
    mx = 0
    for flname in os.listdir("testfiles/random/"):
        mx = max(int(flname.removesuffix(".gfa").split("_")[-1]),mx)
    return mx


findx = get_max_rfile_nr()
for nmarkers in [25,50]:
    for conn_prob in np.arange(0.1,1.1,0.1):
        for self_conn_prob in np.arange(0.1,1.1,0.1):
            for i in range(SAMPLES):
                findx+=1
                flname = "test_{}.gfa".format(findx)
                with open(os.path.join("testfiles","random",flname),"w") as f:
                    rand_gfa(nmarkers,conn_prob,self_conn_prob,1,10,f)


