import matplotlib.pyplot as plt
from argparse import ArgumentParser
from matplotlib import colormaps
import numpy as np
from math import floor

parser = ArgumentParser()
parser.add_argument("tsv")
parser.add_argument("--num-buckets",type=int,default=100)
parser.add_argument("--color-ns",nargs='*',type=int)
parser.add_argument("--xlabel",default="CARP measure in vicinity of node")
args = parser.parse_args()


weights = []
elems = []
with open(args.tsv) as fl:
    for line in fl:
        entries = [int(x) for x in line.strip().split()]
        elems.append(entries[0])
        weights.append(entries[1])

#plt.hist(elems,weights=weights,bins=args.num_buckets)

buckets = [0 for _ in range(args.num_buckets)]

mn = min(elems)
mx = max(elems)

bucketnum = lambda x : int(floor((x-mn)/(mx-mn+1) * (args.num_buckets)))
bucket_centers = [((i+0.5)*mx/args.num_buckets+mn) for i in range(args.num_buckets)]

for e,w in zip(elems,weights):
    b = bucketnum(e)
    buckets[b]+=w








if args.color_ns:
    n_colors = len(args.color_ns)
else:
    n_colors = 1
cmap = plt.get_cmap('inferno')
colors = [cmap(i) for i in np.linspace(0, 1, n_colors)]


plt.bar(bucket_centers,buckets,width=(mx-mn)/args.num_buckets,color=colors[0])

if args.color_ns:
    for lo,hi,cl in zip(args.color_ns,args.color_ns[1::],colors[1::]):
        plt.bar(bucket_centers[lo:hi],buckets[lo:hi],width=(mx-mn)/args.num_buckets,color=cl)

curr_t = list(range(5000,38000,5000))
curr_l = curr_t
plt.xticks(curr_t+[100,500,1500,35000],labels=curr_l+["A","B","C","D"])
plt.xlabel(args.xlabel)
plt.ylabel("Number of Nodes")
#plt.yscale("log")
plt.show()