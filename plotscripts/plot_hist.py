import matplotlib.pyplot as plt
from argparse import ArgumentParser
from math import floor

parser = ArgumentParser()
parser.add_argument("tsv")
parser.add_argument("--num-buckets",type=int,default=100)
parser.add_argument("--color-ns",nargs='*',type=int)
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







plt.bar(bucket_centers,buckets,width=(mx-mn)/args.num_buckets)

if args.color_ns:
    for (lo,hi) in zip(args.color_ns,args.color_ns[1::]):
        plt.bar(bucket_centers[lo:hi],buckets[lo:hi],width=(mx-mn)/args.num_buckets)

curr_t = list(range(0,12000,500))
curr_l = curr_t
plt.xticks(curr_t+[50,140,800,10000],labels=curr_l+["A","B","C","D"])
plt.xlabel("CARP measure in Vicinity of 1000 BP")
plt.ylabel("Number of Nodes")
#plt.yscale("log")
plt.show()