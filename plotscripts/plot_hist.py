import matplotlib.pyplot as plt
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("tsv")
parser.add_argument("--num-buckets",type=int,default=100)
args = parser.parse_args()


weights = []
elems = []
with open(args.tsv) as fl:
    for line in fl:
        entries = [int(x) for x in line.strip().split()]
        elems.append(entries[0])
        weights.append(entries[1])

plt.hist(elems,weights=weights,bins=args.num_buckets)
plt.show()