from argparse import ArgumentParser
import matplotlib.pyplot as plt

CARP_PREFIX="crp:i:"
COLOR_PREFIX="CL:z:"


plt.rcParams.update({'font.size': 20})



mins = dict()
maxs = dict()
parser = ArgumentParser()

parser.add_argument("gfa")

args = parser.parse_args()

with open(args.gfa) as f:
    for line in f:
        if CARP_PREFIX not in line or COLOR_PREFIX not in line:
            continue
        entries = line.split()
        carp = int([e.removeprefix(CARP_PREFIX) for e in entries if e.startswith(CARP_PREFIX)][0])
        color = [e.removeprefix(COLOR_PREFIX) for e in entries if e.startswith(COLOR_PREFIX)][0]
        mins[color]=min(mins.get(color,carp),carp)
        maxs[color]=max(maxs.get(color,carp),carp)




colors = sorted(mins.keys(),key=lambda x: -mins[x])


for color in colors:
    plt.bar(0,height=maxs[color],color=color)

plt.yscale("log")
plt.xticks([])
plt.ylabel("SCJ-CARP measure")
plt.show()