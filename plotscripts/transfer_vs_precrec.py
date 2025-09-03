import matplotlib.pyplot as plt

from argparse import ArgumentParser
import matplotlib.patches as mp
import math
get_measure = lambda e: float(e[1])

get_prec = lambda e: float(e[2])
get_rec = lambda e: float(e[3])


X = 0
YS = ['prec','recall']
YFUN = [get_prec,get_rec]
plt.rcParams.update({'font.size': 20})




parser = ArgumentParser()

parser.add_argument("csv")

args = parser.parse_args()


data = {}
with open(args.csv) as f:
    for line in f:
        entries = line.strip().split(',')
        x = entries[X]
        if not x in data:
            data[x]={}
            for y in YS:
                data[x][y]=[]
        for i,y in enumerate(YS):
            data[x][y].append(YFUN[i](entries))


#for i in range(3):
#    plt.plot(-0.1,-1)

colors = []
for y in YS:
    positions = []
    datapoints = []
    for x, datapointdict in data.items():
        datapoints.append(datapointdict[y])
        positions.append(math.log(float(x)))
    bds = plt.violinplot(datapoints,positions,widths=0.3,showmedians=True)
    colors.append(bds['bodies'][0].get_facecolor())

xs = sorted([x for x in data])

plt.legend([mp.Patch(color=colors[0]),mp.Patch(color=colors[1])],['Precision','Recall'],loc=3)

#plt.legend([mp.Patch(color=colors[0])],['Precision'],loc=3)
tikps = sorted(xs,key=lambda x: float(x))[::2]
plt.xticks([math.log(float(x)) for x in tikps],[str(x) for x in tikps])
plt.yticks([i/10 for i in range(4,11,2)])
plt.xlabel("Transfer Rate")
#plt.xscale('log')
#plt.ylabel("SCJ-CARP Measure")
plt.show()
