import matplotlib.pyplot as plt

from argparse import ArgumentParser
import matplotlib.patches as mp

get_measure = lambda e: float(e[3])



X = 0
YS = ['measure']
YFUN = [get_measure]

plt.rcParams.update({'font.size': 14})




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


for i in range(3):
    plt.plot(-0.1,-1)

colors = []
for y in YS:
    positions = []
    datapoints = []
    for x, datapointdict in data.items():
        datapoints.append(datapointdict[y])
        positions.append(float(x))
    bds = plt.violinplot(datapoints,positions,widths=0.05,showmedians=True)
    colors.append(bds['bodies'][0].get_facecolor())


#plt.legend([mp.Patch(color=colors[0])],['Precision'],loc=3)
plt.xticks([i/10 for i in range(1,11)])
plt.xlabel("Rearrangement Scale")
plt.ylabel("SCJ-CARP Measure")
plt.show()