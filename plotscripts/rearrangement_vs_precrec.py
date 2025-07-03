import matplotlib.pyplot as plt

from argparse import ArgumentParser
import matplotlib.patches as mp

get_precision = lambda e: float(e[1])
get_recall = lambda e: float(e[2])


get_f1 = lambda e: 2*float(e[1])*float(e[2])/(float(e[1])+float(e[2]))

X = 0
YS = ['prec','recall']#,'f1']
YFUN = [get_precision,get_recall]#,get_f1]
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

colors = []
for y in YS:
    positions = []
    datapoints = []
    for x, datapointdict in data.items():
        datapoints.append(datapointdict[y])
        positions.append(float(x))
    bds = plt.violinplot(datapoints,positions,widths=0.05,showmedians=True)
    colors.append(bds['bodies'][0].get_facecolor())


plt.legend([mp.Patch(color=colors[0]),mp.Patch(color=colors[1])],['Precision','Recall'],loc=3)
plt.xticks([i/10 for i in range(1,11,3)])
plt.yticks([i/10 for i in range(4,11,2)])
plt.xlabel("Rearrangement Scale")
plt.show()