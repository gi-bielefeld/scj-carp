from argparse import ArgumentParser
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.font_manager import FontProperties
from scipy.stats import spearmanr,pearsonr
import numpy as np
from matplotlib.ticker import LogLocator

parser = ArgumentParser()

parser.add_argument("results")
parser.add_argument("--resultfile")
parser.add_argument("--norm-size",action="store_true")
args = parser.parse_args()


tools = dict()
plt.rcParams['font.size'] = 14

COLORS = ["lightcoral","mediumslateblue","plum"]


#plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(5))
#plt.gca().yaxis.set_minor_locator(ticker.MultipleLocator(1))

ALIAS = {
    'klebsellia10r':'K. pneumoniae',
    'plasmodium' :'P. vivax',
    'cerevisiae10r' : 'S. cerevisiae',
    'mel10r' : 'D. melanogaster',
    'arr10r' : 'A. thaliana',
    'ypestis10r' : 'Y. pestis',
    'ecoli10r' : 'E. coli',
    'carp5r' : 'C. carpio'
}

runsizes = dict()
with open(args.results) as f:
    for line in f:
        if len(line.strip())==0:
            continue
        run,tool,runsize,mr,ca = line.strip().split()
        if not tool in tools:
            tools[tool] = dict()
        if not run in tools[tool]:
            tools[tool][run] = dict()
        tools[tool][run]=(int(mr),int(ca))
        runsizes[run]=int(runsize)

tool_order = sorted([x for x in tools.keys() if not x.startswith("minigraph")])


if not args.norm_size:
    sorting = lambda r,n: -(tools[tool_order[n]][r][1])
else:
    sorting = lambda r,n: -(tools[tool_order[n]][r][1]/(tools[tool_order[n]][r][0]))/runsizes[r]
run_order = sorted(list(runsizes.keys()),key=lambda x: sorting(x,0))
#run_order.remove('carp5r')
#run_order.remove('arr10r')
#run_order.remove('mel10r')

delta_offset = 1/(len(tools)+1)
offset = 0
width=1
rects = []
#plt.grid(axis='y',zorder=0)

if not args.norm_size:
    eva = lambda tool,run: (tools[tool].get(run,(1,0))[1])
else:
    eva = lambda tool,run: (tools[tool].get(run,(1,0))[1]/tools[tool].get(run,(1,0))[0])/runsizes[run]

font = FontProperties(style='italic',size=12)

t1 = tool_order[1]
t2 = tool_order[2]
ys = [eva(t2,run) for run in run_order]
xs = [eva(t1,run) for run in run_order]
for run in run_order:
    plt.annotate(ALIAS[run],xy=(eva(t1,run),eva(t2,run)),xytext=(1, 1),textcoords="offset points",font=font)
plt.scatter(xs,ys,color='black',marker="x")

font = FontProperties(style='italic',size=12)

TOOL_ALIASES = {'cactus_perm':'Cactus','sibeliaz_ra2':'SibeliaZ (strong dup. filter)', 'sibeliaz_ra20' : 'SibeliaZ (weak dup. filter)'}

#plt.legend(rects,[TOOL_ALIASES.get(t,t) for t in tool_order])
#plt.xticks([(i+0.5*(len(tools)-1)/len(tools)) for i in range(len(run_order))],[ALIAS.get(run,run)+" ({})".format(runsizes[run]) for run in run_order],fontproperties=font)
plt.yscale('log')
plt.xscale('log')
#plt.ylabel("CARP measure")
print("Spearman: ",spearmanr(xs,ys))
print("Pearspn: ",pearsonr(xs,ys))
plt.xlabel(TOOL_ALIASES[t1])
plt.ylabel(TOOL_ALIASES[t2])
plt.xticks(np.logspace(1,5,num=5,base=10))
plt.yticks(np.logspace(1,5,num=5,base=10))

plt.show()