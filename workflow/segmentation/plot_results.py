from argparse import ArgumentParser
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.font_manager import FontProperties
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

tool_order = sorted([x for x in tools.keys()])


if not args.norm_size:
    sorting = lambda r,n: -(tools[tool_order[n]][r][1])
else:
    sorting = lambda r,n: -(tools[tool_order[n]][r][1]/(tools[tool_order[n]][r][0]))/runsizes[r]
run_order = sorted(list(runsizes.keys()),key=lambda x: sorting(x,0))
delta_offset = 1/(len(tools)+1)
offset = 0
width=1
rects = []
plt.grid(axis='y',zorder=0)

if not args.norm_size:
    eva = lambda tool,run: (tools[tool].get(run,(1,0))[1])
else:
    eva = lambda tool,run: (tools[tool].get(run,(1,0))[1]/tools[tool].get(run,(1,0))[0])/runsizes[run]

for i,tool in enumerate(tool_order):
    ys = [eva(tool,run) for run in run_order]
    xs = [i+offset for i in range(len(run_order))]
    boxs= plt.bar(xs,ys,width=delta_offset*width,align='edge',color=COLORS[i],zorder=2)
    rects.append(boxs)
    offset+=delta_offset
    for bar in boxs:
        height = bar.get_height()
        plt.annotate(f"{height}",xy=(bar.get_x()+bar.get_width()/2,height),fontsize=10,ha="center")
    #plt.title(tool)
    #plt.show()

font = FontProperties(style='italic',size=12)

TOOL_ALIASES = {'cactus_perm':'Cactus','sibeliaz_ra2':'SibeliaZ (strong dup. filter)', 'sibeliaz_ra20' : 'SibeliaZ (weak dup. filter)'}

plt.legend(rects,[TOOL_ALIASES.get(t,t) for t in tool_order])
plt.xticks([(i+0.5*(len(tools)-1)/len(tools)) for i in range(len(run_order))],[ALIAS.get(run,run)+" ({})".format(runsizes[run]) for run in run_order],fontproperties=font)
plt.yscale('log')
plt.ylabel("CARP measure")
plt.show()