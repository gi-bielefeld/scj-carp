from argparse import ArgumentParser

def read_adjacencies(fname):
    adjs=[]
    with open(fname) as f:
        for line in f:
            uv = line.strip().split('\t')
            u = tuple(uv[0].split())
            v = tuple(uv[1].split())
            #dismiss telomeres
            if 'o' in u[1] or 'o' in v[1]:
                continue
            adjs.append((u,v))
    return adjs
            

def makedictl(lst):
    d=dict()
    for u,v, in lst:
        if not u in d:
            d[u] = set()
        if not v in d:
            d[v] = set()
        d[u].add(v)
        d[v].add(u)
    return d

parser=ArgumentParser()
parser.add_argument('--ground')
parser.add_argument('--result')

args = parser.parse_args()

ground = read_adjacencies(args.ground)
result = read_adjacencies(args.result)

grounddict = makedictl(ground)
resultdict = makedictl(result)

TP = [(u,v) for u,v in set(result) if v in grounddict.get(u,set())]
FP = [(u,v) for u,v in set(result) if v not in grounddict.get(u,set())]
FN = [(u,v) for u,v in set(ground) if v not in resultdict.get(u,set())]

precision=TP/(TP+FP)
recall=TP/(TP+FN)
print("%s,%s"%(precision,recall))