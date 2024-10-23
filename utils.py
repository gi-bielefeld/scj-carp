import networkx as nx

CHR_CIRCULAR = ')'
CHR_LINEAR = '|'
ORIENT_POSITIVE = '+'
ORIENT_NEGATIVE = '-'
EXTREMITY_TAIL = 't'
EXTREMITY_HEAD = 'h'
EXTREMITY_TELOMERE = '$'
#all telomeres are equivalent concerning the distance
TELOMERE_ID = 'TL'
HEAD = 'h'
TAIL = 't'
TELO = 'o'
TELOGENE = ''


ETYPE_MARKER = 'm'
ETYPE_ADJACENCY = 'a'



#util for dict with lists
def insertl(d, k, v):
    if k not in d:
        d[k] = []
    d[k].append(v)

def getl(d,k):
    if not k in d:
        return []
    return d[k]

def dadd(d,k,v):
    if k not in d:
        d[k] = 0
    d[k]+=v

def dgetc(d,k):
    if not k in d:
        return 0
    return d[k]


def get_core(gnms):
    core = None
    for nm, chrs in gnms:
        marker_set = set()
        for tp, markers in chrs:
            for orient,marker in markers:
                marker_set.add(marker)
        if core is None:
            core=marker_set
        else:
            core.intersection_update(marker_set)
    return core
    
def project_to_gene_set(gnms,geneset):
    gnms_=[]
    for nm, chrs in gnms:
        chrs_ = []
        for tp, markers in chrs:
            markers_=[]
            for orient, marker in markers:
                if marker in geneset:
                    markers_.append((orient,marker))
            if len(markers_)>0:
                chrs_.append((tp,markers_))
        gnms_.append((nm,chrs_))
    return gnms_



'''
This function reads a unimog (genome format) file, but does not yet assign markers unique ids!
'''
def readGenomes(data, genomesOnly=None):
    """Read genome in UniMoG format
    (https://bibiserv.cebitec.uni-bielefeld.de/dcj?id=dcj_manual)"""

    res = list()

    # helper function for parsing each individual gene
    str2gene = lambda x: x.startswith(ORIENT_NEGATIVE) and (ORIENT_NEGATIVE, \
            x[1:]) or (ORIENT_POSITIVE, x.lstrip(ORIENT_POSITIVE))
    # process each line, assuming that the file is well-formatted
    skip = False
    for line in data:
        line = line.strip()
        if line:
            if line.startswith('>'):
                genomeName = line[1:].strip()
                if genomesOnly == None or genomeName in genomesOnly:
                    skip = False
                    res.append((genomeName, list()))
                elif genomesOnly:
                    skip = True
            elif line[-1] not in (CHR_CIRCULAR, CHR_LINEAR):
                raise Exception('Invalid format, expected chromosome to ' + \
                        'end with either \'%s\' or \'%s\'' %(CHR_CIRCULAR, \
                        CHR_LINEAR))
            elif not skip:
                res[-1][1].append((line[-1], list(map(str2gene, line[:-1].split()))))
    return res


def raw_extremities(chrm):
    exts = []
    for o,g in chrm[1]:
        xs = [(g, e) for e in [TAIL, HEAD]]
        if o == ORIENT_NEGATIVE:
            xs = xs[::-1]
        exts.extend(xs)
    if chrm[0] == CHR_LINEAR:
        t1 = (TELO,TELO)
        t2 = (TELO,TELO)
        exts = [t1]+exts+[t2]
    return exts


def build_multi_bp_graph(gnms):
    mbpg = nx.Graph()
    for color, chrs in gnms:
        for chrm in chrs:
            xts = raw_extremities(chrm)
            for x in xts:
                if mbpg.has_node(x):
                    mbpg.nodes[x]['colors'].add(color)
                else:
                    mbpg.add_node(x,colors=set([color]))
            if chrm[0] == CHR_LINEAR:
                start=1
                last = xts[0]
            else:
                start=0
                last=xts[-1]
            etype=ETYPE_ADJACENCY
            for x in xts[start::]:
                if etype==ETYPE_ADJACENCY:
                    if mbpg.has_edge(last,x):
                        mbpg[last][x]['colors'].add(color)
                    else:
                        mbpg.add_edge(last,x,colors=set([color]))
                last=x
                etype=ETYPE_MARKER if etype==ETYPE_ADJACENCY else ETYPE_ADJACENCY
    return mbpg
            
