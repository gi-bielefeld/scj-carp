import yaml
import os
import json
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

WORKDIR = "config2tableworkdir"
if not os.path.isdir(WORKDIR):
    os.mkdir(WORKDIR)
with open('config.yaml') as yf:
    X=yaml.safe_load(yf)
    runs = sorted(X["runs"].keys())
    genomefiles = dict([(r,X["runs"][r]["genomeFiles"]) for r in X["runs"]])
    ks = dict([(r,X["runs"][r]["k_sibeliaz"]) for r in X["runs"]])
    references = dict([(r,X["runs"][r].get("reference",min(X["runs"][r]["genomeFiles"].keys()))) for r in X["runs"]])


print("\\begin{longtable}{llrlrrr}")
print("\\toprule")
print("Taxon&\\texttt{Accession.Version}&\#Genomes & Assembly lvl. &$k$ (SbZ) & $a$ (SbZs) & $a$ (SbZw)\\\\")

def get_accession_version(flname):
    return "_".join(os.path.basename(flname).split("_")[0:2])


def sanitize_underscore(s):
    return "\\_".join(s.split("_"))

def get_assembly_level(run,files_acc_version):
    filename = os.path.join(WORKDIR,run+".txt")
    os.system("datasets summary genome accession %s > %s"%(" ".join(files_acc_version),filename))
    with open(filename) as f:
        jsondata = json.load(f)
    levels = set([ r["assembly_info"]["assembly_level"] for r in jsondata['reports']])
    return levels

def pretty_level(levels):
    if levels == {'Chromosome','Complete Genome'}:
        return '\\texttt{chromosome+}'
    elif levels == {'Chromosome'}:
        return '\\texttt{chromosome}'
    elif levels == {'Complete Genome'}:
        return '\\texttt{complete}'
    else:
        return ", ".join(levels)

for run in runs:
    print("\\midrule")
    reference = "\\texttt{%s}${}^\star$"%get_accession_version(genomefiles[run][references[run]])
    a1 = len(genomefiles[run])*2
    a2 = len(genomefiles[run])*2*10
    files_acc_version = [get_accession_version(fl) for x,fl in genomefiles[run].items() if x != references[run]]
    ass_level = get_assembly_level(run,files_acc_version)
    files_acc_version = [sanitize_underscore(f) for f in files_acc_version]
    print("\\textit{%s}&%s&%i&%s&%i&%i&%i\\\\"%(ALIAS.get(run,run),sanitize_underscore(reference),len(genomefiles[run]),pretty_level(ass_level),ks[run],a1,a2))
    for fl in files_acc_version:
        print("&\\texttt{%s}\\\\"%fl)

print("\\bottomrule")
print("\\end{longtable}")