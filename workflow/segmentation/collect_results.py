import yaml
import os



TRIGGER_MEASURE="Carp index: "
TRIGGER_MARKER_A=" markers"
TRIGGER_MARKER_B="On "
def get_run_result(filename):
    fn = os.path.basename(filename)
    nm = fn.removesuffix(".txt")
    #toolname,markers,CARP
    with open(filename) as f:
        for line in f:
            if line.startswith(TRIGGER_MEASURE):
                ca = int(line.strip().removeprefix(TRIGGER_MEASURE))
            if line.startswith(TRIGGER_MARKER_B) and line.endswith(TRIGGER_MARKER_A):
                mr = int(line.strip().removeprefix(TRIGGER_MARKER_B).removesuffix(TRIGGER_MARKER_A))
    return nm,mr,ca



print("Run\tTool\tngenomes\tMarkers\tCARP")
with open('config.yaml') as yf:
    X=yaml.safe_load(yf)
    runs = X["runs"].keys()
    ngenomes = dict([(r,len(X["runs"][r]["genomeFiles"])) for r in X["runs"]])
    

for r in runs:
    ng = ngenomes[r]
    for fname in os.listdir(os.path.join(r,'carp')):
        
        if fname.endswith(".txt"):
            nm,mr,ca= get_run_result(os.path.join(r,'carp',fname))
            if fname.startswith("sibeliaz"):
                bn, a, s = fname.split("_")
                aparam = int(a.strip("a"))
                rel_a = str(int(aparam/ng))
                nm=bn+"_ra"+rel_a
            print(f"{r}\t{nm}\t{ng}\t{mr}\t{ca}")
