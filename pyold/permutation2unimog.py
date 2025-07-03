from argparse import ArgumentParser
import sys

parser=ArgumentParser()
parser.add_argument("permutation_file")
parser.add_argument("fastafiles",nargs='+')
args = parser.parse_args()


seq_map = {}
for f in args.fastafiles:
    with open(f) as fl:
        for line in fl:
            if line.startswith('>'):
                seq_id = line[1::].split()[0]
                if seq_id in seq_map:
                    print("Warning, sequence already seen: {}".format(seq_id),file=sys.stderr)
                seq_map[seq_id]=f

gnms = dict([(g,[]) for g in seq_map.values()])

with open(args.permutation_file) as f:
    genome=None
    for line in f:
        if line.startswith('>'):
            genome=seq_map[line[1::].strip()]
        else:
            line=line.replace('+','')
            #Todo: allow circular fragments
            line=line.replace('$','|')
            gnms[genome].append(line)
    
for gnm,chrs in gnms.items():
    print(">%s"%gnm)
    for chr in chrs:
        print(chr.strip())
            