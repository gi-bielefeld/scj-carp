import sys

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("--name")
parser.add_argument("--files",nargs='*',required=True)
args = parser.parse_args()
if args.name:
	run_name=args.name
else:
	x = list(set(["/".join(x.split("/")[-2:-1]) for x in args.files]))
	if len(x) == 1:
		run_name=x[0]
	else:
		run_name='???? #TODO:change'

print(f" {run_name}:")
print("  k_sibeliaz: 25 #TODO:change")
print("  genomeFiles:")
for x in args.files:
	nm = x.split("/")[-1].split(".")[0]
	print("   ",nm,":","'"+x+"'")
