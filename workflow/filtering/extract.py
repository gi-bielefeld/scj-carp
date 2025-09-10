from argparse import ArgumentParser
import os


parser = ArgumentParser()

parser.add_argument("resultdir")

args = parser.parse_args()

measures = dict()
marker_nums = dict()
times = dict()
def tool_size_filter(fname):
    fname = fname.removesuffix(".txt")
    tool,filters =  fname.split("_s")
    return tool, int(filters)

for f in os.listdir(os.path.join(args.resultdir,"measures")):
    tool, s = tool_size_filter(f)
    measures[tool]=dict()
    marker_nums[tool] = dict()
    times[tool]=dict()
    with open(os.path.join(args.resultdir,"measures",f)) as fl:
        lines = fl.readlines()
        measure = [int(l.removeprefix("Carp index: ").strip()) for l in lines if l.startswith("Carp index: ")][0]
        nmarkers = [int(l.removeprefix("Number of markers: ").strip()) for l in lines if l.startswith("Number of markers: ")][0]
        measures[tool][s]=measure
        marker_nums[tool][s]= nmarkers
    with open(os.path.join(args.resultdir,"bench",f)) as fl:
        tool, s = tool_size_filter(f)
        with open(os.path.join(args.resultdir,"bench",f)) as fl:
            seconds = float(fl.readlines()[-1].split())
            times[tool][s]=seconds

for tool in measures:
    for s in measures[tool]:
        time = times[tool][s]
        measure = measures[tool][s]
        nmarker = marker_nums[tool][s]
        print(f"{tool}\t{s}\t{nmarker}\t{measure}\t{time}")