# SCJ-CARP

SCJ-CARP is a prototype solving CARP, an ancestral reconstruction problem for pangenomes, under the SCJ-model.
SCJ-CARP allows to calculate the SCJ-CARP measure, an index of pangenome structural complexity.

A necessary preprocessing step for the current version is some form of collinear block detection (for example SibeliaZ) and translating the block orders for each genome into unimog format.

## How to run SCJ-CARP

`python3 scj_carp.py <unimog>` where `<unimog>` contains all genomes in unimog format (https://bibiserv.cebitec.uni-bielefeld.de/dcj).

### Optional Parameters



| parameter  | purpose |
| ------ | ------ |
| `--write-measure <file>` | Write the SCJ-CARP measure to `<file>` |
| `--write-carp-adjacencies <file2>` | Write one possible set of adjacencies of an SCJ-CARP ancestor to `<file2>` |
| `--core` | Project the set of markers to only those occurring in all extant genomes. (Currently not recommended) |



## Snakemake Workflow

To replicate the experiments from the manuscript, first change the variable `ZOMBI_HOME` in the snakefile to the installation of ZOMBI on your machine.
Then navigate in your command line to the top directory of this repository.

The two experiments can be run with the following rules:
| Experiment | Rule | Result File | Plot Scripts (`./plotscripts/`) |
| ------| ------ | ----- | ---- |
| Scaling rates of Adjacency modifying operations (Duplications, Losses, Inversions, Transpositions, Originations) | `all_treescale_results` | `treescale/treescale.txt` | `rearrangement_vs_measure.py`, `rearrangement_vs_precrec.py`|
| Scaling rates of Non-Adjacency modifying operations (Horizontal (Replacement) Transfers) | `all_transfer` |`transferexp/transfercarp.txt` | `transfer_vs_measure.py`, `transfer_vs_precrec.py` |

To run the workflow, simply type `snakemake <rule> --cores <n-cores>`, where `<rule>` is one of the above rules and `<n-cores>` is the number of cores you want to use.

For example, `snakemake all_transfer --cores 4`.

The table of results (SCJ-CARP measure) and precision/recall can then be found in the corresponding result file.
To replicate one of the pots, simply run the corresponding script, inputting the result file as an argument, i.e.

`python3 plotscripts/rearrangement_vs_precrec.py treescale/treescale.txt`.
