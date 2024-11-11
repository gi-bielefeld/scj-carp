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



