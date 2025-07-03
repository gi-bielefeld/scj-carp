# SCJ-CARP

## How to build

`cargo build --release`

The resulting binaries (`carp`,`carpscan`,`extract`) are then found in `./target/release/`.

## How to run


All binaries need to be run with either `--gfa <your-file> ` to read `<your-file>` in [gfa format](https://gfa-spec.github.io/GFA-spec/GFA1.html) or with `--unimog <your-file> ` to read `<your-file>` in [unimog format](https://bibiserv.cebitec.uni-bielefeld.de/dcj).

Currently we do not support the inferrence of telomeres from paths in gfa files. This may lead to small differences in the CARP measure.

Optionally the binaries support the following parameters:

| Parameter | Behavior |
| ------ | ------ |
| `-s`/`--size-thresh <st>`| Filter out all nodes smaller than `<st>`. Note: Since unimog files do not support node lengths, this will filter all nodes in a graph from a unimog file |
| `-t`/`--num-threads <t>`       | Use `<t>` threads for the main computation of the program. This currently does not apply to file reading or graph timming.       |
| `-h`/`--help`       | Displays a help text for the given program |

### `carp`

This program calculates the SCJ CARP measure for the given pangenome and outputs it to the command line. 

`-m`/`--write-measure <p>` writes the CARP measure to file `<p>`.

`-a`/`--write-ancestor <p>`  writes one potential set of  ancestral adjacencies to file `<p>`.

<details><summary>Example</summary>

`./target/release/carp --gfa test.gfa -m test_measure.txt -a test_ancestor.txt `

</details>

### `carpscan`

This program calculates the SCJ CARP measure for the environment of each node in the (trimmed) graph.
By default it outputs the 1% nodes with the most complex environment as identified by the SCJ CARP measure, but it can also color the graph by complexity and output a histogram of complexities.

`-c`, `--context-len <c>` Defines the context length `<c>` in base pairs that will be regarded around each node. Note that since unimog does not support node lengths, for unimog files this is instead the number of nodes in the context.

`--colored-gfa <f>`         Outputs an annotated gfa to `<f>` visualizing complexities. Can be opened in bandage.

`--output-histogram <f>`    Outputs counts for a histogram of complexities. Use `plotscripts/plot_hist.py` to visualize it.

`--lower-percentile <lo>`   Output node ids that lie between the lower and higher percentile to standard output. Default 0.99.

`--higher-percentile <hi>`  Output node ids that lie between the lower and higher percentile to standard output. Default 1.00.

<details><summary>Example</summary>

`./target/release/carpscan --gfa test.gfa  --context-len 2000 --lower-percentile 0.49 --higher-percentile 0.51 --output-histogram test.hist --colored-gfa test_colored.gfa  > test_average_nodes.txt `

View the histogram with: ` python3 plotscripts/plot_hist.py test.hist  --num-buckets 1000`

Open `test_colored.gfa` in bandage for a visualization of node complexities.
</details>

### `extract`

This program allows to ectract the surrounding graph of a node to gfa. This gfa file is output to stdout and needs to be piped into a file.

`-n, --start-node <n>`    Id `<n>` of the start node.

`-d, --max-dist <d>`    Defines the context length `<d>` in base pairs that will be regarded around each node. Note that since unimog does not support node lengths, for unimog files this is instead the number of nodes in the context.


