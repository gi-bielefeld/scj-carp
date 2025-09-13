
## Introduction


This workflow evaluates the CARP measure for different filterings of a graph, as described in the "Comparing Human Pangenome Graphs using CARP" section of the article manuscript.

## Configuration


To use this workflow, you need to adapt the `config.yaml` file to your needs. The config file should contain the following parameters:

* `input_dir`: input directory containing the graph files
* `input_graphs`: list of graph files to process
* `output_dir`: output directory for the results
* `size_thresholds`: list of size thresholds to apply

## Workflow Overview


The target rule `all_measure` generates the measure files for all genomes and size thresholds. The results can be found in the `measures` directory within the specified `output_dir`.

## Post-processing


The results can be extracted and processed using the `extract.py` script. Additionally, the results can be visualized using the `plot_results.py` script.
