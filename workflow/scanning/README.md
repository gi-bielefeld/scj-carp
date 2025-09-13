## Introduction

This workflow evaluates the CARP measure for subgraphs of a certain basepair length around each node, as described in Section "Scoring Region Complexities in a Bacillus Subtilis de Bruijn Graph" of the article manuscript.

## Configuration

To use this workflow, you need to adapt the `config.yaml` file to your needs. The config file should contain the following parameters:

* `input_graph`: input graph file
* `output_dir`: output directory for the results
* `contexts`: list of context lengths to evaluate

## Workflow Overview


The target rule `all_histogram` generates the histogram files for all specified context lengths. The results can be found in the `hist` directory within the specified `output_dir`.

## Post-processing

The resulting histograms can be plotted using the `plotscripts/plot_hist.py` script located at the root of the repository.
